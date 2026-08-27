function results = fcst_SVo_oneorigin(thisT, Tdata, data, dates_trans, ...
    N, p, MCMCdraws, burnin, MCMCreps, Ndraws, fcstNdraws, fcstNhorizons, ...
    theta, mnPrior, doRobustPrior, check_stationarity, setQuantiles, Nlogtwopi)
% FCST_SVO_ONEORIGIN  SVo BVAR estimation + forecast at one origin
%
%   SVo (SV + common scalar outlier) from CCMM (2022, REStat) supplement:
%     v_t = obar_t * A^{-1} Lambda_t^{0.5} e_t
%     obar_t: SCALAR common outlier in {1,...,maxscale}
%     P(obar_t = 1) = 1 - pi, P(obar_t = k>1) = pi/(maxscale-1)
%     h_t = h_{t-1} + eta_t, eta_t ~ N(0,Phi)
%
%   MCMC steps (follows mcmcVARSVobar.m ordering):
%     1. PAI via CTA (sqrtht includes obar)
%     2. A via equation-by-equation (sqrtht includes obar)
%     3. obar via Gibbs (resid standardized by pure SV)
%     4. h via KSC (resid standardized by obar)
%     5. Phi via IW

%% ---- Prepare data ----
thisdata = data(1:thisT, :);
Nobs = thisT;

lags_tmp = zeros(Nobs, N*p);
for l = 1:p
    lags_tmp(p+1:Nobs, (N*(l-1)+1):N*l) = thisdata(p+1-l:Nobs-l, 1:N);
end
X = [ones(Nobs-p,1), lags_tmp(p+1:Nobs,:)];
Y = thisdata(p+1:end,:);
[T_est, K] = size(X);
Klagreg    = K - 1;
ndxKlagreg = 1 + (1:Klagreg);

Xjumpoff    = zeros(K, 1);
Xjumpoff(1) = 1;
for l = 1:p
    Xjumpoff(1+(l-1)*N+(1:N)) = thisdata(Nobs-(l-1), 1:N);
end

yrealized = NaN(N, fcstNhorizons);
for hh = 1:fcstNhorizons
    if thisT + hh <= Tdata
        yrealized(:,hh) = data(thisT+hh, :)';
    end
end

%% ---- Minnesota prior (sparse) ----
np = 12;
if doRobustPrior
    trainingT = sum(dates_trans(p+1:thisT) < datenum(2020,3,1));
    trainingT = max(trainingT, N+20);
else
    trainingT = T_est;
end

ARresid = NaN(trainingT-1, N);
ARrho   = NaN(2, N);
for i = 1:N
    yt_0 = [ones(trainingT-1,1), Y(1:trainingT-1,i)];
    yt_1 = Y(2:trainingT, i);
    ARrho(:,i)    = (yt_0'*yt_0)\(yt_0'*yt_1);
    ARresid(:,i)  = yt_1 - yt_0*ARrho(:,i);
end
AR_s2_vec = diag(ARresid'*ARresid) ./ (trainingT-2);
AR_s2_vec(AR_s2_vec <= 0) = 1e-6;

ARresid_full = NaN(T_est-1, N);
for i = 1:N
    yt_0 = [ones(T_est-1,1), Y(1:T_est-1,i)];
    yt_1 = Y(2:T_est, i);
    ARresid_full(:,i) = yt_1 - yt_0 * ARrho(:,i);
end

Pi_pm      = zeros(N*Klagreg, 1);
Pi_pv_diag = ones(N*Klagreg, 1);
co = 0;
sigma_const = NaN(1,N);
for i = 1:N
    sigma_const(i) = AR_s2_vec(i) * theta(3);
    for l = 1:p
        for j = 1:N
            co = co + 1;
            if (i==j)
                if l==1, Pi_pm(co) = mnPrior(i); end
                Pi_pv_diag(co) = theta(1)/(l^theta(4));
            else
                Pi_pv_diag(co) = (AR_s2_vec(i)/AR_s2_vec(j))*theta(1)*theta(2)/(l^theta(4));
            end
        end
    end
end

omega_diag = [sigma_const; reshape(Pi_pv_diag, Klagreg, N)];
omega_diag = omega_diag(:);
iV         = sparse(1:K*N, 1:K*N, 1./omega_diag);
MU_pai     = [zeros(1,N); reshape(Pi_pm,Klagreg,N)];
iVb_prior  = iV * vec(MU_pai);

%% ---- SV priors ----
d_PHI = N + 3;
s_PHI = d_PHI * (0.15*eye(N)) * 12/np;

Vol_0mean    = zeros(N,1);
Vol_0vcvsqrt = 10*speye(N);

%% ---- Obar prior: ~1 outlier every 4 years ----
SVobarmaxscale = 20;   % max integer scale (CCMM default)
SVobaralpha    = 1/(4*np) * 10*np;   % ~2.5
SVobarbeta     = 10*np - SVobaralpha; % ~117.5

% Obar state grid
obarStates.Ngrid          = SVobarmaxscale - 1;
obarStates.uniformkernel  = ones(1, SVobarmaxscale-1) / (SVobarmaxscale-1);
obarStates.values         = 1:SVobarmaxscale;
obarStates.squaredvalues  = (1:SVobarmaxscale).^2;
obarStates.log2values     = 2*log(1:SVobarmaxscale);

EYEn = eye(N);
comp = [eye(N*(p-1)), zeros(N*(p-1),N)];

%% ---- KSC mixture values ----
[gridKSC, gridKSCt, logy2offset] = getKSC7values(T_est, N);

%% ---- MCMC sampler ----
rndStream = RandStream('threefry', 'Seed', thisT);

all_ydraws   = NaN(N, fcstNhorizons, fcstNdraws);
joint_lden   = NaN(fcstNdraws, 1);
draw_counter = 0;

% Initialize
A_         = EYEn;
invA_      = EYEn;
warnstate  = warning('off', 'MATLAB:rankDeficientMatrix');
PAI        = X(1:min(trainingT,T_est),:) \ Y(1:min(trainingT,T_est),:);
warning(warnstate);

sqrtht_init = sqrt([ARresid_full(1,:).^2; ARresid_full.^2]);
sqrtht_init = max(sqrtht_init, 1e-6);
sqrtht      = sqrtht_init;            % T x N (includes obar, initially obar=1)
Vol_states  = 2*log(sqrtht_init)';    % N x T (pure log-vol)
sqrtPHI_    = sqrt(0.0001)*eye(N);
SVobarprob  = 0.1;                    % scalar probability
SVobarscale = ones(1, T_est);         % 1 x T (scalar per time period)

m = 0;
while m < MCMCreps
    m = m + 1;

    % Step 1: Draw PAI via CTA
    stationary = 0;
    while stationary == 0
        PAI = CTA(Y, X, N, K, T_est, A_, sqrtht, iV, iVb_prior, PAI, rndStream);
        if (check_stationarity==0 || max(abs(eig([PAI(ndxKlagreg,:)'; comp]))) < 1)
            stationary = 1;
        end
    end
    RESID = Y - X*PAI;

    % Step 2: Draw A (diffuse prior = OLS)
    for ii = 2:N
        y_spread_adj = RESID(:,ii) ./ sqrtht(:,ii);
        X_spread_adj = RESID(:,1:ii-1) ./ sqrtht(:,ii);
        ZZ = X_spread_adj'*X_spread_adj;
        Zz = X_spread_adj'*y_spread_adj;
        try
            Valpha_post = ZZ \ eye(ii-1);
            alpha_post  = Valpha_post * Zz;
            alphadraw   = alpha_post + chol(Valpha_post,'lower')*randn(rndStream,ii-1,1);
            A_(ii,1:ii-1) = -alphadraw';
        catch
        end
    end
    invA_ = A_ \ EYEn;

    % Step 3: Draw obar (scalar common outlier)
    % Structural residuals standardized by pure SV (removing obar from sqrtht)
    yresid = RESID * A_';                                  % T x N
    zresid = yresid ./ (sqrtht ./ SVobarscale');           % T x N (pure SV standardized)

    % Inline obar Gibbs draw (avoids checkdiff dependency)
    obarPrior = [1 - SVobarprob, SVobarprob * obarStates.uniformkernel];
    ssr_obar  = sum(zresid.^2, 2) ./ obarStates.squaredvalues;  % T x maxscale
    pdfKernel = exp(-0.5 * ssr_obar) ./ (obarStates.values.^N) .* obarPrior;
    cdf_obar  = cumsum(pdfKernel, 2);
    cdf_obar(:,1:end-1) = cdf_obar(:,1:end-1) ./ cdf_obar(:,end);
    cdf_obar(:,end) = 1;
    ndx_obar  = sum(rand(rndStream, T_est, 1) > cdf_obar, 2) + 1;
    SVobarscale = obarStates.values(ndx_obar);
    SVobarscale = SVobarscale(:)';   % 1 x T

    % Update obar probability
    Noutlier   = sum(ndx_obar > 1);
    SVobarprob = betarnd(SVobaralpha + Noutlier, SVobarbeta + (T_est - Noutlier));

    % Step 4: Draw h (KSC, conditioned on current obar)
    % Remove obar from structural residuals before SV sampler
    ytilderesid = yresid' ./ SVobarscale;    % N x T (structural resid / obar)
    logy2       = log(ytilderesid.^2 + logy2offset);

    try
        [Vol_states, ~, eta] = StochVolKSCcorrsqrt(logy2, Vol_states, sqrtPHI_, ...
            Vol_0mean, Vol_0vcvsqrt, gridKSC, gridKSCt, N, T_est, rndStream);
    catch
        eta = diff([zeros(N,1), Vol_states], 1, 2);
    end

    % Recombine SV + obar
    sqrtht = (exp(Vol_states/2) .* SVobarscale)';  % T x N

    % Step 5: Draw Phi
    Zdraw       = randn(rndStream, N, T_est + d_PHI);
    sqrtPHIpost = chol(s_PHI + eta*eta', 'lower');
    sqrtZZ      = chol(Zdraw * Zdraw');
    sqrtPHI_    = sqrtPHIpost / sqrtZZ;

    %% Post burn-in: forecast
    if m > burnin

        h_last  = Vol_states(:, end);    % N x 1 pure log-vol
        mu_fcst = PAI' * Xjumpoff;       % N x 1

        % Build obar forecast CDF (scalar)
        obarPdf = [1 - SVobarprob, repmat(SVobarprob / obarStates.Ngrid, 1, obarStates.Ngrid)];
        obarCdf = cumsum(obarPdf);

        for nn = 1:Ndraws
            draw_counter = draw_counter + 1;

            fcstX = Xjumpoff;
            h_cur = h_last;

            for hh = 1:fcstNhorizons
                % Propagate SV
                h_cur = h_cur + sqrtPHI_ * randn(rndStream, N, 1);

                % Draw SCALAR common outlier
                u_obar = rand(rndStream);
                obar_draw = obarStates.values(sum(u_obar > obarCdf) + 1);

                % Total scale: obar * exp(h/2)
                total_vol = obar_draw * exp(h_cur/2);  % N x 1

                % Draw forecast
                EY    = PAI' * fcstX;
                ydraw = EY + invA_ * (total_vol .* randn(rndStream, N, 1));

                all_ydraws(:, hh, draw_counter) = ydraw;

                % Joint log score at h=1
                if hh == 1 && ~any(isnan(yrealized(:,1)))
                    dev       = yrealized(:,1) - mu_fcst;
                    e_struct  = A_ * dev;   % structural residuals
                    % |Sigma| = obar^{2N} * exp(sum(h))
                    logdetSig = N * obarStates.log2values(obar_draw) + sum(h_cur);
                    % Sigma^{-1} quadratic = sum(e^2 * exp(-h)) / obar^2
                    quad = sum(e_struct.^2 .* exp(-h_cur)) / obar_draw^2;
                    joint_lden(draw_counter) = -0.5 * (Nlogtwopi + logdetSig + quad);
                end

                % Update state
                fcstX = [1; ydraw; fcstX(2:end-N)];
            end
        end
    end
end

%% ---- Compute forecast summaries ----
results.yrealized  = yrealized;
results.yhat       = mean(all_ydraws, 3);
results.ymedian    = median(all_ydraws, 3);
results.yquantiles = prctile(all_ydraws, setQuantiles, 3);

results.crps     = NaN(N, fcstNhorizons);
results.pit      = NaN(N, fcstNhorizons);
results.logscore = NaN(N, fcstNhorizons);

for hh = 1:fcstNhorizons
    for j = 1:N
        if ~isnan(yrealized(j,hh))
            draws_jh = squeeze(all_ydraws(j, hh, :));
            results.crps(j,hh)  = crpsDraws(yrealized(j,hh), draws_jh);
            results.pit(j,hh)   = mean(draws_jh <= yrealized(j,hh));
            try
                results.logscore(j,hh) = log(ksdensity(draws_jh, yrealized(j,hh)));
            catch
                results.logscore(j,hh) = NaN;
            end
        end
    end
end

valid_lden = joint_lden(~isnan(joint_lden));
if ~isempty(valid_lden)
    maxlden = max(valid_lden);
    results.joint_logscore = log(mean(exp(valid_lden - maxlden))) + maxlden;
else
    results.joint_logscore = NaN;
end

results.ydraws_h1 = squeeze(all_ydraws(:, 1, 1:10:end));

end
