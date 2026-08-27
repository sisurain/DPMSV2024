function results = fcst_SVO_oneorigin(thisT, Tdata, data, dates_trans, ...
    N, p, MCMCdraws, burnin, MCMCreps, Ndraws, fcstNdraws, fcstNhorizons, ...
    theta, mnPrior, doRobustPrior, check_stationarity, setQuantiles, Nlogtwopi)
% FCST_SVO_ONEORIGIN  SVO BVAR estimation + forecast at one origin
%
%   SV + Stock-Watson outlier model from CCMM (2022, REStat):
%     v_t = A^{-1} O_t Lambda_t^{0.5} e_t
%     o_{j,t} in {1,...,10}, P(o>1) = pi_j ~ Beta(alpha,beta)
%     h_t = h_{t-1} + eta_t, eta_t ~ N(0,Phi)

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

% Jumpoff state
Xjumpoff    = zeros(K, 1);
Xjumpoff(1) = 1;
for l = 1:p
    Xjumpoff(1+(l-1)*N+(1:N)) = thisdata(Nobs-(l-1), 1:N);
end

% Realized future values
yrealized = NaN(N, fcstNhorizons);
for h = 1:fcstNhorizons
    if thisT + h <= Tdata
        yrealized(:,h) = data(thisT+h, :)';
    end
end

%% ---- Minnesota prior (sparse for N=42) ----
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

% Full-sample AR residuals for SV initialization
ARresid_full = NaN(T_est-1, N);
for i = 1:N
    yt_0 = [ones(T_est-1,1), Y(1:T_est-1,i)];
    yt_1 = Y(2:T_est, i);
    ARresid_full(:,i) = yt_1 - yt_0 * ARrho(:,i);
end

% Minnesota prior (sparse)
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

%% ---- Outlier prior: ~1 outlier every 4 years ----
SVOpriorobs = 10 * np;                           % 120 months
SVOalpha    = 1/(4*np) * SVOpriorobs;             % ~2.5
SVOalpha    = repmat(SVOalpha, N, 1);
SVObeta     = SVOpriorobs - SVOalpha;
SVOmaxscale = 10;

SVOstates.Ngrid      = SVOmaxscale - 1;           % 9 outlier states
SVOstates.values     = 1:SVOmaxscale;             % {1,2,...,10}
SVOstates.log2values = 2*log(SVOstates.values);

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
sqrtht     = sqrtht_init;         % T x N (includes outlier, initially = pure SV)
Vol_states = 2*log(sqrtht_init)'; % N x T (pure log-vol, no outlier)
sqrtPHI_   = sqrt(0.0001)*eye(N);
SVOprob    = repmat(0.1, N, 1);   % N x 1
SVOlog2    = zeros(N, T_est);     % N x T

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
            % Keep previous A row
        end
    end
    invA_ = A_ \ EYEn;

    % Step 3: Draw h AND outlier states jointly
    logy2 = log((RESID*A_').^2 + logy2offset);
    try
        [Vol_states, ~, eta, SV_combined, SVOlog2, SVOprob, SVOscale] = ...
            StochVolOutlierKSCcorrsqrt(logy2', Vol_states, sqrtPHI_, ...
            Vol_0mean, Vol_0vcvsqrt, SVOlog2, SVOprob, SVOalpha, SVObeta, ...
            SVOstates, gridKSC, gridKSCt, N, T_est, rndStream);
        sqrtht = SV_combined';  % T x N (includes outlier scale)
    catch
        eta = diff([zeros(N,1), Vol_states], 1, 2);
    end

    % Step 4: Draw Phi
    Zdraw       = randn(rndStream, N, T_est + d_PHI);
    sqrtPHIpost = chol(s_PHI + eta*eta', 'lower');
    sqrtZZ      = chol(Zdraw * Zdraw');
    sqrtPHI_    = sqrtPHIpost / sqrtZZ;

    %% Post burn-in: forecast
    if m > burnin
        thisdraw = m - burnin;

        h_last  = Vol_states(:, end);   % N x 1 pure log-vol (no outlier)
        mu_fcst = PAI' * Xjumpoff;      % N x 1

        for nn = 1:Ndraws
            draw_counter = draw_counter + 1;

            fcstX = Xjumpoff;
            h_cur = h_last;

            for hh = 1:fcstNhorizons
                % Propagate SV
                h_cur = h_cur + sqrtPHI_ * randn(rndStream, N, 1);

                % Draw outlier states
                u = rand(rndStream, N, 1);
                o_draw = ones(N, 1);   % default: no outlier
                isOutlier = u < SVOprob;
                if any(isOutlier)
                    % Uniform on {2,...,SVOmaxscale}
                    nout = sum(isOutlier);
                    o_draw(isOutlier) = floor(rand(rndStream, nout, 1) * (SVOmaxscale-1)) + 2;
                end

                % Total scale: sqrt(exp(h)) * o
                total_vol = exp(h_cur/2) .* o_draw;

                % Draw y
                EY    = PAI' * fcstX;
                ydraw = EY + invA_ * (total_vol .* randn(rndStream, N, 1));

                all_ydraws(:, hh, draw_counter) = ydraw;

                % Joint log score at h=1
                if hh == 1 && ~any(isnan(yrealized(:,1)))
                    dev       = yrealized(:,1) - mu_fcst;
                    e_struct  = A_ * dev;
                    total_var = exp(h_cur) .* (o_draw.^2);
                    logdetSig = sum(log(total_var));
                    joint_lden(draw_counter) = -0.5 * (Nlogtwopi + logdetSig + ...
                        sum(e_struct.^2 ./ total_var));
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

% Joint log score (log-sum-exp)
valid_lden = joint_lden(~isnan(joint_lden));
if ~isempty(valid_lden)
    maxlden = max(valid_lden);
    results.joint_logscore = log(mean(exp(valid_lden - maxlden))) + maxlden;
else
    results.joint_logscore = NaN;
end

% Subsampled h=1 draws
results.ydraws_h1 = squeeze(all_ydraws(:, 1, 1:10:end));

end
