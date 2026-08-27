function results = fcst_SVt_oneorigin(thisT, Tdata, data, dates_trans, ...
    N, p, MCMCdraws, burnin, MCMCreps, Ndraws, fcstNdraws, fcstNhorizons, ...
    theta, mnPrior, doRobustPrior, check_stationarity, setQuantiles, Nlogtwopi)
% FCST_SVT_ONEORIGIN  BVAR-t-CSV estimation + forecast at one origin
%
%   Chan (2020, JBES) model:
%     eps_t ~ N(0, lambda_t * exp(h_t) * Sigma)
%     h_t = rho * h_{t-1} + eta_t,  eta_t ~ N(0, sigh2)   [scalar AR(1)]
%     lambda_t ~ IG(nu/2, nu/2)                              [scalar t-scaling]
%     Sigma: full N×N covariance (Kronecker structure)
%
%   Uses Normal-Inverse-Wishart conjugate sampling (NO CTA needed).
%   Much faster per iteration than Cholesky SV models.

%% ---- Prepare data ----
thisdata = data(1:thisT, :);
Nobs = thisT;
tmpY = thisdata;  % (p+T_est) x N

X_lags = zeros(Nobs-p, N*p);
for lag = 1:p
    X_lags(:, (lag-1)*N+1 : lag*N) = tmpY(p-lag+1 : end-lag, :);
end
X = [ones(Nobs-p, 1), X_lags];   % T_est x K
shortY = tmpY(p+1:end, :);       % T_est x N
[T_est, K] = size(X);

% Jumpoff state for forecasting
Xjumpoff    = zeros(K, 1);
Xjumpoff(1) = 1;
for l = 1:p
    Xjumpoff(1+(l-1)*N+(1:N)) = thisdata(Nobs-(l-1), 1:N);
end

% Realized future values
yrealized = NaN(N, fcstNhorizons);
for hh = 1:fcstNhorizons
    if thisT + hh <= Tdata
        yrealized(:,hh) = data(thisT+hh, :)';
    end
end

%% ---- Prior setup ----
if doRobustPrior
    trainingT = sum(dates_trans(p+1:thisT) < datenum(2020,3,1));
    trainingT = max(trainingT, N+20);
else
    trainingT = T_est;
end

% AR(1) residual variances for Minnesota calibration
sig2 = zeros(N, 1);
p_ar = min(4, size(tmpY,1) - 2);
for i = 1:N
    T_ar = min(trainingT, size(tmpY,1) - p_ar - 1);
    Z_ar = ones(T_ar, 1);
    for lag = 1:p_ar
        Z_ar = [Z_ar, tmpY(p_ar+1-lag : p_ar+T_ar-lag, i)]; %#ok<AGROW>
    end
    y_ar = tmpY(p_ar+2 : p_ar+1+T_ar, i);
    tmpb = (Z_ar'*Z_ar) \ (Z_ar'*y_ar);
    sig2(i) = max(mean((y_ar - Z_ar*tmpb).^2), 1e-10);
end

% Minnesota prior on A (Kronecker structure)
% VA0(row) = c1 / (l^2 * sig2(j))  where c1 = theta(1)*theta(2)
% Maps CCMM's cross-variable shrinkage to Kronecker framework
c1 = theta(1) * theta(2);   % 0.05
c2 = theta(3);              % 100 (intercept)

A0  = zeros(K, N);
VA0 = zeros(K, 1);
VA0(1) = c2;
for i = 2:K
    l   = ceil((i-1) / N);
    idx = mod(i-1, N); if idx == 0, idx = N; end
    VA0(i) = c1 / (l^2 * sig2(idx));
end

% Minnesota prior mean: own lag-1 coefficient
for i = 1:N
    A0(1+i, i) = mnPrior(i);   % row 1+i = variable i at lag 1
end

% IW prior on Sigma (data-informed, proper)
nu0 = N + 3;
S0  = diag(sig2);   % prior mean of Sigma ≈ AR variances

% SV priors (Chan 2020)
nuh0  = 5;
Sh0   = 0.01 * (nuh0 - 1);
rho0  = 0.9;
Vrho  = 0.2^2;
nuub  = 100;

%% ---- MCMC sampler ----
rndStream = RandStream('threefry', 'Seed', thisT);

all_ydraws   = NaN(N, fcstNhorizons, fcstNdraws);
joint_lden   = NaN(fcstNdraws, 1);
draw_counter = 0;

% Initialize chain
h     = zeros(T_est, 1);
nu    = 5;
rho   = 0.8;
sigh2 = 0.1;
lam   = 1 ./ gamrnd(nu/2, 2/nu, T_est, 1);

% Pre-compute sparse prior precision
iVA0 = spdiags(1./VA0, 0, K, K);

m = 0;
while m < MCMCreps
    m = m + 1;

    %% Step 1: Sample Sigma and A jointly (Normal-IW conjugate)
    iOm  = spdiags(exp(-h) ./ lam, 0, T_est, T_est);
    XiOm = X' * iOm;
    KA   = iVA0 + XiOm * X;
    Ahat = KA \ (iVA0 * A0 + XiOm * shortY);
    Shat = S0 + A0' * iVA0 * A0 + shortY' * iOm * shortY - Ahat' * KA * Ahat;
    Shat = (Shat + Shat') / 2;

    try
        Sig  = iwishrnd(Shat, nu0 + T_est);
        CSig = chol(Sig, 'lower');
        A    = Ahat + (chol(KA, 'lower')' \ randn(rndStream, K, N)) * CSig';
    catch
        % Keep previous Sig/A if iwishrnd fails (very rare)
        if m == 1
            Sig  = diag(sig2);
            CSig = chol(Sig, 'lower');
            A    = A0;
        end
    end

    %% Step 2: Sample h_t (scalar common SV, Newton-MH)
    U   = shortY - X * A;
    tmp = U / CSig';                % T x N: standardised residuals
    s2  = sum(tmp.^2, 2) ./ lam;   % T x 1
    [h, ~] = sample_h(s2, rho, sigh2, h, N);

    %% Step 3: Sample lambda_t
    s2_lam = sum(tmp.^2, 2) ./ exp(h);
    lam    = 1 ./ gamrnd((N + nu)/2, 2 ./ (s2_lam + nu));

    %% Step 4: Sample nu (MH via sample_nu.m)
    [nu, ~, ~] = sample_nu(lam, nu, nuub);

    %% Step 5: Sample sigh2
    eh    = [h(1)*sqrt(1-rho^2); h(2:end) - rho*h(1:end-1)];
    sigh2 = 1 / gamrnd(nuh0 + T_est/2, 1 / (Sh0 + sum(eh.^2)/2));

    %% Step 6: Sample rho (truncated Normal MH)
    Krho   = 1/Vrho + sum(h(1:T_est-1).^2) / sigh2;
    rhohat = Krho \ (rho0/Vrho + h(1:T_est-1)' * h(2:T_est) / sigh2);
    rhoc   = rhohat + sqrt(Krho)' \ randn(rndStream);
    grho   = @(x) -0.5*log(sigh2./(1-x.^2)) - 0.5*(1-x.^2)/sigh2 * h(1)^2;
    if abs(rhoc) < 0.9999
        alpMH = exp(grho(rhoc) - grho(rho));
        if alpMH > rand(rndStream)
            rho = rhoc;
        end
    end

    %% Post burn-in: forecast
    if m > burnin
        thisdraw = m - burnin;

        h_last = h(end);            % scalar: log-vol at time T
        mu_fcst = (Xjumpoff' * A)'; % N x 1: conditional mean at h=1

        for nn = 1:Ndraws
            draw_counter = draw_counter + 1;

            fcstX = Xjumpoff;
            h_cur = h_last;

            for hh = 1:fcstNhorizons
                % Propagate scalar SV: AR(1)
                h_cur = rho * h_cur + sqrt(sigh2) * randn(rndStream);

                % Draw scalar lambda
                lam_fcst = 1 / gamrnd(nu/2, 2/nu);

                % Total scale: sqrt(lambda * exp(h))
                total_scale = sqrt(lam_fcst * exp(h_cur));

                % Draw forecast: y = X*A + CSig * total_scale * e
                EY    = (fcstX' * A)';  % N x 1
                ydraw = EY + total_scale * CSig * randn(rndStream, N, 1);

                all_ydraws(:, hh, draw_counter) = ydraw;

                % Joint log score at h=1
                if hh == 1 && ~any(isnan(yrealized(:,1)))
                    dev = yrealized(:,1) - mu_fcst;
                    e_std = CSig \ dev;   % Sigma^{-1/2} * dev
                    logdetSig = 2 * sum(log(diag(CSig)));
                    total_var = lam_fcst * exp(h_cur);
                    joint_lden(draw_counter) = -0.5 * (Nlogtwopi ...
                        + N*log(total_var) + logdetSig ...
                        + (1/total_var) * (e_std' * e_std));
                end

                % Update lag state
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
