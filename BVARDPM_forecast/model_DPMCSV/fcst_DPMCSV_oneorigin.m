function results = fcst_DPMCSV_oneorigin(thisT, Tdata, data, dates_trans, ...
    N, p, MCMCdraws, burnin, MCMCreps, Ndraws, fcstNdraws, fcstNhorizons, ...
    setQuantiles, Nlogtwopi)
% FCST_DPMCSV_ONEORIGIN  DPM-CSV BVAR estimation + forecast at one origin
%
%   DPM + Common Scalar Stochastic Volatility:
%     eps_t ~ N(mu_{z_t}, lambda_{z_t} * exp(h_t) * Sigma)
%     h_t = rho * h_{t-1} + eta_t ~ N(0, sigh2)    [scalar AR(1) CSV]
%     (mu_k, lambda_k) ~ DPM with NIG base measure
%     Sigma ~ IW(S0, nu0)                            [Kronecker structure]
%     A ~ Minnesota prior (no intercept)
%
%   nsim=10000, burnin=5000, Ndraws=1 -> 10,000 total draws

%% ---- Prepare data ----
thisdata = data(1:thisT, :);
Y0 = thisdata(1:p, :);          % p x N pre-sample
shortY = thisdata(p+1:end, :);  % T_est x N
T_est = size(shortY, 1);
k = N * p;                       % NO intercept
tmpY = thisdata;

% Construct X (no intercept)
X = zeros(T_est, k);
for i = 1:p
    X(:, (i-1)*N+1 : i*N) = tmpY(p-i+1:end-i, :);
end

% Jumpoff state (no intercept)
Xjumpoff = zeros(k, 1);
for l = 1:p
    Xjumpoff((l-1)*N+(1:N)) = thisdata(thisT-(l-1), 1:N);
end

% Realized future values
yrealized = NaN(N, fcstNhorizons);
for hh = 1:fcstNhorizons
    if thisT + hh <= Tdata
        yrealized(:,hh) = data(thisT+hh, :)';
    end
end

%% ---- Minnesota prior for A (no intercept) ----
[A0, VA0, ~] = construct_prior_A_func(Y0, shortY, p, false);

%% ---- IW prior on Sigma ----
nu0 = N + 3;
S0  = eye(N);

%% ---- SV priors (Chan 2020) ----
nuh0  = 200;
Sh0   = 0.199;
rho0  = 0.9;
Vrho  = 0.2^2;

%% ---- DPM priors ----
a0_dp  = 20;
b0_dp  = 8;
mu0_dp = zeros(N, 1);
V0_dp  = 100 * eye(N);
nu0_dp = 5;
S0_dp  = 45;

%% ---- Pre-compute prior terms ----
iVA0_sparse  = spdiags(1./VA0, 0, k, k);
prior_term_A = iVA0_sparse * A0;
prior_term_S = S0 + A0' * iVA0_sparse * A0;

%% ---- DPM base measure ----
prior_dp = nigSV(mu0_dp, V0_dp, nu0_dp, S0_dp);

%% ---- MCMC sampler ----
rndStream = RandStream('threefry', 'Seed', thisT);

all_ydraws   = NaN(N, fcstNhorizons, fcstNdraws);
joint_lden   = NaN(fcstNdraws, 1);
draw_counter = 0;

% Initialize chain
h     = zeros(T_est, 1);
rho   = 0.9;
sigh2 = 0.001;
lam   = 1 ./ gamrnd(nu0_dp, 1/S0_dp, T_est, 1);
mu    = zeros(T_est, N);
alpha = 2;

% Initial cluster assignments
[z, Theta, nk] = oneDPM(shortY, h, S0, alpha, prior_dp);

% Initialize Sigma
Sig  = eye(N);
CSig = eye(N);

m = 0;
while m < MCMCreps
    m = m + 1;

    %% Step 1: Sample Sigma and A jointly (Normal-IW conjugate)
    iOm_diag = exp(-h) ./ lam;
    iOm      = spdiags(iOm_diag, 0, T_est, T_est);
    XiOm     = X' * iOm;
    Y_adj    = shortY - mu;

    KA   = iVA0_sparse + XiOm * X;
    Ahat = KA \ (prior_term_A + XiOm * Y_adj);
    Shat = prior_term_S + Y_adj' * iOm * Y_adj - Ahat' * KA * Ahat;
    Shat = (Shat + Shat') / 2;

    try
        Sig  = iwishrnd(Shat, nu0 + T_est);
        CSig = chol(Sig, 'lower');
        A    = Ahat + (chol(KA, 'lower')' \ randn(rndStream, k, N)) * CSig';
    catch
        if m == 1
            Sig = eye(N); CSig = eye(N); A = A0;
        end
    end

    %% Step 2: Sample h_t (scalar common SV, Newton-MH)
    U   = shortY - X * A - mu;
    tmp = U / CSig';
    s2  = sum(tmp.^2, 2) ./ lam;
    [h, ~] = sample_h(s2, rho, sigh2, h, N);

    %% Step 3: DPM collapsed Gibbs for (lam, mu, z)
    X_dp = shortY - X * A;

    % Rebuild Theta
    K = max(z); if isempty(K), K = 0; end
    Theta = cell(1, K); nk = zeros(1, K);
    for j = 1:K
        idxT = find(z == j);
        nk(j) = length(idxT);
        if nk(j) > 0
            Theta{j} = prior_dp.clone();
            Theta{j}.addData(X_dp(idxT,:), h(idxT), Sig);
        end
    end

    % Sequential collapsed Gibbs sweep
    for ii = randperm(T_est)
        x  = X_dp(ii,:);
        hi = h(ii);
        kk = z(ii);

        Theta{kk}.delSample(x, hi, Sig);
        nk(kk) = nk(kk) - 1;

        if ~isempty(nk)
            Pk = log(nk) + cellfun(@(t) t.logPredPdf(x, hi, Sig), Theta);
        else
            Pk = [];
        end
        P0 = log(alpha) + prior_dp.logPredPdf(x, hi, Sig);
        pp = [Pk, P0]; pp = pp - max(pp); pp = exp(pp); pp = pp./sum(pp);
        kk = randmn(pp);

        if kk == numel(Theta) + 1
            Theta{kk} = prior_dp.clone();
            Theta{kk}.addSample(x, hi, Sig);
            nk = [nk, 1]; %#ok<AGROW>
        else
            Theta{kk}.addSample(x, hi, Sig);
            nk(kk) = nk(kk) + 1;
        end
        z(ii) = kk;

        % Remove empty clusters
        idx0 = find(nk == 0);
        if ~isempty(idx0)
            Theta(idx0) = []; nk(idx0) = [];
            for jj = idx0
                z(z > jj) = z(z > jj) - 1;
            end
        end
    end

    % Sample cluster parameters (mu_k, lambda_k)
    K = length(Theta);
    lamK = zeros(K,1); muK = zeros(K,N);
    for ii = 1:K
        nui = Theta{ii}.nu_; iVi = Theta{ii}.iV_;
        mui = Theta{ii}.mu_; Ui  = Theta{ii}.U_;
        S_post = Ui - mui'*iVi*mui/2;
        if S_post <= 0, S_post = S0_dp; nui = nu0_dp; end

        % CORRECT: IG(nu_k, S_k) — full shape, not halved
        lamK(ii) = 1 / gamrnd(nui, 1/S_post);

        [L_chol, flag_chol] = chol(iVi/lamK(ii), 'lower');
        if flag_chol == 0
            muK(ii,:) = (mui + L_chol'\randn(rndStream, N, 1))';
        else
            muK(ii,:) = mui';
        end
    end

    if iscolumn(z), z = z'; end
    zK  = sparse(1:T_est, z, ones(T_est,1), T_est, K);
    lam = full(zK * lamK);
    mu  = full(zK * muK);

    %% Step 4: Sample sigh2
    eh    = [h(1)*sqrt(1-rho^2); h(2:end) - rho*h(1:end-1)];
    sigh2 = 1 / gamrnd(nuh0 + T_est/2, 1 / (Sh0 + sum(eh.^2)/2));

    %% Step 5: Sample rho (MH)
    Krho   = 1/Vrho + sum(h(1:T_est-1).^2) / sigh2;
    rhohat = Krho \ (rho0/Vrho + h(1:T_est-1)' * h(2:T_est) / sigh2);
    rhoc   = rhohat + sqrt(1/Krho) * randn(rndStream);
    grho   = @(x) -0.5*log(sigh2./(1-x.^2)) - 0.5*(1-x.^2)/sigh2 * h(1)^2;
    if abs(rhoc) < 0.9999
        alpMH = grho(rhoc) - grho(rho);
        if alpMH > log(rand(rndStream))
            rho = rhoc;
        end
    end

    %% Step 6: Sample alpha (Escobar-West)
    xi_ew    = betarnd(alpha+1, T_est);
    pi_alpha = (a0_dp+K-1) / ((a0_dp+K-1) + T_est*(b0_dp-log(xi_ew)));
    if rand(rndStream) < pi_alpha
        alpha = gamrnd(a0_dp+K, 1/(b0_dp-log(xi_ew)));
    else
        alpha = gamrnd(a0_dp+K-1, 1/(b0_dp-log(xi_ew)));
    end

    %% Post burn-in: forecast
    if m > burnin
        draw_counter = draw_counter + 1;

        h_last = h(end);              % scalar common log-vol
        logdetSig_fixed = 2*sum(log(diag(CSig)));  % log|Sigma|

        % DPM predictive weights
        dpm_weights = [nk, alpha];
        dpm_weights = dpm_weights / sum(dpm_weights);
        dpm_cdf     = cumsum(dpm_weights);

        for nn = 1:Ndraws
            fcstX = Xjumpoff;
            h_cur = h_last;

            for hh = 1:fcstNhorizons
                % Propagate scalar SV: AR(1)
                h_cur = rho * h_cur + sqrt(sigh2) * randn(rndStream);

                % Draw from DPM predictive
                u_dpm = rand(rndStream);
                kk_fcst = sum(u_dpm > dpm_cdf) + 1;

                if kk_fcst <= K
                    mu_fcst  = muK(kk_fcst, :)';
                    lam_fcst = lamK(kk_fcst);
                else
                    % New cluster from NIG prior
                    lam_fcst = 1 / gamrnd(nu0_dp, 1/S0_dp);
                    L_prior  = chol(V0_dp * lam_fcst, 'lower');
                    mu_fcst  = mu0_dp + L_prior * randn(rndStream, N, 1);
                end

                % Total scale: sqrt(lambda * exp(h))
                total_scale = sqrt(lam_fcst * exp(h_cur));

                % Draw y: y = A'*X + mu + total_scale * CSig * e
                EY    = A' * fcstX + mu_fcst;
                ydraw = EY + total_scale * CSig * randn(rndStream, N, 1);

                all_ydraws(:, hh, draw_counter) = ydraw;

                % Joint log score at h=1
                if hh == 1 && ~any(isnan(yrealized(:,1)))
                    dev   = yrealized(:,1) - A' * Xjumpoff - mu_fcst;
                    e_std = CSig \ dev;
                    total_var = lam_fcst * exp(h_cur);
                    joint_lden(draw_counter) = -0.5 * (Nlogtwopi ...
                        + N*log(total_var) + logdetSig_fixed ...
                        + (1/total_var) * (e_std' * e_std));
                end

                % Update state (no intercept)
                fcstX = [ydraw; fcstX(1:end-N)];
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
