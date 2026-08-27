function results = fcst_DPMOISV_oneorigin(thisT, Tdata, data, dates_trans, ...
    N, p, MCMCdraws, burnin, MCMCreps, Ndraws, fcstNdraws, fcstNhorizons, ...
    setQuantiles, Nlogtwopi)
% FCST_DPMOISV_ONEORIGIN  DPM-OISV BVAR estimation + forecast at one origin
%
%   DPM + Order-Invariant Multivariate SV:
%     y_t = X_t * A + mu_{z_t} + u_t
%     u_t ~ N(0, lambda_{z_t} * Sigma_t)
%     Sigma_t = B0^{-1} * D_t * (B0^{-1})'
%     D_t = diag(exp(h_{1,t}), ..., exp(h_{n,t}))
%     h_{i,t} = phi_i * h_{i,t-1} + eta_{i,t}
%     (mu_k, lambda_k) ~ DPM with NIG base measure
%     B0 sampled via Waggoner-Zha-Villani (order-invariant)
%     A ~ Horseshoe-Minnesota prior (no intercept, absorbed by DPM mu)
%
%   nsim=10000, burnin=5000, Ndraws=1 -> 10,000 total draws

%% ---- Prepare data ----
thisdata = data(1:thisT, :);
Y0 = thisdata(1:p, :);          % p x N pre-sample
shortY = thisdata(p+1:end, :);  % T_est x N
T_est = size(shortY, 1);
k = N * p;                       % NO intercept (absorbed by DPM mu)

% Construct X (no intercept)
tmpY = thisdata;
X = zeros(T_est, k);
for i = 1:p
    X(:, (i-1)*N+1 : i*N) = tmpY(p-i+1:end-i, :);
end

% Jumpoff state (no intercept): [y_T, y_{T-1}, ..., y_{T-p+1}]
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

%% ---- Horseshoe-Minnesota prior for A ----
sig2 = get_resid_var_OISV_FE(Y0, shortY);
[C, idx_kappa1_flat, idx_kappa2_flat] = get_C_OISV_FE(N, p, sig2);
C_A = reshape(C(1:k*N), k, N);

idx_kappa1 = cell(N,1);
idx_kappa2 = cell(N,1);
for j = 1:N
    offset = (j-1)*k;
    idx_kappa1{j} = idx_kappa1_flat(idx_kappa1_flat > offset & ...
                    idx_kappa1_flat <= offset+k) - offset;
    idx_kappa2{j} = idx_kappa2_flat(idx_kappa2_flat > offset & ...
                    idx_kappa2_flat <= offset+k) - offset;
end
A0 = zeros(k, N);  % prior mean = 0 (no unit-root; horseshoe handles shrinkage)

%% ---- SV priors ----
Hyper.phi0 = 0.9*ones(N,1);
Hyper.Vphi = (0.2^2)*ones(N,1);
Hyper.nuh  = 200*ones(N,1);
Hyper.Sh   = 0.199*ones(N,1);

%% ---- DPM priors ----
a0_dp   = 20;
b0_dp   = 8;
mu0_dp  = zeros(N, 1);
V0_dp   = 100 * eye(N);
nu0_dpm = 5;
S0_dpm  = 45;

%% ---- B0 prior ----
b0_B = eye(N);
iV_B = eye(N);

%% ---- DPM base measure ----
prior_dp = nigMSV(mu0_dp, V0_dp, nu0_dpm, S0_dpm);

%% ---- MCMC sampler ----
rndStream = RandStream('threefry', 'Seed', thisT);

all_ydraws   = NaN(N, fcstNhorizons, fcstNdraws);
joint_lden   = NaN(fcstNdraws, 1);
draw_counter = 0;

% Initialize A
warnstate = warning('off', 'MATLAB:rankDeficientMatrix');
if rcond(X'*X) < 1e-15
    A = (X'*X + 1e-6*eye(k)) \ (X'*shortY);
else
    A = X \ shortY;
end
warning(warnstate);

% Horseshoe hyperparameters
kappa1 = 1; kappa2 = 1; z_k1 = 1; z_k2 = 1;
psi = ones(k, N); z_psi = ones(k, N);

% B0 and SV
B0     = eye(N);
phi    = 0.9*ones(N,1);
omega2 = 0.001*ones(N,1);

% H initialization
resid_init = shortY - X*A;
H = log(var(resid_init, 0, 1) + 1e-6) .* ones(T_est, N);

% DPM initialization
alpha = 2;
X_dp_init = shortY - X*A;
[z, Theta, nk] = oneDPM_MSV(X_dp_init, B0, H, alpha, prior_dp);
K = length(Theta);
if isrow(z), z = z'; end
lam = ones(T_est, 1);
mu  = zeros(T_est, N);

%% ---- MCMC loop ----
m = 0;
while m < MCMCreps
    m = m + 1;

    %% Step 1: Horseshoe shrinkage hyperparameters
    A_center = A - A0;
    for ii = 1:N
        Vbeta_base_i = zeros(k,1);
        Vbeta_base_i(idx_kappa1{ii}) = C_A(idx_kappa1{ii}, ii) * kappa1;
        Vbeta_base_i(idx_kappa2{ii}) = C_A(idx_kappa2{ii}, ii) * kappa2;
        rate_psi = 1./z_psi(:,ii) + (A_center(:,ii).^2)./(2*Vbeta_base_i);
        psi(:,ii) = 1./gamrnd(1, 1./rate_psi);
        rate_z_psi = 1 + 1./psi(:,ii);
        z_psi(:,ii) = 1./gamrnd(1, 1./rate_z_psi);
    end
    rate_k1 = 1/z_k1;
    for ii = 1:N
        idx1 = idx_kappa1{ii};
        rate_k1 = rate_k1 + sum((A_center(idx1,ii).^2)./(2*psi(idx1,ii).*C_A(idx1,ii)));
    end
    N_k1 = sum(cellfun(@length, idx_kappa1));
    kappa1 = 1/gamrnd(0.5 + 0.5*N_k1, 1/rate_k1);
    z_k1   = 1/gamrnd(1, 1/(1+1/kappa1));

    rate_k2 = 1/z_k2;
    for ii = 1:N
        idx2 = idx_kappa2{ii};
        rate_k2 = rate_k2 + sum((A_center(idx2,ii).^2)./(2*psi(idx2,ii).*C_A(idx2,ii)));
    end
    N_k2 = sum(cellfun(@length, idx_kappa2));
    kappa2 = 1/gamrnd(0.5 + 0.5*N_k2, 1/rate_k2);
    z_k2   = 1/gamrnd(1, 1/(1+1/kappa2));

    %% Step 2: Sample A equation-by-equation
    iV_A = cell(N,1); iVA0A0 = cell(N,1);
    for j = 1:N
        Vbeta_j = zeros(k,1);
        Vbeta_j(idx_kappa1{j}) = C_A(idx_kappa1{j}, j) * kappa1;
        Vbeta_j(idx_kappa2{j}) = C_A(idx_kappa2{j}, j) * kappa2;
        Vbeta_j = Vbeta_j .* psi(:,j);
        iV_A{j}   = spdiags(1./Vbeta_j, 0, k, k);
        iVA0A0{j} = iV_A{j} * A0(:,j);
    end

    Y_adj  = shortY - mu;
    A_draw = A;
    for ii = 1:N
        A_i0 = A_draw; A_i0(:,ii) = 0;
        Z_i = (Y_adj - X*A_i0) * B0';
        B0_i_sq = B0(:,ii).^2;
        B0_i    = B0(:,ii);
        S_it = (1./lam) .* sum(exp(-H) .* B0_i_sq', 2);
        M_it = (1./lam) .* sum(Z_i .* exp(-H) .* B0_i', 2);
        X_weighted_S = X .* sqrt(S_it);
        K_alpha_i    = iV_A{ii} + X_weighted_S' * X_weighted_S;
        sum_term     = iVA0A0{ii} + X' * M_it;
        [L_chol, flag_chol] = chol(K_alpha_i, 'lower');
        if flag_chol == 0
            A_hat_i      = L_chol' \ (L_chol \ sum_term);
            A_draw(:,ii) = A_hat_i + L_chol'\randn(rndStream, k, 1);
        end
    end
    A = A_draw;

    %% Step 3: Sample B0 (Waggoner-Zha-Villani)
    U_star  = shortY - X*A - mu;
    B0_draw = B0;
    for ii = 1:N
        weights_i  = 1./(lam .* exp(H(:,ii)));
        U_weighted = U_star .* sqrt(weights_i);
        K_bi       = iV_B + U_weighted' * U_weighted;
        b0_i    = b0_B(ii,:)';
        b_hat_i = K_bi \ (iV_B * b0_i);
        [Ci_L, flag_chol] = chol(K_bi/T_est, 'lower');
        if flag_chol ~= 0; continue; end
        B0i = B0_draw; B0i(ii,:) = [];
        v_perp = null(B0i);
        if isempty(v_perp); continue; end
        v1     = Ci_L \ v_perp;
        v1     = v1 / norm(v1);
        V_orth = null(v1');
        xi_hat = (b_hat_i' * Ci_L * [v1, V_orth])';
        xi_1      = anormrnd_local(xi_hat(1), 1/T_est);
        xi_others = xi_hat(2:end) + sqrt(1/T_est)*randn(rndStream, N-1, 1);
        b_i_draw = (Ci_L') \ ([v1, V_orth] * [xi_1; xi_others]);
        if b_i_draw(ii) < 0, b_i_draw = -b_i_draw; end
        B0_draw(ii,:) = b_i_draw';
    end
    B0 = B0_draw;

    %% Step 4: Sample H (individual SV per variable)
    U_orth = (shortY - X*A - mu) * B0';  % structural residuals
    s2     = U_orth.^2 ./ lam;
    for ii = 1:N
        try
            H(:,ii) = sample_SV_OISV(log(s2(:,ii)+1e-8), ...
                        H(:,ii), phi(ii), omega2(ii), 0);
        catch
        end
    end

    %% Step 5: Sample SV parameters (phi, omega2)
    [phi, omega2, ~] = sample_SV0para_OISV(H, phi, Hyper);

    %% Step 6: DPM collapsed Gibbs
    X_dp = shortY - X*A;

    % Reconstruct Theta
    K = max(z); if isempty(K), K = 0; end
    Theta = cell(1, K); nk = zeros(1, K);
    for j = 1:K
        idxT = find(z == j);
        nk(j) = length(idxT);
        if nk(j) > 0
            Theta{j} = prior_dp.clone();
            Theta{j}.addData(X_dp(idxT,:), B0, H(idxT,:));
        end
    end

    for ii = randperm(T_est)
        x     = X_dp(ii,:);
        h_vec = H(ii,:);
        kk    = z(ii);
        Theta{kk}.delSample(x, B0, h_vec);
        nk(kk) = nk(kk) - 1;
        if nk(kk) == 0
            Theta(kk) = []; nk(kk) = [];
            K = K - 1; z(z > kk) = z(z > kk) - 1;
        end
        if K > 0
            Pk = log(nk) + cellfun(@(t) t.logPredPdf(x, B0, h_vec), Theta);
        else
            Pk = [];
        end
        P0 = log(alpha) + prior_dp.logPredPdf(x, B0, h_vec);
        pp = [Pk, P0]; pp = pp - max(pp); pp = exp(pp); pp = pp./sum(pp);
        kk = randmn(pp);
        if kk == K + 1
            K = K + 1; Theta{kk} = prior_dp.clone();
            Theta{kk}.addSample(x, B0, h_vec); nk(kk) = 1;
        else
            Theta{kk}.addSample(x, B0, h_vec); nk(kk) = nk(kk) + 1;
        end
        z(ii) = kk;
    end

    % Sample cluster parameters (mu_k, lambda_k)
    lamK = zeros(K,1); muK = zeros(K,N);
    for ii = 1:K
        nui = Theta{ii}.nu_; iVi = Theta{ii}.iV_;
        mui = Theta{ii}.mu_; Ui  = Theta{ii}.U_;
        S_post = Ui - mui'*iVi*mui/2;
        if S_post <= 0, S_post = S0_dpm; nui = nu0_dpm; end
        lamK(ii) = 1 / gamrnd(nui, 1/S_post);
        [L_chol, flag_chol] = chol(iVi/lamK(ii), 'lower');
        if flag_chol == 0
            muK(ii,:) = (mui + L_chol'\randn(rndStream, N, 1))';
        else
            muK(ii,:) = mui';
        end
    end

    if isrow(z), z = z'; end
    if K > 0
        zK  = sparse((1:T_est)', z, ones(T_est,1), T_est, K);
        lam = full(zK * lamK);
        mu  = full(zK * muK);
    end

    %% Step 7: Sample alpha (Escobar-West)
    K_adj = max(K, 1);
    xi_ew = betarnd(alpha+1, T_est);
    pi_alpha = (a0_dp+K_adj-1) / ((a0_dp+K_adj-1) + T_est*(b0_dp-log(xi_ew)));
    if rand(rndStream) < pi_alpha
        alpha = gamrnd(a0_dp+K_adj, 1/(b0_dp-log(xi_ew)));
    else
        alpha = gamrnd(a0_dp+K_adj-1, 1/(b0_dp-log(xi_ew)));
    end

    %% Post burn-in: forecast
    if m > burnin
        draw_counter = draw_counter + 1;

        h_last = H(T_est, :)';   % N x 1 individual log-vols
        iB0    = B0 \ eye(N);    % precompute for forecast
        logdetB0 = log(abs(det(B0)));

        % DPM predictive weights
        dpm_weights = [nk, alpha];
        dpm_weights = dpm_weights / sum(dpm_weights);
        dpm_cdf     = cumsum(dpm_weights);

        for nn = 1:Ndraws
            fcstX = Xjumpoff;
            h_cur = h_last;

            for hh = 1:fcstNhorizons
                % Propagate individual SV: h_{i,T+h} = phi_i*h_{i,T+h-1} + sqrt(omega2_i)*z
                h_cur = phi .* h_cur + sqrt(omega2) .* randn(rndStream, N, 1);

                % Draw from DPM predictive
                u_dpm = rand(rndStream);
                kk_fcst = sum(u_dpm > dpm_cdf) + 1;

                if kk_fcst <= K
                    mu_fcst  = muK(kk_fcst, :)';
                    lam_fcst = lamK(kk_fcst);
                else
                    % New cluster: draw from NIG prior
                    lam_fcst = 1 / gamrnd(nu0_dpm, 1/S0_dpm);
                    L_prior  = chol(V0_dp * lam_fcst, 'lower');
                    mu_fcst  = mu0_dp + L_prior * randn(rndStream, N, 1);
                end

                % Draw y: y = X*A + mu + sqrt(lam)*iB0*diag(exp(h/2))*e
                EY    = A' * fcstX + mu_fcst;
                ydraw = EY + sqrt(lam_fcst) * iB0 * (exp(h_cur/2) .* randn(rndStream, N, 1));

                all_ydraws(:, hh, draw_counter) = ydraw;

                % Joint log score at h=1
                if hh == 1 && ~any(isnan(yrealized(:,1)))
                    dev      = yrealized(:,1) - A' * Xjumpoff - mu_fcst;
                    e_struct = B0 * dev;   % structural residuals
                    logdetSig = N*log(lam_fcst) + sum(h_cur) - 2*logdetB0;
                    quad = (1/lam_fcst) * sum(e_struct.^2 .* exp(-h_cur));
                    joint_lden(draw_counter) = -0.5 * (Nlogtwopi + logdetSig + quad);
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

%% ===================== Local Helper Function ============================
function draw = anormrnd_local(mu, rho)
    % Samples from absolute-normal distribution (positive half)
    w = 1/(1+exp(2*mu/rho));
    if w > rand
        mu1 = mu/2 - sqrt(mu^2+4)/2;
        sig21 = mu1^2*rho/(1+mu1^2);
        draw = mu1 + sqrt(sig21)*randn;
    else
        mu2 = mu/2 + sqrt(mu^2+4)/2;
        sig22 = mu2^2*rho/(1+mu2^2);
        draw = mu2 + sqrt(sig22)*randn;
    end
end
