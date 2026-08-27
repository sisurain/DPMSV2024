% =========================================================================
% BVAR_DPM_OIMSV_FE.m
%
% Full-sample estimation of the Bayesian VAR with Dirichlet Process
% Mixture innovations and Order-Invariant Multivariate Stochastic
% Volatility (DPM-OIMSV). No global intercept: the DPM cluster means
% mu_k absorb the location.
%
% Error structure:
%   y_t - X_t*A - mu_t ~ N(0, lambda_t * Sigma_t)
%   where Sigma_t = B0^{-1} * D_t * (B0^{-1})'
%         D_t = diag(exp(h_{1,t}), ..., exp(h_{n,t}))
%   and (mu_t, lambda_t) ~ DPM(alpha, G0) with G0 = NIG base measure.
%
% CHANGES (relative to BVAR_DPM_IOMSV_FE_01_fix.m):
%   1. Added storage for H (individual SV paths), SV parameters, mu, z
%   2. Updated for new dataset (CCMM16_data.mat), p=12, nsims=10000
%
% Priors:
%   A      ~ Horseshoe-Minnesota (Chan-Koop-Yu 2021)
%   B0     ~ N(I, I) with sign normalization
%   h_{it} ~ AR(1) with mean 0 (for identification), KSC sampler
%   phi_i  ~ N(0.95, 0.0025) I(|phi|<1)
%   omega2_i ~ IG(3, 0.1)
%   alpha  ~ Gamma(a0_dp, b0_dp)   [Escobar-West]
%   (mu_k, lambda_k) ~ NIG(mu0, V0, nu0_dpm, S0_dpm)
%
% Required files (in same folder or on path):
%   get_resid_var_OISV_FE.m, get_C_OISV_FE.m,
%   nigMSV.m, oneDPM_MSV.m, randmn.m,
%   sample_SV_OISV.m, sample_SV0para_OISV.m,
%   construct_Sigt_DPM_OISV.m
%
% =========================================================================

clear; clc; close all;

rng('default'); rng('shuffle');

%% ===================== Configuration ====================================
nsims   = 10000;   % Post-burn-in draws
burnin  = 5000;    % Burn-in period

%% ===================== Load Data ========================================
if exist('CCMM16_data.mat','file')
    load('CCMM16_data.mat', 'Y0', 'shortY', 'p', 'n', 'T', ...
         'dates_shortY', 'labels', 'minnesotaPriorMean');
    fprintf('Loaded CCMM16_data.mat: n=%d, T=%d, p=%d\n', n, T, p);
else
    error('CCMM16_data.mat not found. Run STEP0_prepare_CCMM_data.m first.');
end

k = n*p;
tmpY = [Y0; shortY];

%% ===================== Set Prior Parameters =============================

% 1. Priors for A (Horseshoe-Minnesota, Chan-Koop-Yu 2021)
sig2 = get_resid_var_OISV_FE(Y0, shortY);
[C, idx_kappa1_flat, idx_kappa2_flat] = get_C_OISV_FE(n, p, sig2);
C_A = reshape(C(1:k*n), k, n);

idx_kappa1 = cell(n,1);
idx_kappa2 = cell(n,1);
for j = 1:n
    offset = (j-1)*k;
    idx_kappa1{j} = idx_kappa1_flat(idx_kappa1_flat > offset & idx_kappa1_flat <= offset+k) - offset;
    idx_kappa2{j} = idx_kappa2_flat(idx_kappa2_flat > offset & idx_kappa2_flat <= offset+k) - offset;
end
A0 = zeros(k, n);

% 2. Priors for B0
b0_B = eye(n);
iV_B = eye(n);

% 3. Priors for SV process (mean-zero, for identification)
%    Tight prior on omega2 (prior mean 1e-3)
Hyper.phi0 = 0.9*ones(n,1);
Hyper.Vphi = (0.2^2)*ones(n,1);
Hyper.nuh  = 200*ones(n,1);
Hyper.Sh   = 0.199*ones(n,1);

% 4. Priors for DPM (NIG base measure)
a0_dp   = 20;
b0_dp   = 8;
mu0_dp  = zeros(n,1);
V0_dp   = 100*eye(n);
nu0_dpm = 5;
S0_dpm  = 45;
prior_dp = nigMSV(mu0_dp, V0_dp, nu0_dpm, S0_dpm);

%% ===================== Construct Regressor Matrix X =====================
X = zeros(T, k);
for i = 1:p
    X(:, (i-1)*n+1 : i*n) = tmpY(p-i+1:end-i, :);
end

%% ===================== Initialize Storage ===============================
store_A          = zeros(k, n, nsims);
store_B0         = zeros(n, n, nsims);
store_H          = zeros(T, n, nsims);
store_lam        = zeros(T, nsims);
store_mu         = zeros(T, n, nsims);
store_z          = zeros(T, nsims);
store_DPM_params = zeros(nsims, 2);      % [alpha, K]
store_SV_params  = zeros(nsims, 2*n);    % [phi(1..n), omega2(1..n)]

%% ===================== Initialize the Chain =============================

% Initialize A (robust OLS/Ridge)
if rcond(X'*X) < 1e-15
    A = (X'*X + 1e-6*eye(k)) \ (X'*shortY);
else
    A = X \ shortY;
end



% Horseshoe hyperparameters
kappa1 = 1; kappa2 = 1; z_k1 = 1; z_k2 = 1;
psi = ones(k, n); z_psi = ones(k, n);

% B0 and SV parameters
B0     = eye(n);
phi    = 0.9*ones(n,1);
omega2 = 0.001*ones(n,1);

% Initialize H from empirical residual variance
resid_init = shortY - X*A;
H = log(var(resid_init, 0, 1) + 1e-6) .* ones(T, n);

% DPM initialization
alpha = 2;
X_dp_init = shortY - X*A;
[z, Theta, nk] = oneDPM_MSV(X_dp_init, B0, H, alpha, prior_dp);
K = length(Theta);
if isrow(z), z = z'; end
lam = ones(T,1); mu = zeros(T,n);

%% ===================== MCMC =============================================
disp('Starting MCMC for BVAR-DPM-OIMSV ...');
start_time = tic;

for isim = 1:(nsims + burnin)

    %----------------------------------------------------------------------
    % Step 1: Sample Horseshoe Shrinkage Hyperparameters
    %----------------------------------------------------------------------
    A_center = A - A0;
    for i = 1:n
        Vbeta_base_i = zeros(k,1);
        Vbeta_base_i(idx_kappa1{i}) = C_A(idx_kappa1{i}, i) * kappa1;
        Vbeta_base_i(idx_kappa2{i}) = C_A(idx_kappa2{i}, i) * kappa2;
        rate_psi = 1./z_psi(:,i) + (A_center(:,i).^2)./(2*Vbeta_base_i);
        psi(:,i) = 1./gamrnd(1, 1./rate_psi);
        rate_z_psi = 1 + 1./psi(:,i);
        z_psi(:,i) = 1./gamrnd(1, 1./rate_z_psi);
    end
    % kappa1
    rate_k1 = 1/z_k1;
    for i = 1:n
        idx1 = idx_kappa1{i};
        rate_k1 = rate_k1 + sum((A_center(idx1,i).^2)./(2*psi(idx1,i).*C_A(idx1,i)));
    end
    N_k1 = sum(cellfun(@length, idx_kappa1));
    kappa1 = 1/gamrnd(0.5 + 0.5*N_k1, 1/rate_k1);
    z_k1   = 1/gamrnd(1, 1/(1+1/kappa1));
    % kappa2
    rate_k2 = 1/z_k2;
    for i = 1:n
        idx2 = idx_kappa2{i};
        rate_k2 = rate_k2 + sum((A_center(idx2,i).^2)./(2*psi(idx2,i).*C_A(idx2,i)));
    end
    N_k2 = sum(cellfun(@length, idx_kappa2));
    kappa2 = 1/gamrnd(0.5 + 0.5*N_k2, 1/rate_k2);
    z_k2   = 1/gamrnd(1, 1/(1+1/kappa2));

    %----------------------------------------------------------------------
    % Step 2: Sample A equation-by-equation
    %----------------------------------------------------------------------
    iV_A = cell(n,1); iVA0A0 = cell(n,1);
    for j = 1:n
        Vbeta_j = zeros(k,1);
        Vbeta_j(idx_kappa1{j}) = C_A(idx_kappa1{j}, j) * kappa1;
        Vbeta_j(idx_kappa2{j}) = C_A(idx_kappa2{j}, j) * kappa2;
        Vbeta_j = Vbeta_j .* psi(:,j);
        iV_A{j}   = spdiags(1./Vbeta_j, 0, k, k);
        iVA0A0{j} = iV_A{j} * A0(:,j);
    end

    Y_adj  = shortY - mu;
    A_draw = A;
    for i = 1:n
        A_i0 = A_draw;
        A_i0(:,i) = 0;

        Z_i = (Y_adj - X*A_i0) * B0';

        B0_i_sq = B0(:,i).^2;
        B0_i    = B0(:,i);

        S_it = (1./lam) .* sum(exp(-H) .* B0_i_sq', 2);
        M_it = (1./lam) .* sum(Z_i .* exp(-H) .* B0_i', 2);

        X_weighted_S = X .* sqrt(S_it);
        K_alpha_i    = iV_A{i} + X_weighted_S' * X_weighted_S;
        sum_term     = iVA0A0{i} + X' * M_it;
        
        [L, flag_chol] = chol(K_alpha_i, 'lower');
        if flag_chol == 0
            % Solve for mean via triangular systems (no warning)
            A_hat_i = L' \ (L \ sum_term);
            A_draw(:,i) = A_hat_i + L'\randn(k,1);
        end
    end
    A = A_draw;

    %----------------------------------------------------------------------
    % Step 3: Sample B0 (Waggoner-Zha-Villani with absolute-normal)
    %----------------------------------------------------------------------
    U_star  = shortY - X*A - mu;
    B0_draw = B0;
    for i = 1:n
        weights_i  = 1./(lam .* exp(H(:,i)));
        U_weighted = U_star .* sqrt(weights_i);
        K_bi       = iV_B + U_weighted' * U_weighted;

        b0_i    = b0_B(i,:)';
        b_hat_i = K_bi \ (iV_B * b0_i);

        [Ci_L, flag_chol] = chol(K_bi/T, 'lower');
        if flag_chol ~= 0; continue; end

        B0i = B0_draw;  B0i(i,:) = [];
        v_perp = null(B0i);
        if isempty(v_perp); continue; end

        v1     = Ci_L \ v_perp;
        v1     = v1 / norm(v1);
        V_orth = null(v1');

        xi_hat = (b_hat_i' * Ci_L * [v1, V_orth])';

        xi_1      = anormrnd_internal(xi_hat(1), 1/T);
        xi_others = xi_hat(2:end) + sqrt(1/T)*randn(n-1,1);

        b_i_draw = (Ci_L') \ ([v1, V_orth] * [xi_1; xi_others]);

        if b_i_draw(i) < 0
            b_i_draw = -b_i_draw;
        end
        B0_draw(i,:) = b_i_draw';
    end
    B0 = B0_draw;

    %----------------------------------------------------------------------
    % Step 4: Sample H (individual stochastic volatilities)
    %----------------------------------------------------------------------
    U_orth = (shortY - X*A - mu) * B0';
    s2     = U_orth.^2 ./ lam;

    mu_sv = 0;   % mean fixed at 0 for identification
    for i = 1:n
        % Argument order: (ystar, h, rho=phi, sig2=omega2, mu)
        H(:,i) = sample_SV_OISV(log(s2(:,i)+1e-8), H(:,i), phi(i), omega2(i), mu_sv);
    end

    %----------------------------------------------------------------------
    % Step 5: Sample SV parameters (phi, omega2)
    %----------------------------------------------------------------------
    [phi, omega2, ~] = sample_SV0para_OISV(H, phi, Hyper);

    %----------------------------------------------------------------------
    % Step 6: DPM collapsed Gibbs for (lam, mu, z)
    %----------------------------------------------------------------------
    X_dp = shortY - X*A;

    % Reconstruct Theta from scratch
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

    % Sequential collapsed Gibbs sweep
    for i = randperm(T)
        x     = X_dp(i,:);
        h_vec = H(i,:);
        kk    = z(i);

        Theta{kk}.delSample(x, B0, h_vec);
        nk(kk) = nk(kk) - 1;

        if nk(kk) == 0
            Theta(kk) = [];
            nk(kk)    = [];
            K = K - 1;
            z(z > kk) = z(z > kk) - 1;
        end

        if K > 0
            Pk = log(nk) + cellfun(@(t) t.logPredPdf(x, B0, h_vec), Theta);
        else
            Pk = [];
        end
        P0 = log(alpha) + prior_dp.logPredPdf(x, B0, h_vec);
        pp = [Pk, P0];

        pp = pp - max(pp);
        pp = exp(pp);
        pp = pp ./ sum(pp);
        kk = randmn(pp);

        if kk == K + 1
            K = K + 1;
            Theta{kk} = prior_dp.clone();
            Theta{kk}.addSample(x, B0, h_vec);
            nk(kk) = 1;
        else
            Theta{kk}.addSample(x, B0, h_vec);
            nk(kk) = nk(kk) + 1;
        end
        z(i) = kk;
    end

    % Sample cluster parameters (mu_k, lambda_k)
    %   lambda_k | data ~ IG(nu_post, S_post), sampled as 1/gamrnd(nu_post, 1/S_post)
    lamK = zeros(K,1); muK = zeros(K,n);
    for i = 1:K
        Ui = Theta{i}.U_; iVi = Theta{i}.iV_;
        mui = Theta{i}.mu_; nui = Theta{i}.nu_;

        S_post = Ui - mui'*iVi*mui/2;
        if S_post <= 0
            S_post = S0_dpm; nui = nu0_dpm;
        end

        lamK(i) = 1 ./ gamrnd(nui, 1/S_post);

        [L, flag_chol] = chol(iVi/lamK(i), 'lower');
        if flag_chol == 0
            muK(i,:) = (mui + L'\randn(n,1))';
        else
            muK(i,:) = mui';
        end
    end

    % Map to time dimension
    if isrow(z), z = z'; end
    if K > 0
        zK  = sparse((1:T)', z, ones(T,1), T, K);
        lam = zK * lamK;
        mu  = zK * muK;
    end

    %----------------------------------------------------------------------
    % Step 7: Sample alpha (Escobar-West)
    %----------------------------------------------------------------------
    K_adj = max(K, 1);
    xi = betarnd(alpha+1, T);
    pi_alpha = (a0_dp+K_adj-1) / ((a0_dp+K_adj-1) + T*(b0_dp-log(xi)));
    if rand < pi_alpha
        alpha = gamrnd(a0_dp+K_adj,   1/(b0_dp-log(xi)));
    else
        alpha = gamrnd(a0_dp+K_adj-1, 1/(b0_dp-log(xi)));
    end

    %% --- Store draws -----------------------------------------------------
    if isim > burnin
        isave = isim - burnin;
        store_A(:,:,isave)       = A;
        store_B0(:,:,isave)      = B0;
        store_H(:,:,isave)       = H;
        store_lam(:,isave)       = lam;
        store_mu(:,:,isave)      = mu;
        store_z(:,isave)         = z;
        store_DPM_params(isave,:)= [alpha, K];
        store_SV_params(isave,:) = [phi', omega2'];
    end

    if mod(isim, 500) == 0
        fprintf('  Iter %5d / %d  (K=%d)  elapsed %.1f sec\n', ...
                isim, nsims+burnin, K, toc(start_time));
    end
end

elapsed_time = toc(start_time);
fprintf('\nDPM-OIMSV MCMC completed in %.1f seconds (%.1f min).\n', ...
        elapsed_time, elapsed_time/60);

%% ===================== Save Results =====================================
save('results_DPM_OIMSV.mat', ...
     'store_A', 'store_B0', 'store_H', 'store_lam', ...
     'store_mu', 'store_z', 'store_DPM_params', 'store_SV_params', ...
     'n', 'T', 'p', 'k', 'nsims', 'burnin', ...
     'elapsed_time', 'dates_shortY', 'labels', ...
     'Hyper', 'a0_dp', 'b0_dp', 'mu0_dp', 'V0_dp', 'nu0_dpm', 'S0_dpm', ...
     '-v7.3');
fprintf('Results saved to results_DPM_OIMSV.mat\n');

%% ===================== Local Helper Function ============================
function draw = anormrnd_internal(mu, rho)
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
