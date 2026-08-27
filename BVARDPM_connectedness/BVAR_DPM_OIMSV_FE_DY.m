% =========================================================================
% BVAR_DPM_OIMSV_FE_DY.m
%
% MCMC: 10,000 posterior draws after 5,000 burn-in. Heavy storage arrays
% are thinned by a factor of 10 (1,000 stored samples each);
% store_C_total is kept unthinned (10,000 draws).
%
% DPM-OISV estimation for the Diebold-Yilmaz bank-network connectedness
% application. Adapts the full-sample sampler (BVAR_DPM_OIMSV_FE.m) to:
%
%   1. Load the top-30 DDLY bank panel (log-variance series)
%   2. Compute GFEVD-based connectedness summaries INLINE during the MCMC
%      loop, dropping the heavy storage of H/mu/lam arrays
%
% Reference papers:
%   - Pesaran-Shin (1998) for the generalized FEVD
%   - Diebold-Yilmaz (2012, 2014) for the connectedness measures
%   - Demirer-Diebold-Liu-Yilmaz (2018) for the bank network application
%   - Chan-Yu (2022) for the SV-based dynamic connectedness approach
%
% MODEL (DPM-OISV, no global intercept):
%   y_t - X_t*A - mu_t ~ N(0, lambda_t * Sigma_t)
%   Sigma_t = B0^{-1} * D_t * (B0^{-1})'
%   D_t = diag(exp(h_{1,t}), ..., exp(h_{n,t}))
%   (mu_t, lambda_t) ~ DPM(alpha, NIG)
%
% Required helper files (must be on the MATLAB path):
%   compute_GFEVD.m, get_resid_var_OISV_FE.m, get_C_OISV_FE.m,
%   nigMSV.m, oneDPM_MSV.m, randmn.m,
%   sample_SV_OISV.m, sample_SV0para_OISV.m
%
% =========================================================================

clear; clc; close all;
rng('default'); rng('shuffle');

%% ===================== Configuration ====================================
% MCMC settings (paper-quality: matches GaR DPM convention)
nsims  = 10000;
burnin = 5000;

% Thinning for heavy arrays. We keep every thin-th retained draw for
% storage. store_C_total is kept UNTHINNED (it is small).
thin   = 10;
nsims_thin = nsims / thin;   % 1000 at nsims=10000, thin=10
if mod(nsims, thin) ~= 0
    error('nsims must be divisible by thin');
end

% Connectedness settings
gfevd_horizon         = 10;     % H in DY framework, standard value
gfevd_subsample_freq  = 20;     % monthly for dev pass; set to 5 for weekly
                                % paper pass

% Snapshot dates for full n x n GFEVD storage. Each is mapped to the
% closest actual trading day in dates_shortY at runtime.
snapshot_target_dates = {
    '2007-08-09',  'BNP_liquidity_crisis';
    '2008-09-15',  'Lehman_bankruptcy';
    '2010-05-06',  'Flash_Crash_Greek_crisis';
    '2011-08-08',  'US_downgrade_Euro_crisis';
    '2013-06-19',  'Taper_Tantrum';
};

% Output filename
out_path = 'results/results_DPM_OIMSV_DY.mat';

%% ===================== Load Data ========================================
data_path = 'data/ddly_top30_data.mat';
if ~exist(data_path, 'file')
    error('BVAR_DPM_OIMSV_FE_DY:DataNotFound', ...
        'Cannot find %s. Run load_ddly_data.m first.', data_path);
end

load(data_path, 'Y0', 'shortY', 'dates_shortY', 'labels', 'bank_info', ...
                'n', 'T', 'p', 'k', 'minnesotaPriorMean');

% Defensive casts: when v7.3 .mat files are loaded by Octave (or some older
% MATLAB versions), scalars stored via HDF5 can come back as a non-standard
% subclass that built-ins like eye() do not accept. Forcing a double cast
% is a no-op in modern MATLAB and a safety net everywhere else.
n = double(n);
T = double(T);
p = double(p);
k = double(k);

fprintf('Loaded %s\n', data_path);
fprintf('  n = %d banks, T = %d, p = %d, k = %d\n', n, T, p, k);
fprintf('  Sample: %s to %s\n', ...
    datestr(dates_shortY(1),  'yyyy-mm-dd'), ...
    datestr(dates_shortY(end),'yyyy-mm-dd'));

tmpY = [Y0; shortY];

% Make sure output folder exists
[out_dir, ~, ~] = fileparts(out_path);
if ~isempty(out_dir) && ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

%% ===================== Map snapshot dates to trading days ===============
N_snap = size(snapshot_target_dates, 1);
snapshot_t       = zeros(N_snap, 1);
snapshot_actual  = cell(N_snap, 1);
fprintf('\nSnapshot dates for full n x n GFEVD storage:\n');
for s = 1:N_snap
    target = datenum(snapshot_target_dates{s, 1}, 'yyyy-mm-dd');
    [~, idx] = min(abs(dates_shortY - target));
    snapshot_t(s) = idx;
    snapshot_actual{s} = datestr(dates_shortY(idx), 'yyyy-mm-dd');
    fprintf('  %s (%s) -> %s (t = %d)\n', ...
        snapshot_target_dates{s, 1}, snapshot_target_dates{s, 2}, ...
        snapshot_actual{s}, idx);
end

%% ===================== Subsampled time-grid for GFEVD ===================
% Index into rows of shortY (1..T). Use ceil so we always include t=1.
t_subsample = (1 : gfevd_subsample_freq : T)';
T_sub       = length(t_subsample);
fprintf('\nGFEVD time-grid:\n');
fprintf('  subsample frequency = %d trading days\n', gfevd_subsample_freq);
fprintf('  T_sub = %d evaluation points (out of T=%d)\n', T_sub, T);
fprintf('  GFEVD horizon H = %d steps\n', gfevd_horizon);

%% ===================== Set Prior Parameters =============================

% 1. Priors for A (Horseshoe-Minnesota, Chan-Koop-Yu 2021)
% Note: get_resid_var_OISV_FE defaults to AR(4) for the residual-variance
% prior calibration. With p = 3 (daily DY-standard lag length) we only have
% 3 pre-sample rows in Y0, so we explicitly pass p_ar = p to use AR(p)
% instead. This is fine for prior calibration purposes — the AR order is
% only there to give a sensible variance scale, not to fit the data well.
sig2 = get_resid_var_OISV_FE(Y0, shortY, p);
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

% 3. Priors for SV process (matching the corrected forecast file)
Hyper.phi0 = 0.9*ones(n,1);
Hyper.Vphi = (0.2^2)*ones(n,1);
Hyper.nuh  = 200*ones(n,1);
Hyper.Sh   = 0.199*ones(n,1);

% 4. Priors for DPM (NIG base measure, matching corrected forecast file)
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
% Paper-quality thinning scheme:
%   - store_C_total is kept UNTHINNED (small: T_sub x nsims, ~10 MB at
%     T_sub=134, nsims=10000). Needed for fine-grained time-series plots.
%   - All other "heavy" arrays are thinned by a factor of `thin`, giving
%     nsims_thin = nsims/thin retained samples each.
%
% Heavy arrays (thinned by factor `thin`):
store_A          = zeros(k, n, nsims_thin);
store_B0         = zeros(n, n, nsims_thin);
store_DPM_params = zeros(nsims_thin, 2);      % [alpha, K]
store_SV_params  = zeros(nsims_thin, 2*n);    % [phi(1..n), omega2(1..n)]

% Connectedness summaries (thinning applied to C_from and C_to only;
% C_total is kept unthinned for time-series plotting):
store_C_total = zeros(T_sub, nsims);          % UNTHINNED
store_C_from  = zeros(T_sub, n, nsims_thin);  % thinned
store_C_to    = zeros(T_sub, n, nsims_thin);  % thinned

% Full GFEVD matrices at the snapshot dates only (thinned):
store_GFEVD_snap = zeros(n, n, N_snap, nsims_thin);

% Full SV log-volatility paths (thinned):
store_h = zeros(T, n, nsims_thin);

% Spectral-radius diagnostic (thinned):
store_spec_radius = zeros(nsims_thin, 1);

%% ===================== Initialize the Chain =============================
warnstate = warning('off', 'MATLAB:rankDeficientMatrix');
if rcond(X'*X) < 1e-15
    A = (X'*X + 1e-6*eye(k)) \ (X'*shortY);
else
    A = X \ shortY;
end
warning(warnstate);

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
lam = ones(T,1);
mu  = zeros(T,n);

%% ===================== MCMC =============================================
fprintf('\n----- Starting MCMC for BVAR-DPM-OIMSV-DY -----\n');
fprintf('  total iterations: %d (burnin = %d, retained = %d)\n', ...
        nsims+burnin, burnin, nsims);
start_time = tic;
n_unstable = 0;

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
    rate_k1 = 1/z_k1;
    for i = 1:n
        idx1 = idx_kappa1{i};
        rate_k1 = rate_k1 + sum((A_center(idx1,i).^2)./(2*psi(idx1,i).*C_A(idx1,i)));
    end
    N_k1 = sum(cellfun(@length, idx_kappa1));
    kappa1 = 1/gamrnd(0.5 + 0.5*N_k1, 1/rate_k1);
    z_k1   = 1/gamrnd(1, 1/(1+1/kappa1));

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
        [L_chol, flag_chol] = chol(K_alpha_i, 'lower');
        if flag_chol == 0
            A_hat_i = L_chol' \ (L_chol \ sum_term);
            A_draw(:,i) = A_hat_i + L_chol'\randn(k,1);
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
    mu_sv  = 0;
    for i = 1:n
        try
            H(:,i) = sample_SV_OISV(log(s2(:,i)+1e-8), H(:,i), phi(i), omega2(i), mu_sv);
        catch
            % Keep previous H(:,i) if sampler fails
        end
    end

    %----------------------------------------------------------------------
    % Step 5: Sample SV parameters (phi, omega2)
    %----------------------------------------------------------------------
    [phi, omega2, ~] = sample_SV0para_OISV(H, phi, Hyper);

    %----------------------------------------------------------------------
    % Step 6: DPM collapsed Gibbs for (lam, mu, z)
    %----------------------------------------------------------------------
    X_dp = shortY - X*A;
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

    for i = randperm(T)
        x     = X_dp(i,:);
        h_vec = H(i,:);
        kk    = z(i);
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
        pp = [Pk, P0]; pp = pp - max(pp); pp = exp(pp); pp = pp ./ sum(pp);
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

    % Sample cluster parameters (mu_k, lambda_k) — uses CORRECTED IG sampler
    lamK = zeros(K,1); muK = zeros(K,n);
    for i = 1:K
        Ui  = Theta{i}.U_;  iVi = Theta{i}.iV_;
        mui = Theta{i}.mu_; nui = Theta{i}.nu_;
        S_post = Ui - mui'*iVi*mui/2;
        if S_post <= 0
            S_post = S0_dpm; nui = nu0_dpm;
        end
        lamK(i) = 1 ./ gamrnd(nui, 1/S_post);
        [L_chol, flag_chol] = chol(iVi/lamK(i), 'lower');
        if flag_chol == 0
            muK(i,:) = (mui + L_chol'\randn(n,1))';
        else
            muK(i,:) = mui';
        end
    end

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

    %----------------------------------------------------------------------
    % Step 8 (post-burnin only): Compute GFEVD-based connectedness inline
    %----------------------------------------------------------------------
    if isim > burnin
        isave = isim - burnin;                  % unthinned index (1..nsims)
        is_thin_step = (mod(isave, thin) == 0); % true every thin-th draw
        isave_thin   = isave / thin;            % thinned index (1..nsims_thin)

        % Thinned parameter storage
        if is_thin_step
            store_A(:,:,isave_thin)         = A;
            store_B0(:,:,isave_thin)        = B0;
            store_DPM_params(isave_thin,:)  = [alpha, K];
            store_SV_params(isave_thin,:)   = [phi', omega2'];
            store_h(:,:,isave_thin)         = H;
        end

        % Build reduced-form B_lags from the (k x n) coefficient matrix A.
        % In the DPM-OISV file, X(:, (l-1)*n+1 : l*n) holds y_{t-l}, so the
        % column-i row-(l-1)*n+j element of A is the coefficient on y_{j,t-l}
        % in the equation for y_{i,t}. The transpose aligns with the
        % standard y_t = sum_l B_l y_{t-l} convention.
        B_lags = zeros(n, n, p);
        for l = 1:p
            B_lags(:,:,l) = A((l-1)*n+1 : l*n, :)';
        end

        % Spectral-radius diagnostic on companion form (thinned)
        if is_thin_step
            if p == 1
                comp = B_lags(:,:,1);
            else
                comp_top = reshape(B_lags, n, n*p);
                comp = [comp_top; eye(n*(p-1)) zeros(n*(p-1), n)];
            end
            sr = max(abs(eig(comp)));
            store_spec_radius(isave_thin) = sr;
            if sr >= 1
                n_unstable = n_unstable + 1;
            end
        end

        % Precompute iB0 once per draw — Sigma_t depends on it for every t
        iB0 = B0 \ eye(n);

        % --- Inline GFEVD on the subsampled time grid ---
        for t_idx = 1:T_sub
            t = t_subsample(t_idx);

            % Sigma_t = lambda_t * iB0 * diag(exp(H_t)) * iB0'
            expHt   = exp(H(t,:));               % 1 x n row vector
            scale_t = lam(t) * expHt;            % 1 x n
            Sigma_t = (iB0 .* scale_t) * iB0';
            Sigma_t = 0.5 * (Sigma_t + Sigma_t');

            G = compute_GFEVD(B_lags, Sigma_t, gfevd_horizon);

            % Diebold-Yilmaz summary measures (in percent).
            % C_total is stored every draw; C_from and C_to are thinned.
            diag_G = diag(G);
            store_C_total(t_idx, isave) = (sum(G(:)) - sum(diag_G)) / n * 100;
            if is_thin_step
                store_C_from(t_idx, :, isave_thin) = (sum(G, 2)' - diag_G') / n * 100;
                store_C_to(t_idx, :, isave_thin)   = (sum(G, 1)  - diag_G') / n * 100;
            end
        end

        % --- Full GFEVD matrices at snapshot dates (thinned) ---
        if is_thin_step
            for s = 1:N_snap
                t = snapshot_t(s);
                expHt   = exp(H(t,:));
                scale_t = lam(t) * expHt;
                Sigma_t = (iB0 .* scale_t) * iB0';
                Sigma_t = 0.5 * (Sigma_t + Sigma_t');
                store_GFEVD_snap(:,:,s,isave_thin) = compute_GFEVD(B_lags, Sigma_t, gfevd_horizon);
            end
        end
    end

    if mod(isim, 100) == 0
        fprintf('  iter %5d / %d  (K=%d, elapsed %.1f sec)\n', ...
                isim, nsims+burnin, K, toc(start_time));
    end
end

elapsed_time = toc(start_time);
fprintf('\n----- DPM-OIMSV-DY MCMC completed -----\n');
fprintf('  elapsed: %.1f sec (%.1f min)\n', elapsed_time, elapsed_time/60);
fprintf('  mean spectral radius (thinned):  %.4f\n', mean(store_spec_radius));
fprintf('  unstable draws (sr >= 1): %d / %d (%.1f%%) [thinned]\n', ...
        n_unstable, nsims_thin, 100*n_unstable/nsims_thin);
fprintf('  mean cluster count K: %.1f\n', mean(store_DPM_params(:,2)));
fprintf('  mean alpha: %.3f\n', mean(store_DPM_params(:,1)));

% Quick summary of total connectedness
fprintf('\n  Total connectedness (posterior mean):\n');
C_total_mean = mean(store_C_total, 2);
fprintf('    full sample mean: %.2f%%\n', mean(C_total_mean));
fprintf('    sample range:     [%.2f, %.2f]%%\n', min(C_total_mean), max(C_total_mean));

%% ===================== Save Results =====================================
% Note: gfevd_subsample_freq is also saved so post-processing knows which
% calendar dates the T_sub time grid corresponds to.
dates_subsample = dates_shortY(t_subsample);

save(out_path, ...
     'store_A', 'store_B0', ...
     'store_DPM_params', 'store_SV_params', ...
     'store_C_total', 'store_C_from', 'store_C_to', ...
     'store_GFEVD_snap', 'store_spec_radius', 'store_h', ...
     'snapshot_target_dates', 'snapshot_t', 'snapshot_actual', ...
     't_subsample', 'dates_subsample', 'T_sub', ...
     'gfevd_horizon', 'gfevd_subsample_freq', ...
     'n', 'T', 'p', 'k', 'nsims', 'burnin', 'thin', 'nsims_thin', ...
     'elapsed_time', 'n_unstable', ...
     'dates_shortY', 'labels', 'bank_info', ...
     'Hyper', 'a0_dp', 'b0_dp', 'mu0_dp', 'V0_dp', 'nu0_dpm', 'S0_dpm', ...
     '-v7.3');
fprintf('\nResults saved to %s\n', out_path);

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
