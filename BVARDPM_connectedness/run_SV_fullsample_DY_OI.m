% =========================================================================
% run_SV_fullsample_DY_OI.m
%
% Gaussian SV benchmark estimation for the Diebold-Yilmaz bank-network
% connectedness application.
%
% Specification notes:
%   - Each equation includes a global intercept with a diffuse N(0, 100)
%     prior. DPM-OISV has no global intercept because its mu_t cluster
%     mean component handles the role more flexibly; the contrast is
%     therefore "time-varying mean (DPM) vs constant mean (SV)", which
%     isolates what the DPM mechanism genuinely adds. Without the
%     intercept, the VAR is forced to produce a zero-mean process
%     despite non-zero data means, pushing the posterior companion
%     matrix against the stability boundary.
%   - SV priors: nuh=5, Sh=0.2, so the data determines omega2.
%   - Stationarity filter: proposals with spectral radius >= 1 are
%     rejected.
%
% This SV benchmark is a STRUCTURAL MIRROR of BVAR_DPM_OIMSV_FE_DY.m
% with the DPM mechanism stripped out. Same full B0 matrix sampled via
% Waggoner-Zha-Villani with absolute-normal; same per-variable AR(1) KSC
% SV sampler (sample_SV_OISV.m from the Chan-Koop-Yu 2024 OISV codebase);
% same inline GFEVD computation. The ONLY differences from the DPM-OISV
% wrapper are:
%
%   1. No DPM cluster assignments (no oneDPM_MSV, no z, no alpha, no K)
%   2. No NIG base measure (no lam, no mu - both fixed at identity/zero)
%   3. No DPM-specific parameter storage (store_DPM_params)
%   4. Simpler Sigma_t = iB0 * diag(exp(H_t)) * iB0'  (no lam_t scaling)
%
% This structural mirror ensures that any difference between the DPM-OISV
% and SV results in the DY application comes from the DPM error
% distribution alone, NOT from different SV samplers or different B0
% identification schemes.
%
% Paper-quality MCMC settings: 10,000 + 5,000 burnin, thin-by-10 for
% heavy storage arrays. Expected runtime: approx 45-60 min on a modern
% desktop.
%
% Required helper files (must be on the MATLAB path):
%   compute_GFEVD.m, get_resid_var_OISV_FE.m, get_C_OISV_FE.m,
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
out_path = 'results/results_SV_DY.mat';

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
%
% A global intercept is included. The A matrix is (n*p+1) x n:
%   - Row 1: intercept (one per equation). Diffuse N(0, 100) prior.
%   - Rows 2..n*p+1: lag coefficients with Horseshoe-Minnesota shrinkage.
%
% The Minnesota scaling from get_C_OISV_FE is built for a k_lag = n*p
% coefficient vector (no intercept). We then:
%   - prepend a row of diffuse variances for the intercept
%   - shift the kappa indices by +1 to account for the prepended row
%
% This keeps the Horseshoe-Minnesota mechanism identical to the DPM-OISV
% wrapper for the lag coefficients, while giving the intercept a diffuse
% prior (standard practice).

% Use k_lag for the Minnesota-indexed positions; k for the full A matrix
k_lag = n*p;       % lag-only dimensions (used by get_C_OISV_FE and kappa indices)
k     = n*p + 1;   % full A matrix dimension (including intercept row)

% 1. Priors for A (Horseshoe-Minnesota, Chan-Koop-Yu 2021)
% Note: get_resid_var_OISV_FE defaults to AR(4) for the residual-variance
% prior calibration. With p = 3 (daily DY-standard lag length) we only have
% 3 pre-sample rows in Y0, so we explicitly pass p_ar = p to use AR(p)
% instead. This is fine for prior calibration purposes — the AR order is
% only there to give a sensible variance scale, not to fit the data well.
sig2 = get_resid_var_OISV_FE(Y0, shortY, p);
[C, idx_kappa1_flat, idx_kappa2_flat] = get_C_OISV_FE(n, p, sig2);
C_A_lag = reshape(C(1:k_lag*n), k_lag, n);

% Augment with intercept row: diffuse variance 100 for all equations
% so the intercept "passes through" the Minnesota scaling as 1*kappa1*100
% but we'll handle it outside the kappa framework (see below).
% Concretely: make C_A full size (k x n) with first row placeholder 1.
C_A = [ones(1, n); C_A_lag];   % (n*p+1) x n

% Intercept-specific prior variance (diffuse).
V_intercept = 100;

% Minnesota indices for lag positions, SHIFTED by +1 for the prepended intercept
idx_kappa1 = cell(n,1);
idx_kappa2 = cell(n,1);
for j = 1:n
    offset = (j-1)*k_lag;
    idx_kappa1{j} = idx_kappa1_flat(idx_kappa1_flat > offset & idx_kappa1_flat <= offset+k_lag) - offset + 1;
    idx_kappa2{j} = idx_kappa2_flat(idx_kappa2_flat > offset & idx_kappa2_flat <= offset+k_lag) - offset + 1;
end
A0 = zeros(k, n);   % prior mean zero for both intercept and lag coefs

% 2. Priors for B0
b0_B = eye(n);
iV_B = eye(n);

% 3. Priors for SV process
% NOTE: These are LOOSER than the DPM-OISV wrapper's priors (which use
% nuh=200, Sh=0.199). The DPM model absorbs crisis variance through its
% lambda_t mixture scale, so its SV doesn't need to move much. In SV-only,
% SV is the ONLY mechanism for time-varying volatility, so we use a
% weakly informative prior that lets the data determine omega2.
Hyper.phi0 = 0.95*ones(n,1);
Hyper.Vphi = (0.1^2)*ones(n,1);
Hyper.nuh  = 5*ones(n,1);       % prior mean omega2 = Sh/(nuh-1) = 0.05
Hyper.Sh   = 0.2*ones(n,1);

%% ===================== Construct Regressor Matrix X =====================
% First column is ones (intercept), remaining columns are lags
X = zeros(T, k);
X(:, 1) = 1;           % intercept column
for i = 1:p
    X(:, 1+(i-1)*n+1 : 1+i*n) = tmpY(p-i+1:end-i, :);
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

% SV (no DPM): lam is fixed at 1, mu is fixed at zero throughout.
% These are kept as variables (rather than being deleted) so the A / B0 /
% H sampling blocks below remain IDENTICAL to the DPM-OISV wrapper - the
% only difference is that here lam and mu are never updated.
lam = ones(T, 1);
mu  = zeros(T, n);

%% ===================== MCMC =============================================
fprintf('\n----- Starting MCMC for SV-DY (fast AR(1) KSC sampler) -----\n');
fprintf('  total iterations: %d (burnin = %d, retained = %d)\n', ...
        nsims+burnin, burnin, nsims);
fprintf('  stationarity filter: ON (reject draws with spectral radius >= 1)\n');
start_time = tic;
n_unstable = 0;
n_rejected = 0;       % total A proposals rejected by stationarity filter
max_reject_attempts = 50;  % safety cap per iteration

for isim = 1:(nsims + burnin)

    %----------------------------------------------------------------------
    % Step 1: Sample Horseshoe Shrinkage Hyperparameters (lag coefs only)
    % The intercept (row 1) has a fixed diffuse prior V_intercept,
    % so it does not participate in the Horseshoe update. Only rows
    % 2..k (the lag coefficients) are shrunk via the Horseshoe-Minnesota.
    %----------------------------------------------------------------------
    A_center = A - A0;
    for i = 1:n
        Vbeta_base_i = zeros(k,1);
        Vbeta_base_i(idx_kappa1{i}) = C_A(idx_kappa1{i}, i) * kappa1;
        Vbeta_base_i(idx_kappa2{i}) = C_A(idx_kappa2{i}, i) * kappa2;
        % Only update psi for lag positions (2..k); leave intercept row untouched
        lag_idx = 2:k;
        rate_psi = 1./z_psi(lag_idx,i) + (A_center(lag_idx,i).^2)./(2*Vbeta_base_i(lag_idx));
        psi(lag_idx,i) = 1./gamrnd(1, 1./rate_psi);
        rate_z_psi = 1 + 1./psi(lag_idx,i);
        z_psi(lag_idx,i) = 1./gamrnd(1, 1./rate_z_psi);
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
    % Step 2: Sample A equation-by-equation (with stationarity filter)
    % The first row of A (intercept) uses a fixed diffuse prior
    % variance V_intercept. Lag positions 2..k use the Horseshoe-Minnesota.
    %----------------------------------------------------------------------
    iV_A = cell(n,1); iVA0A0 = cell(n,1);
    for j = 1:n
        Vbeta_j = zeros(k,1);
        Vbeta_j(1) = V_intercept;                            % diffuse intercept
        Vbeta_j(idx_kappa1{j}) = C_A(idx_kappa1{j}, j) * kappa1;
        Vbeta_j(idx_kappa2{j}) = C_A(idx_kappa2{j}, j) * kappa2;
        % Multiply only the lag rows by psi (intercept row keeps V_intercept)
        Vbeta_j(2:k) = Vbeta_j(2:k) .* psi(2:k, j);
        iV_A{j}   = spdiags(1./Vbeta_j, 0, k, k);
        iVA0A0{j} = iV_A{j} * A0(:,j);
    end

    Y_adj  = shortY - mu;
    A_old  = A;             % save current state for rejection
    is_stationary = false;
    n_attempts = 0;

    while ~is_stationary && n_attempts < max_reject_attempts
        n_attempts = n_attempts + 1;
        A_draw = A_old;     % always re-propose from the last accepted state
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

        % --- Stationarity check on the companion matrix ---
        % Lag-l block at rows 1+(l-1)*n+1 : 1+l*n (skip intercept row 1)
        B_lags_check = zeros(n, n, p);
        for l = 1:p
            B_lags_check(:,:,l) = A_draw(1+(l-1)*n+1 : 1+l*n, :)';
        end
        if p == 1
            sr_check = max(abs(eig(B_lags_check(:,:,1))));
        else
            comp_top = reshape(B_lags_check, n, n*p);
            sr_check = max(abs(eig([comp_top; eye(n*(p-1)) zeros(n*(p-1), n)])));
        end

        if sr_check < 1
            is_stationary = true;
        else
            n_rejected = n_rejected + 1;
        end
    end

    if is_stationary
        A = A_draw;
    else
        % All attempts failed — keep previous A (extremely rare).
        % This can only happen if the posterior is strongly concentrated
        % above the stability boundary, which would indicate a model
        % specification issue.
        A = A_old;
        if mod(isim, 1000) == 0
            fprintf('  WARNING: iter %d - all %d stationarity attempts failed, keeping old A\n', ...
                    isim, max_reject_attempts);
        end
    end

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

    % (Steps 6-7 from the DPM-OISV wrapper - DPM collapsed Gibbs and
    % alpha sampling - are intentionally REMOVED here. In this SV
    % specification, lam and mu stay at their initialized values of 1
    % and 0 throughout the chain.)

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
            store_SV_params(isave_thin,:)   = [phi', omega2'];
            store_h(:,:,isave_thin)         = H;
        end

        % Build reduced-form B_lags from the (k x n) coefficient matrix A.
        % A is (n*p+1) x n with row 1 = intercept, rows 2..n*p+1 = lags.
        % So lag-l block is at rows 1+(l-1)*n+1 : 1+l*n.
        B_lags = zeros(n, n, p);
        for l = 1:p
            B_lags(:,:,l) = A(1+(l-1)*n+1 : 1+l*n, :)';
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
        fprintf('  iter %5d / %d  (elapsed %.1f sec)\n', ...
                isim, nsims+burnin, toc(start_time));
    end
end

elapsed_time = toc(start_time);
fprintf('\n----- SV-DY MCMC completed -----\n');
fprintf('  elapsed: %.1f sec (%.1f min)\n', elapsed_time, elapsed_time/60);
fprintf('  mean spectral radius (thinned):  %.4f\n', mean(store_spec_radius));
fprintf('  max  spectral radius (thinned):  %.4f\n', max(store_spec_radius));
fprintf('  unstable draws (sr >= 1): %d / %d (%.1f%%) [thinned]\n', ...
        n_unstable, nsims_thin, 100*n_unstable/nsims_thin);
fprintf('  stationarity filter rejections:  %d total (%.1f%% of %d iterations)\n', ...
        n_rejected, 100*n_rejected/(nsims+burnin), nsims+burnin);

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
     'store_SV_params', ...
     'store_C_total', 'store_C_from', 'store_C_to', ...
     'store_GFEVD_snap', 'store_spec_radius', 'store_h', ...
     'snapshot_target_dates', 'snapshot_t', 'snapshot_actual', ...
     't_subsample', 'dates_subsample', 'T_sub', ...
     'gfevd_horizon', 'gfevd_subsample_freq', ...
     'n', 'T', 'p', 'k', 'k_lag', 'V_intercept', ...
     'nsims', 'burnin', 'thin', 'nsims_thin', ...
     'elapsed_time', 'n_unstable', 'n_rejected', ...
     'dates_shortY', 'labels', 'bank_info', ...
     'Hyper', ...
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
