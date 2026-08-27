% =========================================================================
% STEP2_generate_sim_datasets.m
%
% Generates R=20 Monte Carlo replications for each combination of
%   DGP  in {1, 2, 3}
%   T_sim in {300, 800}
%
% using the fixed DGP parameters from SVo full-sample estimation
% (from STEP1_run_SVo_fullsample.m).
%
% THE DESIGN (following Prof. Chan's email exactly):
%
%   FIXED across all replications and all DGPs (from SVo posterior means):
%     A_hat       : K x n  VAR coefficients (incl. intercept)
%     Sigma_hat   : n x n  base innovation covariance
%     phi_hat     : scalar SV AR(1) persistence
%     sigh2_hat   : scalar SV AR(1) innovation variance
%     Y0_fixed    : p x n  real initial conditions (first p rows of data)
%
%   RANDOM across replications (different seed each rep):
%     h_t^(r)     : T x 1  common log-volatility path (AR(1) with fixed params)
%     lambda_t^(r): T x 1  transitory scale (DGP-specific distribution)
%     eps_t^(r)   : T x n  innovations ~ N(0, lambda_t * exp(h_t) * Sigma_hat)
%
% THREE DGPs (only lambda_t differs):
%
%   DGP1 - Fat-tail continuous mixture (favours SV-t):
%     lambda_t ~ (1-q)*IG(nu1/2, nu1/2) + q*kappa*IG(nu2/2, nu2/2)
%     q=0.2, nu1=10, nu2=3, kappa=5
%
%   DGP2 - Discrete outliers (favours SVo):
%     lambda_t = 1        w.p. 1-q
%     lambda_t = o_t^2    w.p. q,   o_t ~ U{2,...,20}
%     q=0.05
%     NOTE: stores o_t^2 (variance scale), NOT o_t (std dev scale)
%
%   DGP3 - Clustered pandemic block (stress test):
%     lambda_t = 1             t not in T_out
%     lambda_t = lambda_star   t in T_out = {t0,...,t0+L-1}
%     t0=floor(0.7*T), L=6, lambda_star=25
%
% OUTPUTS:
%   One .mat file per (DGP, T_sim) combination:
%     simdata_DGP{d}_T{T}_R{R}.mat
%   Each file contains:
%     all_Y     : R x 1 cell, each cell T x n
%     all_lam   : R x 1 cell, each cell T x 1 (true lambda_t)
%     all_h     : R x 1 cell, each cell T x 1 (true h_t)
%     dgp_params: the fixed parameter struct (for reference)
%     R, T_sim, n, p, dgp_id, base_seed
%
% PREREQUISITE:
%   Run STEP0 and STEP1 first to generate CCMM16_data.mat and
%   SVo_fullsample_params.mat.
% =========================================================================

clear; clc;
fprintf('=== STEP 2: Generate Simulation Datasets ===\n\n');

%% ---- 1. Load fixed DGP parameters from SVo estimation ----
load('SVo_fullsample_params.mat', 'dgp_params');
load('CCMM16_data.mat', 'Y0');   % p x n real initial conditions

n         = dgp_params.n;        % 16
p         = dgp_params.p;        % 12
K         = dgp_params.K;        % n*p + 1

% Unpack fixed objects used in every replication
A_hat     = dgp_params.A_hat;    % K x n
Sigma_hat = dgp_params.Sigma_hat; % n x n
phi_hat   = dgp_params.phi_hat;  % scalar
sigh2_hat = dgp_params.sigh2_hat; % scalar
Y0_fixed  = Y0;                  % p x n  (first p obs of real data)

% Cholesky factor of Sigma_hat (precomputed, used inside inner loop)
CSig = chol(Sigma_hat, 'lower'); % n x n

fprintf('Fixed parameters loaded:\n');
fprintf('  n=%d, p=%d, K=%d\n', n, p, K);
fprintf('  phi_hat=%.4f, sigh2_hat=%.4f\n', phi_hat, sigh2_hat);
fprintf('  Y0_fixed: %d x %d (real data pre-sample)\n\n', size(Y0_fixed));

%% ---- 2. Simulation design ----
R         = 20;              % replications per (DGP, T_sim)
T_list    = [300, 800];      % two sample sizes
DGP_list  = [1, 2, 3];      % three DGPs
base_seed = 2024;            % base random seed (replication r uses base_seed + (r-1)*100)

% Burn-in for h_t generation (ensures stationarity regardless of T_sim)
T_burnin  = 200;

%% ---- 3. DGP parameter values ----

% DGP1: Fat-tail continuous mixture
dgp1.q     = 0.2;   % mixing probability for heavy component
dgp1.nu1   = 10;    % df for normal component   -> moderate tails
dgp1.nu2   = 3;     % df for heavy component    -> very fat tails
dgp1.kappa = 5;     % scale multiplier for heavy component

% DGP2: Discrete outliers
dgp2.q     = 0.05;  % 5% outlier probability per period

% DGP3: Clustered pandemic block
dgp3.L            = 6;    % block length
dgp3.lambda_star  = 25;   % variance scale inside block (std dev = 5x normal)
% t0 = floor(0.7 * T_sim), computed inside the loop

%% ---- 4. Main generation loop ----

for T_sim = T_list
    for dgp_id = DGP_list

        fprintf('--- DGP %d | T_sim=%d | R=%d replications ---\n', dgp_id, T_sim, R);

        all_Y   = cell(R, 1);
        all_lam = cell(R, 1);
        all_h   = cell(R, 1);

        for r = 1:R

            rng_seed = base_seed + (r - 1) * 100;
            rng(rng_seed);

            % ---- 4a. Generate common log-volatility path h_t ----
            %
            % h_1 ~ N(0, sigh2/(1-phi^2))   [stationary initialisation]
            % h_t = phi * h_{t-1} + u_t,  u_t ~ N(0, sigh2)
            %
            % Generate T_sim + T_burnin periods, discard first T_burnin.

            T_total  = T_sim + T_burnin;
            h_full   = zeros(T_total, 1);
            h_full(1)= sqrt(sigh2_hat / (1 - phi_hat^2)) * randn;
            for t = 2:T_total
                h_full(t) = phi_hat * h_full(t-1) + sqrt(sigh2_hat) * randn;
            end
            h_true = h_full(T_burnin+1:end);   % T_sim x 1

            % ---- 4b. Generate lambda_t ----

            switch dgp_id

                case 1
                    % DGP1: two-component continuous IG mixture
                    % Component indicator: 1 = heavy, 0 = normal
                    ind        = (rand(T_sim, 1) < dgp1.q);
                    lam_norm   = 1 ./ gamrnd(dgp1.nu1/2, 2/dgp1.nu1, T_sim, 1);
                    lam_heavy  = dgp1.kappa ./ gamrnd(dgp1.nu2/2, 2/dgp1.nu2, T_sim, 1);
                    lam_true   = lam_norm .* (1 - ind) + lam_heavy .* ind;

                case 2
                    % DGP2: discrete outliers
                    % lambda_t = 1 (regular) or o_t^2 (outlier, VARIANCE scale)
                    ind      = (rand(T_sim, 1) < dgp2.q);
                    o_draw   = randi([2, 20], T_sim, 1);   % integer from {2,...,20}
                    lam_true = ones(T_sim, 1);
                    lam_true(ind == 1) = o_draw(ind == 1).^2;   % store o_t^2

                case 3
                    % DGP3: clustered pandemic block
                    t0       = floor(0.7 * T_sim);
                    T_out    = t0 : min(t0 + dgp3.L - 1, T_sim);
                    lam_true = ones(T_sim, 1);
                    lam_true(T_out) = dgp3.lambda_star;
            end

            % ---- 4c. Generate VAR observations ----
            %
            % At each t:
            %   x_t = [1, y_{t-1}', ..., y_{t-p}']'   (K x 1)
            %   eps_t = sqrt(lam_t) * exp(h_t/2) * CSig * z_t,  z_t ~ N(0,I)
            %   y_t = A_hat' * x_t + eps_t        (written as row: y_t' = x_t' * A_hat)
            %
            % Initial conditions: Y0_fixed (real data, fixed for all reps)

            Y      = zeros(T_sim, n);
            allY   = [Y0_fixed; Y];   % (p+T_sim) x n, rows 1..p = fixed Y0

            for t = 1:T_sim
                % Build regressor [1, y_{t-1}', ..., y_{t-p}']
                X_t    = zeros(1, K);
                X_t(1) = 1;    % intercept
                for lag = 1:p
                    X_t(1 + (lag-1)*n + (1:n)) = allY(p + t - lag, :);
                end

                % Draw innovation
                scale_t  = sqrt(lam_true(t)) * exp(h_true(t) / 2);
                eps_t    = scale_t * (CSig * randn(n, 1))';   % 1 x n

                Y(t, :)        = X_t * A_hat + eps_t;
                allY(p+t, :)   = Y(t, :);
            end

            % ---- 4d. Store ----
            all_Y{r}   = Y;
            all_lam{r} = lam_true;
            all_h{r}   = h_true;

            if mod(r, 5) == 0
                fprintf('  Rep %d/%d done\n', r, R);
            end

        end % rep loop

        % ---- 4e. Save this (DGP, T_sim) batch ----
        fname = sprintf('simdata_DGP%d_T%d_R%d.mat', dgp_id, T_sim, R);
        save(fname, 'all_Y', 'all_lam', 'all_h', ...
            'dgp_params', 'Y0_fixed', ...
            'R', 'T_sim', 'n', 'p', 'dgp_id', 'base_seed', ...
            'dgp1', 'dgp2', 'dgp3', '-v7.3');
        fprintf('  Saved: %s\n\n', fname);

    end % dgp loop
end % T_sim loop

fprintf('=== STEP 2 complete. All datasets generated and saved. ===\n\n');
fprintf('Files produced:\n');
for T_sim = T_list
    for dgp_id = DGP_list
        fprintf('  simdata_DGP%d_T%d_R%d.mat\n', dgp_id, T_sim, R);
    end
end
fprintf('\nNext step: run each model (DPM-CSV, SVo, SV-t) on each dataset\n');
fprintf('           then compute MSE_h and MSE_lambda.\n');
