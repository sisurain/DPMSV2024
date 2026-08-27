% =========================================================================
% run_SVo_sim.m
%
% Runs the CCMM SVo model (mcmcVARSVobar) on all 6 simulation datasets
% and saves per-replication results to model_SVo/results/.
%
% MODEL: BVAR-SV with common outlier scaling (SVo)
%   y_t = X_t * PAI + v_t
%   v_t = inv(A) * diag(sqrt(h_{j,t})) * o_t * e_t,  e_t ~ N(0, I)
%   o_t in {1, 2, ..., 20} with Beta-Bernoulli prior on outlier probability
%   h_{j,t} follows correlated random-walk SV (KSC sampler)
%
% OUTPUT per .mat file:
%   lam_hat   T x 1  posterior mean of lambda_t (o_t^2, variance scale)
%   h_hat     T x 1  posterior mean of common log-volatility (mean-centred)
%   lam_true  T x 1  true lambda_t from DGP
%   h_true    T x 1  true h_t from DGP (mean-centred)
%   dgp, T_sim, r : bookkeeping scalars
%
% FILES REQUIRED in model_SVo/:
%   run_SVo_sim.m             (this file)
%   mcmcVARSVobar.m           (CCMM SVo Gibbs sampler)
%   CTA.m                     (triangular algorithm for VAR coefficients)
%   obarGibbsdraw.m           (common outlier scale sampler)
%   StochVolKSCcorrsqrt.m     (KSC stochastic volatility sampler)
%   KSC.m                     (KSC mixture smoother)
%   getKSC7values.m           (KSC mixture constants)
%   rwnoisePrecisionBasedSampler.m  (precision-based state sampler)
%   betadraw.m, checkdiff.m, parid.m, progressbar.m, vec.m
% =========================================================================

clear; clc;
addpath(pwd, '-begin');

%% ---- 0. Paths ----
data_dir    = '../sim_data/';
results_dir = 'results/';
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
    fprintf('Created results directory: %s\n', results_dir);
end

%% ---- 1. MCMC settings ----
MCMCdraws         = 500;     % posterior draws to keep
p                 = 12;      % VAR lag order
np                = 12;      % periods per year (monthly)
SVOmaxscale       = 20;      % max integer outlier scale
SVobaralpha       = 1/(4*np) * 10*np;   % = 2.5
SVobarbeta        = 10*np - SVobaralpha; % = 117.5
doTightPrior      = false;
doRobustPrior     = false;
check_stationarity = 0;
ndxYIELDS         = [];
ELBbound          = [];
fcstNdraws        = MCMCdraws;
fcstNhorizons     = 1;

%% ---- 2. Minnesota prior means (16 CCMM variables) ----
%  Variable order:
%  1:RPI  2:DPCERA  3:INDPRO  4:CUMFNS  5:UNRATE  6:PAYEMS
%  7:CES0600000007  8:CES0600000008  9:WPSFD49207  10:PCEPI
%  11:HOUST  12:SP500  13:EXUSUKx  14:GS5  15:GS10  16:BAAFFM
n = 16;
minnesotaPriorMean = zeros(n, 1);
persistent_vars = [4, 5, 9, 10, 11, 14, 15, 16];
minnesotaPriorMean(persistent_vars) = 1;

%% ---- 3. Loop over DGPs and sample sizes ----
DGP_list = [1, 2, 3];
T_list   = [300, 800];

for dgp = DGP_list
    for T_sim = T_list

        %% Load simulation file
        fname = sprintf('%ssimdata_DGP%d_T%d_R20.mat', data_dir, dgp, T_sim);
        if ~exist(fname, 'file')
            warning('File not found, skipping: %s', fname);
            continue;
        end
        load(fname, 'all_Y', 'all_lam', 'all_h', 'Y0_fixed', 'R');

        fprintf('\n=== DGP%d  T=%d  (R=%d replications) ===\n', dgp, T_sim, R);

        %% Fixed random stream per (DGP, T_sim) for reproducibility
        base_seed = 1000 * dgp + T_sim;
        rndStream = RandStream('mt19937ar', 'Seed', base_seed);

        for r = 1:R
            fprintf('  Rep %2d/%d ... ', r, R);
            t_start = tic;

            %% Build data matrix: [Y0_fixed; Y] of size (p+T_sim) x n
            Y    = all_Y{r};
            data = [Y0_fixed; Y];

            % Integer ydates (mcmcVARSVobar uses these for sample truncation
            % and trainingT calculation; integers work correctly for both)
            ydates = (1 : p + T_sim)';
            thisT  = p + T_sim;

            % Placeholder for forecast output (discarded)
            yrealized = NaN(n, fcstNhorizons);

            %% Run SVo MCMC
            [~, ~, ~, sqrtht_all, ~, SVobarscale_all, ...
             ~, ~, ~, ~, ~, ~] = ...
                mcmcVARSVobar(thisT, MCMCdraws, p, np, data, ydates, ...
                    minnesotaPriorMean, doTightPrior, doRobustPrior, ...
                    SVobaralpha, SVobarbeta, SVOmaxscale, ...
                    ndxYIELDS, ELBbound, check_stationarity, ...
                    yrealized, fcstNdraws, fcstNhorizons, rndStream);

            %% Extract lambda_hat (posterior mean of o_t^2)
            lam_draws = squeeze(SVobarscale_all).^2;
            lam_hat   = mean(lam_draws, 2);

            %% Extract h_hat (posterior mean of common log-volatility)
            % Strip outlier scale: sqrtht_pure = sqrtht / o_t
            obar_expanded = repmat(SVobarscale_all, n, 1, 1);
            sqrtht_pure   = sqrtht_all ./ obar_expanded;

            % Cross-sectional mean of per-variable log-variance
            log_var_pm = mean(2 * log(sqrtht_pure), 3);
            h_hat_raw  = mean(log_var_pm, 1)';
            h_hat      = h_hat_raw - mean(h_hat_raw);

            %% True values (mean-centre h_true)
            lam_true = all_lam{r};
            h_true   = all_h{r};
            h_true   = h_true - mean(h_true);

            %% Save
            out_file = sprintf('%sSVo_DGP%d_T%d_rep%02d.mat', ...
                                results_dir, dgp, T_sim, r);
            save(out_file, ...
                'lam_hat', 'h_hat', ...
                'lam_true', 'h_true', ...
                'dgp', 'T_sim', 'r');

            fprintf('done (%.1f sec) -> %s\n', toc(t_start), out_file);

        end % r loop
    end % T_sim loop
end % dgp loop

fprintf('\n=== SVo simulation complete. Results in %s ===\n', results_dir);
