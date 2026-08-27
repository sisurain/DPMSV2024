% =========================================================================
% STEP1_run_SVo_fullsample.m
%
% Runs the SVo (BVAR with common outlier-bar scaling) model on the full
% CCMM 16-variable sample and extracts the fixed DGP parameters needed
% for the Monte Carlo simulation in Steps 2/3.
%
% This runs mcmcVARSVobar.m ONCE on the FULL sample (thisT = T_total),
% NOT in a quasi-real-time loop like goVARSVobar.m.
% We only need one set of full-sample posterior means.
%
% REQUIRED CCMM FILES (must be on MATLAB path or in same folder):
% -----------------------------------------------------------------
%  Core estimation:
%    mcmcVARSVobar.m       -- main Gibbs sampler (SVo model)
%    CTA.m                 -- triangular algorithm for VAR coefficients
%    obarGibbsdraw.m       -- draws common outlier scale o_t
%    StochVolKSCcorrsqrt.m -- KSC mixture SV sampler (correlated, sqrt)
%    getKSC7values.m       -- KSC 7-component mixture constants
%    KSC.m                 -- KSC offset constant
%    betadraw.m            -- beta distribution draws
%
%  Labels / utilities (used only for diagnostic plots, not estimation):
%    fredMDprettylabel.m   -- pretty variable labels for plots
%
%  NOT needed from matlabtoolbox (we replace inline below):
%    parid, getparpoolsize, initRandStreams,
%    xtickdates, plotCI, initwrap, finishwrap, ltitr, progressbar
%
% OUTPUTS (saved to SVo_fullsample_params.mat):
%   dgp_params  -- struct with all fixed inputs for DGP simulation:
%       .A_hat        K x n   posterior mean VAR coefficients (incl. intercept)
%       .invA_hat     n x n   posterior mean impact matrix A^{-1}
%       .Sigma_hat    n x n   base covariance (constructed from SVo outputs)
%       .d_hat        n x 1   per-variable mean variance (sqrtht^2 time avg)
%       .phi_hat      scalar  AR(1) persistence of common h_t
%       .sigh2_hat    scalar  AR(1) innovation variance of common h_t
%       .h_common     T x 1   common log-volatility path (for diagnostics)
%       .n, .T, .p, .K        dimensions
%   Also saves raw MCMC draws: PAI_all, invA_all, sqrtht_all, PHI_all,
%   SVobarscale_all, SVobarprob_all (for further inspection if needed)
%
% PREREQUISITE:
%   Run STEP0_prepare_CCMM_data.m first to generate CCMM16_data.mat
% =========================================================================

clear; clc;
fprintf('=== STEP 1: SVo Full-Sample Estimation ===\n\n');

%% ---- 0. Add CCMM code to path ----
% Adjust these paths to wherever your CCMM .m files live.
% All files listed in "REQUIRED CCMM FILES" above must be findable.
addpath(pwd);              % current folder (if files are here)
% addpath('../CCMM_code'); % example: adjust if files are in a subfolder

%% ---- 1. Load prepared data ----
load('CCMM16_data.mat');   % loads: Y0, shortY, dates_shortY, n, T, p,
                            %        minnesotaPriorMean, labels, codes

fprintf('Data loaded: n=%d variables, T=%d obs, p=%d lags\n', n, T, p);
fprintf('Sample: %s to %s\n\n', ...
    datestr(dates_shortY(1),'yyyy:mmm'), datestr(dates_shortY(end),'yyyy:mmm'));

% Reconstruct full data matrix in the format mcmcVARSVobar expects:
%   data = [Y0; shortY]  -->  (p+T) x n,   ydates = dates for all rows
data   = [Y0; shortY];              % (p+T) x n
ydates = [dates_Y0; dates_shortY];  % (p+T) x 1

Tdata  = size(data, 1);             % p + T = total rows

%% ---- 2. MCMC settings ----
MCMCdraws         = 2000;    % posterior draws to keep (burn-in is 20% extra)
np                = 12;      % periods per year (monthly)
SVOmaxscale       = 20;      % maximum outlier scale (o_t in {1,...,20})
check_stationarity = 0;       % do not enforce stationarity (0 = faster)
doTightPrior      = false;
doRobustPrior     = false;
fcstNdraws        = MCMCdraws;  % must be a multiple of MCMCdraws (min valid value)
fcstNhorizons     = 1;             % minimal horizon; we discard forecast output anyway

%% ---- 3. SVobar prior parameters (following goVARSVobar.m) ----
%
% p(o_t=1) ~ Beta(alpha, beta) where:
%   alpha = (1/(4*np)) * 10*np = 2.5   (prior mean outlier prob ≈ 0.20)
%   beta  = 10*np - alpha = 117.5
% This is the CCMM standard prior on the common outlier probability.
SVobaralpha = 1 / (4 * np) * 10 * np;    % = 2.5
SVobarbeta  = 10 * np - SVobaralpha;      % = 117.5

%% ---- 4. ELB / yield censoring ----
% CCMM censors yields at ELBbound = 0.25% when making forecasts.
% For parameter extraction (full sample), we do NOT need this, but we
% pass [] to mcmcVARSVobar so the SV sampler is unconstrained.
ELBbound  = [];
ndxYIELDS = [];    % no yield censoring needed for parameter extraction

%% ---- 5. Random stream ----
% mcmcVARSVobar expects a random stream object.
% We use a single default stream (no parallel workers needed here).
rndStream = RandStream('mt19937ar', 'Seed', 2024);
RandStream.setGlobalStream(rndStream);

%% ---- 6. thisT = full sample ----
% In the quasi-real-time loop in goVARSVobar.m, thisT indexes the last
% observation used. Here we use the full sample: thisT = Tdata.
thisT = Tdata;   % = p + T

fprintf('Running SVo MCMC: MCMCdraws=%d, n=%d, T=%d (full sample)\n', ...
    MCMCdraws, n, T);
fprintf('This may take several minutes for n=16, T=%d, p=%d...\n\n', T, p);

%% ---- 7. Construct yrealized (needed by mcmcVARSVobar signature) ----
% For full-sample estimation we have no "out of sample" realized values.
% Pass NaN matrix of the correct size.
yrealized = NaN(n, fcstNhorizons);

%% ---- 8. Run SVo estimation ----
tic;
[PAI_all, PHI_all, invA_all, sqrtht_all, ...
    SVobarprob_all, SVobarscale_all, ...
    ~, ~, ~, ~, ~, ~] = ...
    mcmcVARSVobar(thisT, MCMCdraws, p, np, data, ydates, ...
        minnesotaPriorMean, doTightPrior, doRobustPrior, ...
        SVobaralpha, SVobarbeta, SVOmaxscale, ...
        ndxYIELDS, ELBbound, ...
        check_stationarity, ...
        yrealized, ...
        fcstNdraws, fcstNhorizons, rndStream);
elapsed = toc;
fprintf('\nSVo MCMC completed in %.1f seconds (%.1f minutes)\n\n', ...
    elapsed, elapsed/60);

% Sizes after MCMC:
% PAI_all      : K x n x MCMCdraws   (VAR coefficients including intercept)
% invA_all     : n x n x MCMCdraws   (impact matrix A^{-1})
% sqrtht_all   : n x T x MCMCdraws   (sqrt(h_{j,t}), variable-specific SV)
%                NOTE: sqrtht includes the outlier scale o_t in CCMM's coding,
%                      i.e., sqrtht(:,t,m) = sqrtht_pure(:,t,m) * SVobarscale(:,t,m)
% SVobarscale_all: 1 x T x MCMCdraws (common outlier scale o_t, scalar per t)
% SVobarprob_all : 1 x MCMCdraws     (outlier probability p_o)

[K, ~] = size(PAI_all(:,:,1));

fprintf('Output dimensions:\n');
fprintf('  PAI_all:         %s\n', mat2str(size(PAI_all)));
fprintf('  invA_all:        %s\n', mat2str(size(invA_all)));
fprintf('  sqrtht_all:      %s\n', mat2str(size(sqrtht_all)));
fprintf('  SVobarscale_all: %s\n', mat2str(size(SVobarscale_all)));

%% ---- 9. Extract pure SV (strip out the outlier scale from sqrtht) ----
%
% In mcmcVARSVobar, sqrtht(:,t,m) already includes the outlier scale:
%   sqrtht(j,t,m) = sqrtht_pure(j,t,m) * o_t(m)
% where sqrtht_pure is the pure AR(1) SV component.
%
% For our Sigma_hat construction, we want the PURE SV component
% (the persistent volatility, not contaminated by transitory outliers),
% because in the DGP we model lambda_t separately.
%
% Pure SV: sqrtht_pure = sqrtht ./ SVobarscale
% SVobarscale is 1 x T x MCMCdraws, so broadcast across n.

sqrtht_pure_all = sqrtht_all ./ SVobarscale_all;  % n x T x MCMCdraws
% Now sqrtht_pure_all(j,t,m) = exp(h_{j,t}/2) from the AR(1) process

%% ---- 10. Compute posterior means ----

A_hat        = mean(PAI_all,          3);   % K x n
invA_hat     = mean(invA_all,         3);   % n x n
sqrtht_pm    = mean(sqrtht_pure_all,  3);   % n x T  (posterior mean of pure sqrt(h))
SVobarprob_pm = mean(SVobarscale_all, 3);   % 1 x T  (posterior mean of o_t)

%% ---- 11. Construct Sigma_hat ----
%
% The DGP is:  eps_t ~ N(0, lambda_t * exp(h_t) * Sigma_hat)
% where exp(h_t) is a SCALAR common to all variables.
%
% In SVo, the analogous object is:
%   Var(eps_t) = invA * diag(sqrtht_pure(:,t).^2) * invA'
%             = invA * diag(exp(h_{j,t})) * invA'
%
% To construct Sigma_hat, we use the TIME-AVERAGE of the per-variable
% variances as the "normal-period" baseline:
%   d_hat(j) = (1/T) * sum_t sqrtht_pm(j,t)^2
%
% Then:
%   Sigma_hat = invA_hat * diag(d_hat) * invA_hat'
%
% Interpretation: Sigma_hat is the average reduced-form covariance matrix
% in a "normal" period (no outlier, at the time-average SV level).
% In the DGP, exp(h_t) then captures DEVIATIONS from this average,
% with h_t mean-zero by construction (AR(1) with zero unconditional mean).

d_hat     = mean(sqrtht_pm .^ 2, 2);    % n x 1, time-average variance per variable
Sigma_hat = invA_hat * diag(d_hat) * invA_hat';
Sigma_hat = (Sigma_hat + Sigma_hat') / 2;   % symmetrize for numerics

fprintf('Sigma_hat constructed:\n');
fprintf('  Is PD: %d  (all eigenvalues > 0)\n', all(eig(Sigma_hat) > 0));
fprintf('  Diag range: [%.6f, %.6f]\n', min(diag(Sigma_hat)), max(diag(Sigma_hat)));

%% ---- 12. Construct common h_t path and fit AR(1) ----
%
% SVo has n=16 variable-specific SV processes h_{j,t}.
% The common-SV DGP needs ONE scalar phi_hat and sigh2_hat.
%
% Strategy:
%   (a) Mean-centre each variable's log-variance path:
%       hj_centred(t) = 2*log(sqrtht_pm(j,t)) - log(d_hat(j))
%       => hj_centred has unconditional mean ~0 by construction
%   (b) Take cross-sectional average:
%       h_common(t) = (1/n) * sum_j hj_centred(t)
%   (c) Fit AR(1) by OLS to h_common to get phi_hat, sigh2_hat

% (a) Log variance, mean-centred
log_var_pm = 2 * log(sqrtht_pm);              % n x T  (log of per-variable variance)
h_mat      = log_var_pm - log(d_hat);         % n x T  (mean-centred)

% (b) Cross-sectional average -> common path
h_common = mean(h_mat, 1)';    % T x 1

% (c) OLS AR(1): h_t = phi * h_{t-1} + u_t
h_lag  = h_common(1:end-1);    % (T-1) x 1
h_now  = h_common(2:end);      % (T-1) x 1

phi_hat   = (h_lag' * h_lag) \ (h_lag' * h_now);
resid_h   = h_now - phi_hat * h_lag;
sigh2_hat = var(resid_h, 1);    % MLE variance

fprintf('\nAR(1) fit to common h_t:\n');
fprintf('  phi_hat   = %.4f\n', phi_hat);
fprintf('  sigh2_hat = %.4f\n', sigh2_hat);

if phi_hat >= 1.0
    warning('phi_hat=%.4f >= 1. Clipping to 0.999 for stationarity.', phi_hat);
    phi_hat = 0.999;
end

%% ---- 13. Pack dgp_params struct ----
dgp_params.A_hat      = A_hat;        % K x n   -- VAR coefficients (incl. intercept)
dgp_params.invA_hat   = invA_hat;     % n x n   -- impact matrix
dgp_params.Sigma_hat  = Sigma_hat;    % n x n   -- base covariance for DGP
dgp_params.d_hat      = d_hat;        % n x 1   -- per-variable mean variance
dgp_params.phi_hat    = phi_hat;      % scalar  -- SV AR(1) persistence
dgp_params.sigh2_hat  = sigh2_hat;    % scalar  -- SV innovation variance
dgp_params.h_common   = h_common;     % T x 1   -- common SV path (for diagnostics)
dgp_params.n          = n;            % 16
dgp_params.T          = T;            % estimation sample length
dgp_params.K          = K;            % = n*p + 1
dgp_params.p          = p;            % 12
dgp_params.labels     = labels;       % variable names
dgp_params.dates_shortY = dates_shortY; % for plotting

%% ---- 14. Save ----
save('SVo_fullsample_params.mat', ...
    'dgp_params', ...
    'PAI_all', 'invA_all', 'sqrtht_all', 'sqrtht_pure_all', ...
    'PHI_all', 'SVobarscale_all', 'SVobarprob_all', ...
    'n', 'T', 'p', 'K', 'elapsed', '-v7.3');
fprintf('\nSaved: SVo_fullsample_params.mat\n');
fprintf('Next step: run STEP2_generate_sim_datasets.m\n\n');

%% ---- 15. Diagnostic plots ----

% Plot common h_t path
figure('Name','Common log-volatility path h_t','NumberTitle','off');
plot(dates_shortY, h_common, 'b-', 'LineWidth', 1.2);
datetick('x','yyyy','keeplimits');
xlabel('Date'); ylabel('h_t (common)');
title(sprintf('Common log-volatility: \\phi=%.3f, \\sigma_h^2=%.3f', ...
    phi_hat, sigh2_hat));
box off; grid on;

% Plot individual variable SV paths (posterior means)
figure('Name','Variable-specific SV paths (posterior mean)','NumberTitle','off');
for j = 1:n
    subplot(4,4,j);
    plot(dates_shortY, sqrtht_pm(j,:), 'r-', 'LineWidth', 0.7);
    datetick('x','yy','keeplimits');
    title(labels{j}, 'FontSize', 6, 'Interpreter','none');
    box off;
end
sgtitle('Posterior mean of sqrtht per variable (pure SV, outliers stripped)');

% Plot SVo outlier scale (posterior mean of o_t)
figure('Name','SVo common outlier scale o_t','NumberTitle','off');
plot(dates_shortY, squeeze(mean(SVobarscale_all,3)), 'k-', 'LineWidth', 1.2);
datetick('x','yyyy','keeplimits');
xlabel('Date'); ylabel('o_t (outlier scale)');
title('Posterior mean of SVo common outlier scale o_t'); box off;

fprintf('=== STEP 1 complete ===\n\n');