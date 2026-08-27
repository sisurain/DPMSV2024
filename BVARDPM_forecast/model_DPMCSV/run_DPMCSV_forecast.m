% =========================================================================
% run_DPMCSV_forecast.m
%
% Expanding-window forecasting for DPM-CSV BVAR.
%
% MODEL:
%   eps_t ~ N(mu_{z_t}, lambda_{z_t} * exp(h_t) * Sigma)
%   h_t = rho * h_{t-1} + eta_t ~ N(0, sigh2)    [scalar common SV]
%   (mu_k, lambda_k) ~ DPM with NIG base measure
%   Sigma ~ IW(S0, nu0)                            [Kronecker structure]
%   A ~ Minnesota prior (no intercept)
%
% MCMC: nsim=10000, burnin=5000, Ndraws=1 -> 10,000 total draws
%
% REQUIRED FILES:
%   fcst_DPMCSV_oneorigin.m, construct_prior_A_func.m,
%   nigSV.m, oneDPM.m, sample_h.m, randmn.m,
%   crpsDraws.m
%   ../data/FREDMD42_forecast_data.mat
% =========================================================================

clear; clc; close all;

%% ---- 1. Settings ----
MCMCdraws          = 1e4;
Ndraws             = 1;
% fcstNhorizons is loaded from the data file (set in STEP0)
setQuantiles       = [5, 16, 50, 84, 95];
Nquantiles         = length(setQuantiles);

burnin     = 5000;
MCMCreps   = MCMCdraws + burnin;
fcstNdraws = MCMCdraws * Ndraws;

%% ---- 2. Load data ----
load('../data/FREDMD42_forecast_data.mat');
fprintf('Loaded FREDMD42_forecast_data.mat: n=%d, T_trans=%d, p=%d\n', n, T_trans, p);

N     = n;
data  = data_trans;
Tdata = T_trans;
Nlogtwopi = N * log(2*pi);

%% ---- 3. Forecast origins ----
Tjumpoffs = fcst_firstOrigin : Tdata;
Njumpoffs = length(Tjumpoffs);
fprintf('Forecast origins: %d  (%s to %s)\n', Njumpoffs, ...
    datestr(dates_trans(Tjumpoffs(1)),'yyyy:mmm'), ...
    datestr(dates_trans(Tjumpoffs(end)),'yyyy:mmm'));

%% ---- 4. Allocate output storage ----
fcstYrealized     = NaN(N, fcstNhorizons, Njumpoffs);
fcstYhat          = NaN(N, fcstNhorizons, Njumpoffs);
fcstYmedian       = NaN(N, fcstNhorizons, Njumpoffs);
fcstCRPS          = NaN(N, fcstNhorizons, Njumpoffs);
fcstLogscore      = NaN(N, fcstNhorizons, Njumpoffs);
fcstJointLogscore = NaN(1, Njumpoffs);
fcstYquantiles    = NaN(N, fcstNhorizons, Nquantiles, Njumpoffs);
fcstPIT           = NaN(N, fcstNhorizons, Njumpoffs);
Ndraws_h1_sub     = fcstNdraws / 10;
fcstYdraws_h1     = NaN(N, Ndraws_h1_sub, Njumpoffs);

%% ---- 5. Parallel pool ----
poolobj = gcp('nocreate');
if isempty(poolobj)
    poolobj = parpool;   % default local pool
end
fprintf('Using %d workers\n\n', poolobj.NumWorkers);

%% ---- 6. Main loop ----
fprintf('Starting DPM-CSV forecast: %d origins x %d MCMC iters x %d draws\n', ...
    Njumpoffs, MCMCreps, Ndraws);
totalTimer = tic;

parfor ndxT = 1 : Njumpoffs

    thisT = Tjumpoffs(ndxT);

    res = fcst_DPMCSV_oneorigin(thisT, Tdata, data, dates_trans, ...
        N, p, MCMCdraws, burnin, MCMCreps, Ndraws, fcstNdraws, fcstNhorizons, ...
        setQuantiles, Nlogtwopi);

    fcstYrealized(:,:,ndxT)    = res.yrealized;
    fcstYhat(:,:,ndxT)         = res.yhat;
    fcstYmedian(:,:,ndxT)      = res.ymedian;
    fcstCRPS(:,:,ndxT)         = res.crps;
    fcstLogscore(:,:,ndxT)     = res.logscore;
    fcstJointLogscore(ndxT)    = res.joint_logscore;
    fcstYquantiles(:,:,:,ndxT) = res.yquantiles;
    fcstPIT(:,:,ndxT)          = res.pit;
    fcstYdraws_h1(:,:,ndxT)    = res.ydraws_h1;

    if mod(ndxT, 5) == 0 || ndxT == 1 || ndxT == Njumpoffs
        fprintf('  [%s] Origin %3d/%d (%s)\n', datestr(now,'HH:MM:SS'), ...
            ndxT, Njumpoffs, datestr(dates_trans(thisT),'yyyy:mmm'));
    end
end

elapsedTotal = toc(totalTimer);
fprintf('\nDone: %.1f sec (%.1f min, %.1f hrs)\n', ...
    elapsedTotal, elapsedTotal/60, elapsedTotal/3600);

%% ---- 7. Save ----
modelname = 'DPMCSV';
results_dir = fullfile('..', 'results');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end
save(fullfile(results_dir, sprintf('results_%s_forecast.mat', modelname)), ...
    'fcstYrealized', 'fcstYhat', 'fcstYmedian', ...
    'fcstCRPS', 'fcstLogscore', 'fcstJointLogscore', ...
    'fcstYquantiles', 'fcstPIT', 'fcstYdraws_h1', ...
    'setQuantiles', 'Njumpoffs', 'Tjumpoffs', ...
    'MCMCdraws', 'burnin', 'Ndraws', 'fcstNdraws', 'fcstNhorizons', ...
    'N', 'p', 'n', 'labels', 'codes', ...
    'dates_trans', 'modelname', ...
    '-v7.3');
fprintf('Saved ../results/results_%s_forecast.mat\n', modelname);

%% ---- 8. Quick diagnostic ----
fprintf('\n--- Avg CRPS at h=1 ---\n');
avg_crps = mean(fcstCRPS(:,1,:), 3, 'omitnan');
for j = 1:N
    fprintf('  %-28s: %.4f\n', labels{j}, avg_crps(j));
end
fprintf('\n--- Avg Joint Log Score h=1: %.2f ---\n', ...
    mean(fcstJointLogscore,'omitnan'));
