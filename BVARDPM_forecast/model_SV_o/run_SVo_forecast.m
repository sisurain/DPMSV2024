% =========================================================================
% run_SVo_forecast.m
%
% Expanding-window forecasting for SVo (SV + common scalar outlier) BVAR.
%
% MODEL (CCMM 2022, REStat supplement):
%   v_t = obar_t * A^{-1} Lambda_t^{0.5} e_t
%   obar_t: SCALAR common outlier in {1,...,20}
%   P(obar_t=1) = 1-pi, P(obar_t=k>1) = pi/19
%   h_t = h_{t-1} + eta_t ~ N(0,Phi)
%   pi ~ Beta(alpha,beta)
%
% MCMC: nsim=1000, burnin=200, Ndraws=10 -> 10,000 total draws
%
% REQUIRED FILES:
%   fcst_SVo_oneorigin.m, CTA.m, StochVolKSCcorrsqrt.m,
%   rwnoisePrecisionBasedSampler.m, getKSC7values.m,
%   crpsDraws.m, vec.m
%   ../data/FREDMD42_forecast_data.mat
% =========================================================================

clear; clc; close all;

%% ---- 1. Settings ----
MCMCdraws          = 1e3;
Ndraws             = 10;
fcstNhorizons      = 24;
setQuantiles       = [5, 16, 50, 84, 95];
Nquantiles         = length(setQuantiles);

doTightPrior       = false;
doRobustPrior      = true;
check_stationarity = 0;

if doTightPrior
    theta = [0.05 0.5 100 2];
else
    theta = [0.1  0.5 100 2];
end

burnin     = 2 * ceil(0.1 * MCMCdraws);
MCMCreps   = MCMCdraws + burnin;
fcstNdraws = MCMCdraws * Ndraws;

%% ---- 2. Load data ----
load('../data/FREDMD42_forecast_data.mat');
fprintf('Loaded FREDMD42_forecast_data.mat: n=%d, T_trans=%d, p=%d\n', n, T_trans, p);

N     = n;
data  = data_trans;
Tdata = T_trans;
Nlogtwopi = N * log(2*pi);
mnPrior   = minnesotaPriorMean;

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

%% ---- 5. Parallel pool (12 workers) ----
poolobj = gcp('nocreate');
if isempty(poolobj)
    poolobj = parpool('local', 12);
end
fprintf('Using %d workers\n\n', poolobj.NumWorkers);

%% ---- 6. Main loop ----
fprintf('Starting SVo forecast: %d origins x %d MCMC iters x %d draws\n', ...
    Njumpoffs, MCMCreps, Ndraws);
totalTimer = tic;

parfor ndxT = 1 : Njumpoffs

    thisT = Tjumpoffs(ndxT);

    res = fcst_SVo_oneorigin(thisT, Tdata, data, dates_trans, ...
        N, p, MCMCdraws, burnin, MCMCreps, Ndraws, fcstNdraws, fcstNhorizons, ...
        theta, mnPrior, doRobustPrior, check_stationarity, setQuantiles, Nlogtwopi);

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
modelname = 'SVo';
save(sprintf('results_%s_forecast.mat', modelname), ...
    'fcstYrealized', 'fcstYhat', 'fcstYmedian', ...
    'fcstCRPS', 'fcstLogscore', 'fcstJointLogscore', ...
    'fcstYquantiles', 'fcstPIT', 'fcstYdraws_h1', ...
    'setQuantiles', 'Njumpoffs', 'Tjumpoffs', ...
    'MCMCdraws', 'burnin', 'Ndraws', 'fcstNdraws', 'fcstNhorizons', ...
    'N', 'p', 'n', 'labels', 'codes', ...
    'dates_trans', 'modelname', ...
    '-v7.3');
fprintf('Saved results_%s_forecast.mat\n', modelname);

%% ---- 8. Quick diagnostic ----
fprintf('\n--- Avg CRPS at h=1 ---\n');
avg_crps = mean(fcstCRPS(:,1,:), 3, 'omitnan');
for j = 1:N
    fprintf('  %-28s: %.4f\n', labels{j}, avg_crps(j));
end
fprintf('\n--- Avg Joint Log Score h=1: %.2f ---\n', ...
    mean(fcstJointLogscore,'omitnan'));
