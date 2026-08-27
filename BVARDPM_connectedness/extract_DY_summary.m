% =========================================================================
% extract_DY_summary.m
%
% Extracts a compact summary of the DY results for cross-model
% sanity checking. Reads both result files (DPM-OISV and SV)
% from results/ and writes a single dy_summary.mat file that is
% small enough to upload to Claude for verification.
%
% The summary contains:
%   - MCMC diagnostics (spectral radius, unstable draws, runtime)
%   - DPM parameters (alpha, K) for the DPM model
%   - Total connectedness time series with 16/50/84 quantiles
%   - Per-bank net transmission (averaged over time and draws)
%   - Snapshot GFEVD matrices (averaged over draws) at the 5 crisis dates
%   - Dates, bank labels, and metadata
%
% Usage (run from BVARDPM_connectedness/):
%   >> extract_DY_summary
%
% Output:
%   dy_summary.mat   (approx 2-5 MB, safe to upload)
%
% =========================================================================

clear; clc;

fprintf('========== DY Results Summary Extraction ==========\n\n');

%% ===================== Load results ====================================
results_dir = 'results/';

dpm_file = [results_dir 'results_DPM_OIMSV_DY.mat'];
sv_file  = [results_dir 'results_SV_DY.mat'];

if ~exist(dpm_file, 'file')
    error('Cannot find %s. Run BVAR_DPM_OIMSV_FE_DY.m first.', dpm_file);
end
if ~exist(sv_file, 'file')
    error('Cannot find %s. Run run_SV_fullsample_DY_fast.m first.', sv_file);
end

fprintf('Loading %s ...\n', dpm_file);
D = load(dpm_file, ...
    'store_C_total', 'store_C_from', 'store_C_to', ...
    'store_DPM_params', 'store_spec_radius', ...
    'store_GFEVD_snap', 'store_SV_params', ...
    'dates_subsample', 'snapshot_actual', 'snapshot_target_dates', ...
    't_subsample', 'T_sub', 'labels', 'bank_info', ...
    'n', 'T', 'p', 'k', 'nsims', 'burnin', 'thin', 'nsims_thin', ...
    'elapsed_time', 'n_unstable');

fprintf('Loading %s ...\n', sv_file);
S = load(sv_file, ...
    'store_C_total', 'store_C_from', 'store_C_to', ...
    'store_spec_radius', 'store_GFEVD_snap', 'store_SV_params', ...
    'dates_subsample', 'snapshot_actual', 'snapshot_target_dates', ...
    't_subsample', 'T_sub', 'labels', 'bank_info', ...
    'n', 'T', 'p', 'k', 'nsims', 'burnin', 'thin', 'nsims_thin', ...
    'elapsed_time', 'n_unstable');

% n_rejected may not exist in all result files
try
    tmp = load(sv_file, 'n_rejected');
    S.n_rejected = tmp.n_rejected;
catch
    S.n_rejected = NaN;
end

fprintf('Both files loaded successfully.\n\n');

%% ===================== Consistency checks ==============================
fprintf('---- Consistency checks ----\n');
assert(D.n == S.n, 'n mismatch');
assert(D.T == S.T, 'T mismatch');
assert(D.p == S.p, 'p mismatch');
assert(D.T_sub == S.T_sub, 'T_sub mismatch');
assert(all(D.dates_subsample == S.dates_subsample), 'dates_subsample mismatch');
for ib = 1:D.n
    assert(strcmp(D.labels{ib}, S.labels{ib}), ...
        sprintf('label mismatch at bank %d: %s vs %s', ib, D.labels{ib}, S.labels{ib}));
end
fprintf('  n = %d, T = %d, p = %d, T_sub = %d (matched)\n', D.n, D.T, D.p, D.T_sub);
fprintf('  Bank labels match across both models\n\n');

%% ===================== DPM-OISV summary =================================
fprintf('---- DPM-OISV ----\n');
fprintf('  nsims=%d, burnin=%d, thin=%d, nsims_thin=%d\n', ...
    D.nsims, D.burnin, D.thin, D.nsims_thin);
fprintf('  elapsed: %.1f sec (%.1f min)\n', D.elapsed_time, D.elapsed_time/60);
fprintf('  mean spectral radius (thinned): %.4f\n', mean(D.store_spec_radius));
fprintf('  max  spectral radius (thinned): %.4f\n', max(D.store_spec_radius));
fprintf('  unstable draws (sr >= 1): %d / %d (%.1f%%)\n', ...
    D.n_unstable, D.nsims_thin, 100*D.n_unstable/D.nsims_thin);
fprintf('  mean cluster count K: %.2f\n', mean(D.store_DPM_params(:,2)));
fprintf('  range of K: [%d, %d]\n', ...
    min(D.store_DPM_params(:,2)), max(D.store_DPM_params(:,2)));
fprintf('  mean alpha: %.3f\n', mean(D.store_DPM_params(:,1)));
fprintf('  range of alpha: [%.3f, %.3f]\n', ...
    min(D.store_DPM_params(:,1)), max(D.store_DPM_params(:,1)));

% Total connectedness (median across unthinned draws)
C_med_D = median(D.store_C_total, 2);
C_q16_D = prctile(D.store_C_total, 16, 2);
C_q84_D = prctile(D.store_C_total, 84, 2);
fprintf('\n  Total connectedness C^H_t (median across %d draws):\n', size(D.store_C_total, 2));
fprintf('    sample mean:     %.2f%%\n', mean(C_med_D));
fprintf('    min:             %.2f%% (at %s)\n', ...
    min(C_med_D), datestr(D.dates_subsample(find(C_med_D == min(C_med_D), 1)), 'yyyy-mm-dd'));
fprintf('    max:             %.2f%% (at %s)\n', ...
    max(C_med_D), datestr(D.dates_subsample(find(C_med_D == max(C_med_D), 1)), 'yyyy-mm-dd'));
fprintf('    dynamic range:   %.2f pp\n', max(C_med_D) - min(C_med_D));

%% ===================== SV summary =======================================
fprintf('\n---- SV (AR(1) KSC sampler) ----\n');
fprintf('  nsims=%d, burnin=%d, thin=%d, nsims_thin=%d\n', ...
    S.nsims, S.burnin, S.thin, S.nsims_thin);
fprintf('  elapsed: %.1f sec (%.1f min)\n', S.elapsed_time, S.elapsed_time/60);
fprintf('  mean spectral radius (thinned): %.4f\n', mean(S.store_spec_radius));
fprintf('  max  spectral radius (thinned): %.4f\n', max(S.store_spec_radius));
fprintf('  unstable draws (sr >= 1): %d / %d (%.1f%%)\n', ...
    S.n_unstable, S.nsims_thin, 100*S.n_unstable/S.nsims_thin);
if ~isnan(S.n_rejected)
    fprintf('  stationarity filter rejections: %d (%.1f%% of %d iterations)\n', ...
        S.n_rejected, 100*S.n_rejected/(S.nsims+S.burnin), S.nsims+S.burnin);
end

% SV persistence diagnostic
phi_S = S.store_SV_params(:, 1:S.n);      % nsims_thin x n
omega2_S = S.store_SV_params(:, S.n+1:end);
fprintf('  mean posterior phi (SV persistence): %.4f  [range %.3f to %.3f]\n', ...
    mean(phi_S(:)), min(mean(phi_S,1)), max(mean(phi_S,1)));
fprintf('  mean posterior omega2 (SV innov var): %.4f\n', mean(omega2_S(:)));

% Total connectedness (median)
C_med_S = median(S.store_C_total, 2);
C_q16_S = prctile(S.store_C_total, 16, 2);
C_q84_S = prctile(S.store_C_total, 84, 2);
fprintf('\n  Total connectedness C^H_t (median across %d draws):\n', size(S.store_C_total, 2));
fprintf('    sample mean:     %.2f%%\n', mean(C_med_S));
fprintf('    min:             %.2f%% (at %s)\n', ...
    min(C_med_S), datestr(S.dates_subsample(find(C_med_S == min(C_med_S), 1)), 'yyyy-mm-dd'));
fprintf('    max:             %.2f%% (at %s)\n', ...
    max(C_med_S), datestr(S.dates_subsample(find(C_med_S == max(C_med_S), 1)), 'yyyy-mm-dd'));
fprintf('    dynamic range:   %.2f pp\n', max(C_med_S) - min(C_med_S));

%% ===================== Cross-model comparison ==========================
fprintf('\n---- Cross-model comparison ----\n');

range_D = max(C_med_D) - min(C_med_D);
range_S = max(C_med_S) - min(C_med_S);
fprintf('  Dynamic range: DPM = %.2f pp,  SV = %.2f pp\n', range_D, range_S);
fprintf('    difference  = %+.2f pp\n', range_D - range_S);
fprintf('    ratio DPM/SV = %.2fx\n', range_D / range_S);

% Per-bank net transmission
% store_C_from / store_C_to are (T_sub x n x nsims_thin)
% Average over time (dim 1) and draws (dim 3)
% NOTE: squeeze(mean(mean(...))) produces a 1xn ROW vector in MATLAB.
% Force column vectors with (:) so corr() sees 30 observations of 1 variable.
avg_from_D = squeeze(mean(mean(D.store_C_from, 3), 1));  avg_from_D = avg_from_D(:);
avg_to_D   = squeeze(mean(mean(D.store_C_to,   3), 1));  avg_to_D   = avg_to_D(:);
net_D      = avg_to_D - avg_from_D;

avg_from_S = squeeze(mean(mean(S.store_C_from, 3), 1));  avg_from_S = avg_from_S(:);
avg_to_S   = squeeze(mean(mean(S.store_C_to,   3), 1));  avg_to_S   = avg_to_S(:);
net_S      = avg_to_S - avg_from_S;

fprintf('\n  Per-bank net transmission:\n');
fprintf('    DPM range: [%.2f, %.2f] pp\n', min(net_D), max(net_D));
fprintf('    SV  range: [%.2f, %.2f] pp\n', min(net_S), max(net_S));

% Ranks
[~, ord_D] = sort(net_D, 'descend');
[~, ord_S] = sort(net_S, 'descend');
rk_D = zeros(D.n, 1);  rk_D(ord_D) = 1:D.n;
rk_S = zeros(S.n, 1);  rk_S(ord_S) = 1:S.n;

rho_spear = corr(net_D, net_S, 'Type', 'Spearman');
rho_pear  = corr(net_D, net_S, 'Type', 'Pearson');
fprintf('\n  Rank correlation of net transmission:\n');
fprintf('    Spearman rho: %.4f\n', rho_spear);
fprintf('    Pearson  rho: %.4f\n', rho_pear);

% Top net transmitters under each model
[~, top5_D] = maxk(net_D, 5);
[~, top5_S] = maxk(net_S, 5);
[~, bot5_D] = mink(net_D, 5);
[~, bot5_S] = mink(net_S, 5);

fprintf('\n  Top 5 net transmitters:\n');
fprintf('    DPM: ');
for k = 1:5, fprintf('%s(%+.1f)  ', D.labels{top5_D(k)}, net_D(top5_D(k))); end
fprintf('\n');
fprintf('    SV:  ');
for k = 1:5, fprintf('%s(%+.1f)  ', S.labels{top5_S(k)}, net_S(top5_S(k))); end
fprintf('\n');

fprintf('\n  Bottom 5 net transmitters (top net receivers):\n');
fprintf('    DPM: ');
for k = 1:5, fprintf('%s(%+.1f)  ', D.labels{bot5_D(k)}, net_D(bot5_D(k))); end
fprintf('\n');
fprintf('    SV:  ');
for k = 1:5, fprintf('%s(%+.1f)  ', S.labels{bot5_S(k)}, net_S(bot5_S(k))); end
fprintf('\n');

% Largest rank discrepancies
rank_diff = rk_S - rk_D;
[~, big_disp] = sort(abs(rank_diff), 'descend');
fprintf('\n  Top 10 rank discrepancies (|SV rank - DPM rank|):\n');
fprintf('    %-14s %8s %8s %10s\n', 'ticker', 'DPM rk', 'SV rk', 'diff');
for k = 1:min(10, D.n)
    ib = big_disp(k);
    fprintf('    %-14s %8d %8d %+10d\n', ...
        D.labels{ib}, rk_D(ib), rk_S(ib), rank_diff(ib));
end

%% ===================== Snapshot total connectedness ====================
% Cross-check: total connectedness computed from the avg snapshot GFEVD
% matrix at the 5 crisis dates, compared to the nearest subsampled C_total
fprintf('\n  Total connectedness at crisis snapshots (from avg GFEVD matrix):\n');
fprintf('    %-12s %12s %12s\n', 'Date', 'DPM', 'SV');
N_snap = length(D.snapshot_actual);
snap_C_D = zeros(N_snap, 1);
snap_C_S = zeros(N_snap, 1);
for s = 1:N_snap
    G_D = mean(D.store_GFEVD_snap(:,:,s,:), 4);
    G_S = mean(S.store_GFEVD_snap(:,:,s,:), 4);
    snap_C_D(s) = (sum(G_D(:)) - trace(G_D)) / D.n * 100;
    snap_C_S(s) = (sum(G_S(:)) - trace(G_S)) / S.n * 100;
    fprintf('    %-12s %12.2f %12.2f\n', ...
        D.snapshot_actual{s}, snap_C_D(s), snap_C_S(s));
end

%% ===================== Build compact summary structure ================
fprintf('\n---- Building compact summary file ----\n');

summary = struct();

% Metadata
summary.labels           = D.labels;
summary.bank_info        = D.bank_info;
summary.dates_subsample  = D.dates_subsample;
summary.T_sub            = D.T_sub;
summary.n                = D.n;
summary.T                = D.T;
summary.p                = D.p;

% DPM-OISV metadata and diagnostics
summary.dpm.nsims            = D.nsims;
summary.dpm.burnin           = D.burnin;
summary.dpm.thin             = D.thin;
summary.dpm.nsims_thin       = D.nsims_thin;
summary.dpm.elapsed_time     = D.elapsed_time;
summary.dpm.n_unstable       = D.n_unstable;
summary.dpm.spec_radius_mean = mean(D.store_spec_radius);
summary.dpm.spec_radius_max  = max(D.store_spec_radius);
summary.dpm.K_mean           = mean(D.store_DPM_params(:,2));
summary.dpm.K_range          = [min(D.store_DPM_params(:,2)), max(D.store_DPM_params(:,2))];
summary.dpm.alpha_mean       = mean(D.store_DPM_params(:,1));
summary.dpm.alpha_range      = [min(D.store_DPM_params(:,1)), max(D.store_DPM_params(:,1))];
summary.dpm.store_DPM_params = D.store_DPM_params;
summary.dpm.store_spec_radius = D.store_spec_radius;

% DPM total connectedness time series (3 quantiles, T_sub x 1 each)
summary.dpm.C_total_q16  = C_q16_D;
summary.dpm.C_total_med  = C_med_D;
summary.dpm.C_total_q84  = C_q84_D;

% DPM per-bank net transmission
summary.dpm.net_trans    = net_D;
summary.dpm.from_trans   = avg_from_D;
summary.dpm.to_trans     = avg_to_D;
summary.dpm.rank         = rk_D;

% DPM snapshot GFEVD (avg over draws)
summary.dpm.GFEVD_snap_avg = zeros(D.n, D.n, N_snap);
for s = 1:N_snap
    summary.dpm.GFEVD_snap_avg(:,:,s) = mean(D.store_GFEVD_snap(:,:,s,:), 4);
end
summary.dpm.snapshot_actual = D.snapshot_actual;
summary.dpm.snap_C_total    = snap_C_D;

% SV metadata and diagnostics
summary.sv.nsims            = S.nsims;
summary.sv.burnin           = S.burnin;
summary.sv.thin             = S.thin;
summary.sv.nsims_thin       = S.nsims_thin;
summary.sv.elapsed_time     = S.elapsed_time;
summary.sv.n_unstable       = S.n_unstable;
summary.sv.spec_radius_mean = mean(S.store_spec_radius);
summary.sv.spec_radius_max  = max(S.store_spec_radius);
summary.sv.phi_mean         = mean(phi_S, 1);        % per-variable posterior mean
summary.sv.omega2_mean      = mean(omega2_S, 1);
summary.sv.store_spec_radius = S.store_spec_radius;

% SV total connectedness time series
summary.sv.C_total_q16  = C_q16_S;
summary.sv.C_total_med  = C_med_S;
summary.sv.C_total_q84  = C_q84_S;

% SV per-bank net transmission
summary.sv.net_trans    = net_S;
summary.sv.from_trans   = avg_from_S;
summary.sv.to_trans     = avg_to_S;
summary.sv.rank         = rk_S;

% SV snapshot GFEVD (avg over draws)
summary.sv.GFEVD_snap_avg = zeros(S.n, S.n, N_snap);
for s = 1:N_snap
    summary.sv.GFEVD_snap_avg(:,:,s) = mean(S.store_GFEVD_snap(:,:,s,:), 4);
end
summary.sv.snapshot_actual = S.snapshot_actual;
summary.sv.snap_C_total    = snap_C_S;

% Cross-model summary
summary.cross.range_D       = range_D;
summary.cross.range_S       = range_S;
summary.cross.range_ratio   = range_D / range_S;
summary.cross.spearman_rho  = rho_spear;
summary.cross.pearson_rho   = rho_pear;
summary.cross.rank_diff     = rank_diff;

% Save
save('dy_summary.mat', '-struct', 'summary', '-v7.3');

f = dir('dy_summary.mat');
fprintf('  Saved dy_summary.mat (%.1f MB)\n', f.bytes/1e6);

fprintf('\n========== Done ==========\n');
fprintf('Please upload dy_summary.mat to Claude for verification.\n');
