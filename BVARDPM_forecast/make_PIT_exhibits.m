% =========================================================================
% make_PIT_exhibits.m
%
% R2 revision items R5-M5 (PIT uniformity statistic) and the recommended
% PIT-ECDF appendix figure.
%
% (a) Computes, for each model and variable, the uniformity statistic
%     proposed by referee 5: the sum over the 10 histogram bins of the
%     absolute deviation of the bin density from 1, using the h = 1 PITs
%     over all 250 forecast origins. Prints the table and saves it as CSV.
% (b) Plots the empirical CDFs of the h = 1 PITs against the 45-degree
%     line (uniformity), one panel per variable, all four models overlaid
%     -- the appendix figure responding to the referee's aside.
%
% Conventions match postprocess_forecast_results.m: h = 1, Nbins = 10,
% models SVo, SV-t, DPM-CSV, DPM-OISV; variables IP, PCE Prices,
% Unemployment Rate.
% =========================================================================

clear; clc; close all;

%% ---- Configuration ----
results_dir = 'results/';
output_dir  = 'pit_exhibits/';
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

modelNames  = {'SV_o', 'SVt', 'DPMCSV', 'DPMOISV'};
modelLabels = {'SVo', 'SV-t', 'DPM-CSV', 'DPM-OISV'};
Nmodels = length(modelNames);
Nbins   = 10;

var_search = {
    {'INDPRO', 'Industrial', 'IP ('}
    {'PCEPI', 'PCE Price', 'PCE price'}
    {'UNRATE', 'Unemployment', 'Unemp'}
};

%% ---- Load results ----
R = cell(Nmodels, 1);
for mm = 1:Nmodels
    fname = sprintf('%sresults_%s_forecast.mat', results_dir, modelNames{mm});
    if ~exist(fname, 'file'), error('Missing: %s', fname); end
    R{mm} = load(fname, 'fcstPIT', 'labels');
end
labels = R{1}.labels;

var_idx = NaN(1, length(var_search));
for vv = 1:length(var_search)
    for pp = 1:length(var_search{vv})
        idx_found = find(contains(labels, var_search{vv}{pp}, 'IgnoreCase', true));
        if ~isempty(idx_found), var_idx(vv) = idx_found(1); break; end
    end
end
var_labels = labels(var_idx);

%% ---- (a) Uniformity statistic: sum_b |density_b - 1| at h = 1 ---------
U = NaN(length(var_idx), Nmodels);
for vv = 1:length(var_idx)
    jj = var_idx(vv);
    for mm = 1:Nmodels
        pit = squeeze(R{mm}.fcstPIT(jj, 1, :));   % (var, h, origin), h = 1
        pit = pit(~isnan(pit));
        cnt = histcounts(pit, linspace(0, 1, Nbins+1));
        dens = cnt / numel(pit) * Nbins;
        U(vv, mm) = sum(abs(dens - 1));
    end
end

fprintf('PIT uniformity statistic (sum_b |bin density - 1|), h = 1:\n');
fprintf('%-22s', 'Variable');
fprintf('%12s', modelLabels{:}); fprintf('\n');
for vv = 1:length(var_idx)
    fprintf('%-22s', var_labels{vv});
    fprintf('%12.3f', U(vv, :)); fprintf('\n');
end

fid = fopen(fullfile(output_dir, 'pit_uniformity_stats.csv'), 'w');
fprintf(fid, 'variable'); fprintf(fid, ',%s', modelLabels{:}); fprintf(fid, '\n');
for vv = 1:length(var_idx)
    fprintf(fid, '%s', var_labels{vv});
    fprintf(fid, ',%.4f', U(vv, :)); fprintf(fid, '\n');
end
fclose(fid);

%% ---- (b) PIT ECDF appendix figure ------------------------------------
mcolors = lines(Nmodels);
fig = figure('Name', 'PIT ECDFs (h=1)', 'Position', [100 100 1400 420]);
for vv = 1:length(var_idx)
    jj = var_idx(vv);
    subplot(1, length(var_idx), vv); hold on;
    plot([0 1], [0 1], 'k--', 'LineWidth', 1);        % uniform reference
    hs = gobjects(Nmodels, 1);
    for mm = 1:Nmodels
        pit = squeeze(R{mm}.fcstPIT(jj, 1, :));
        pit = sort(pit(~isnan(pit)));
        n = numel(pit);
        hs(mm) = stairs([0; pit; 1], [0; (1:n)'/n; 1], ...
            'Color', mcolors(mm, :), 'LineWidth', 1.3);
    end
    hold off; axis([0 1 0 1]); axis square; box off;
    xlabel('PIT'); ylabel('Empirical CDF');
    title(var_labels{vv}, 'Interpreter', 'none', 'FontSize', 10);
    if vv == 1
        legend(hs, modelLabels, 'Location', 'southeast', 'Box', 'off');
    end
end
saveas(fig, fullfile(output_dir, 'pit_ecdf_h1.png'));
savefig(fig, fullfile(output_dir, 'pit_ecdf_h1.fig'));
fprintf('\nSaved PIT ECDF figure and statistics to %s\n', output_dir);
