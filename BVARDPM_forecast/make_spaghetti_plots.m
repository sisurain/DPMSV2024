% =========================================================================
% make_spaghetti_plots.m
%
% Standalone script to generate the forecast-path (spaghetti) figures:
%   for each forecast origin, the full h = 1..12 predictive path is
%   plotted on calendar dates T+1..T+12, overlaying all paths, with
%   the realized series in solid black on top.
%
% Two views are produced for each variable: the full evaluation sample
% and a COVID-period zoom (file suffix _covid). The forecast-path
% figures in the paper are the COVID-zoom versions for PCE inflation
% and the unemployment rate.
%
% Saves both .png and .fig versions into spaghetti_figures/.
%
% REQUIRED (in results_dir):
%   results_SV_o_forecast.mat, results_SVt_forecast.mat,
%   results_DPMOISV_forecast.mat
% =========================================================================

clear; clc; close all;

%% ---- Configuration ----
results_dir = 'results/';
output_dir  = 'spaghetti_figures/';
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

modelNames  = {'SV_o', 'SVt', 'DPMOISV'};
modelLabels = {'SVo', 'SV-t', 'DPM-OISV'};
Nmodels = length(modelNames);

H_plot = 12;

var_search = {
    {'INDPRO', 'Industrial', 'IP ('}
    {'PCEPI', 'PCE Price', 'PCE price'}
    {'UNRATE', 'Unemployment', 'Unemp'}
};
% Optional per-variable cap on the common upper y-limit (Inf = no cap).
% Unemployment: one DPM-OISV median path from the May 2020 jumpoff peaks
% at 53.4 and would compress all panels; capped at 30 with a caption
% disclosure (see tex).
var_ycap = [Inf, Inf, 30];

covid_start = datenum(2020, 3, 1);
covid_end   = datenum(2021, 12, 31);

views = {
    struct('name','full',  'xlim',[],                                    'suffix','')
    struct('name','covid', 'xlim',[datenum(2019,1,1) datenum(2023,6,1)], 'suffix','_covid')
};

%% ---- Load results ----
fprintf('Loading results...\n');
R = cell(Nmodels, 1);
for mm = 1:Nmodels
    fname = sprintf('%sresults_%s_forecast.mat', results_dir, modelNames{mm});
    if ~exist(fname, 'file'), error('Missing: %s', fname); end
    R{mm} = load(fname);
    fprintf('  Loaded %s\n', modelNames{mm});
end

N         = R{1}.N;
labels    = R{1}.labels;
dates     = R{1}.dates_trans;
Tjumpoffs = R{1}.Tjumpoffs;
Njumpoffs = R{1}.Njumpoffs;

%% ---- Match variable indices ----
var_idx = NaN(1, length(var_search));
for vv = 1:length(var_search)
    for pp = 1:length(var_search{vv})
        idx_found = find(contains(labels, var_search{vv}{pp}, 'IgnoreCase', true));
        if ~isempty(idx_found)
            var_idx(vv) = idx_found(1);
            fprintf('  Variable %d: "%s" (index %d)\n', vv, labels{var_idx(vv)}, var_idx(vv));
            break;
        end
    end
end

colors = lines(Nmodels);

%% ---- Precompute realized series span (model-independent) ----
first_date_idx = Tjumpoffs(1) + 1;
last_date_idx  = min(Tjumpoffs(end) + H_plot, length(dates));
realized_dates = dates(first_date_idx:last_date_idx);

%% ---- Generate spaghetti plots ----
fprintf('\n=== Generating spaghetti plots (common y-axis per variable) ===\n');

for vw = 1:length(views)
    view_cfg = views{vw};
    fprintf('\n-- View: %s --\n', view_cfg.name);

    for vv = 1:length(var_idx)
        if isnan(var_idx(vv)), continue; end
        jj = var_idx(vv);

        % ---- PASS 1 (R5-M4): collect all values across ALL models ------
        all_vals = [];
        realized_vals = NaN(length(realized_dates), 1);
        for mm = 1:Nmodels
            for tt = 1:Njumpoffs
                end_idx = Tjumpoffs(tt) + H_plot;
                if end_idx > length(dates), continue; end
                path_dates = dates(Tjumpoffs(tt) + (1:H_plot));
                path_vals  = squeeze(R{mm}.fcstYmedian(jj, 1:H_plot, tt));
                if isempty(view_cfg.xlim)
                    all_vals = [all_vals; path_vals(:)]; %#ok<AGROW>
                else
                    in_win = path_dates >= view_cfg.xlim(1) & path_dates <= view_cfg.xlim(2);
                    sel = path_vals(in_win(:));
                    all_vals = [all_vals; sel(:)]; %#ok<AGROW>
                end
                % realized (fill once; identical across models)
                if mm == 1
                    for hh = 1:H_plot
                        idx_date = Tjumpoffs(tt) + hh - first_date_idx + 1;
                        if idx_date >= 1 && idx_date <= length(realized_vals)
                            val = R{mm}.fcstYrealized(jj, hh, tt);
                            if ~isnan(val), realized_vals(idx_date) = val; end
                        end
                    end
                end
            end
        end
        if isempty(view_cfg.xlim)
            vis_real = realized_vals;
        else
            in_win   = realized_dates >= view_cfg.xlim(1) & realized_dates <= view_cfg.xlim(2);
            vis_real = realized_vals(in_win);
        end
        visible = [all_vals; vis_real];
        visible = visible(~isnan(visible));
        pad = 0.05 * (max(visible) - min(visible) + eps);
        common_ylim = [min(visible)-pad, max(visible)+pad];
        common_ylim(2) = min(common_ylim(2), var_ycap(vv));   % optional cap

        % ---- PASS 2: draw all panels with the common y-limits ----------
        fig = figure('Name', sprintf('Spaghetti: %s (%s)', labels{jj}, view_cfg.name), ...
            'Position', [100 100 1400 320*Nmodels]);

        for mm = 1:Nmodels
            subplot(Nmodels, 1, mm);
            hold on;
            for tt = 1:Njumpoffs
                end_idx = Tjumpoffs(tt) + H_plot;
                if end_idx > length(dates), continue; end
                path_dates = dates(Tjumpoffs(tt) + (1:H_plot));
                path_vals  = squeeze(R{mm}.fcstYmedian(jj, 1:H_plot, tt));
                plot(path_dates, path_vals, '-', 'Color', [colors(mm,:) 0.35], ...
                     'LineWidth', 0.6);
            end
            plot(realized_dates, realized_vals, 'k-', 'LineWidth', 1.5);
            hold off;

            if ~isempty(view_cfg.xlim), xlim(view_cfg.xlim); end
            ylim(common_ylim);                        % R5-M4: same scale everywhere

            if isempty(view_cfg.xlim)
                datetick('x','yyyy','keeplimits');
            else
                datetick('x','mmmyy','keeplimits');
            end
            ylabel(labels{jj}, 'Interpreter','none');
            title(sprintf('%s: %s (h=1..%d median paths)', modelLabels{mm}, labels{jj}, H_plot), ...
                  'FontSize', 10);
            box off;

            patch([covid_start covid_end covid_end covid_start], ...
                  [common_ylim(1) common_ylim(1) common_ylim(2) common_ylim(2)], 'r', ...
                  'FaceAlpha', 0.05, 'EdgeColor', 'none');
        end

        safe_label = regexprep(labels{jj}, '[^a-zA-Z0-9]', '_');
        fname_base = sprintf('spaghetti_%s_h%d%s', safe_label, H_plot, view_cfg.suffix);
        exportgraphics(fig, [output_dir fname_base '.png'], 'Resolution', 200);
        savefig(fig, [output_dir fname_base '.fig']);
        fprintf('  Saved %s (.png and .fig)\n', fname_base);
    end
end

fprintf('\nDone. Output in: %s\n', output_dir);
