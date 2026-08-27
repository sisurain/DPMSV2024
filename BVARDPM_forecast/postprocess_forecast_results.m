% =========================================================================
% postprocess_forecast_results.m
%
% Master post-processing script for the forecasting exercise.
% Loads the results_XX_forecast.mat files and produces:
%
%   1. RMSE ratio tables by regime (pre-COVID / COVID / post-COVID / full)
%   2. CRPS ratio tables by regime
%   3. DM-West significance tests (*** p<.001, ** p<.01, * p<.05)
%   4. One-step-ahead joint log predictive scores (table + time series)
%   5. PIT histograms (h=1, selected models)
%   6. LaTeX table output
%
% The forecast-path (spaghetti) figures in the paper are produced by
% the standalone script make_spaghetti_plots.m.
%
% REGIME DEFINITIONS:
%   Pre-COVID:  forecast origins up to 2020:Feb
%   COVID:      2020:Mar to 2021:Dec
%   Post-COVID: 2022:Jan onward
%
% BENCHMARK: DPM-OISV (numerator in all ratios; values < 1 favor DPM-OISV)
%
% HORIZONS REPORTED: h = 1, 6, 12
%
% REQUIRED: results_SV_o_forecast.mat, results_SVO_forecast.mat,
%   results_SVt_forecast.mat, results_DPMCSV_forecast.mat,
%   results_DPMOISV_forecast.mat
% =========================================================================

clear; clc; close all;

%% ---- 0. Configuration ----
results_dir = 'results/';
output_dir  = 'tables_figures/';
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

% Log all console output to a text file
diary([output_dir 'console_output.txt']);
diary on;
fprintf('=== Forecast Post-Processing ===\n');
fprintf('Date: %s\n\n', datestr(now));

% DPM-OISV is the benchmark (numerator in all ratios)
% Comparison models: SVo, SVO, SV-t, DPM-CSV
% Ratios < 1 mean DPM-OISV outperforms (lower RMSE/CRPS)
modelNames  = {'SV_o','SVO','SVt','DPMCSV','DPMOISV'};
modelLabels = {'SVo','SVO','SV-t','DPM-CSV','DPM-OISV'};
Nmodels = length(modelNames);

benchmark_idx = 5;  % DPM-OISV is the benchmark (numerator)
compare_idx   = [1, 2, 3, 4];  % SVo, SVO, SV-t, DPM-CSV
horizons_report = [1, 6, 12];  % horizons to report (as in the paper)
Nhorizons_report = length(horizons_report);

%% ---- 1. Load all results ----
fprintf('Loading results...\n');
R = cell(Nmodels, 1);
for mm = 1:Nmodels
    fname = sprintf('%sresults_%s_forecast.mat', results_dir, modelNames{mm});
    if ~exist(fname, 'file')
        warning('Missing: %s', fname);
        continue;
    end
    R{mm} = load(fname);
    fprintf('  Loaded %s: %d origins\n', modelNames{mm}, R{mm}.Njumpoffs);
end

% Extract common variables from benchmark
N       = R{benchmark_idx}.N;
p       = R{benchmark_idx}.p;
labels  = R{benchmark_idx}.labels;
dates   = R{benchmark_idx}.dates_trans;
Tjumpoffs = R{benchmark_idx}.Tjumpoffs;
Njumpoffs = R{benchmark_idx}.Njumpoffs;

% Print variable labels
fprintf('\nVariable labels (%d variables):\n', N);
for j = 1:N
    fprintf('  %2d: "%s"\n', j, labels{j});
end

%% ---- 2. Define regimes ----
origin_dates = dates(Tjumpoffs);

covid_start = datenum(2020, 3, 1);
covid_end   = datenum(2021, 12, 31);
postcovid_start = datenum(2022, 1, 1);

idx_preCOVID  = origin_dates < covid_start;
idx_COVID     = (origin_dates >= covid_start) & (origin_dates <= covid_end);
idx_postCOVID = origin_dates >= postcovid_start;

regime_names = {'Pre-COVID', 'COVID', 'Post-COVID', 'Full Sample'};
regime_idx   = {idx_preCOVID, idx_COVID, idx_postCOVID, true(Njumpoffs,1)};
Nregimes = length(regime_names);

fprintf('\nRegime breakdown:\n');
fprintf('  Pre-COVID:  %d origins (%s to %s)\n', sum(idx_preCOVID), ...
    datestr(min(origin_dates(idx_preCOVID)),'yyyy:mmm'), ...
    datestr(max(origin_dates(idx_preCOVID)),'yyyy:mmm'));
fprintf('  COVID:      %d origins (%s to %s)\n', sum(idx_COVID), ...
    datestr(min(origin_dates(idx_COVID)),'yyyy:mmm'), ...
    datestr(max(origin_dates(idx_COVID)),'yyyy:mmm'));
fprintf('  Post-COVID: %d origins (%s to %s)\n', sum(idx_postCOVID), ...
    datestr(min(origin_dates(idx_postCOVID)),'yyyy:mmm'), ...
    datestr(max(origin_dates(idx_postCOVID)),'yyyy:mmm'));

%% ---- 3. Compute RMSE and CRPS ----
% Squared errors and CRPS: N x H x Njumpoffs
SE = cell(Nmodels, 1);
CRPS_all = cell(Nmodels, 1);

for mm = 1:Nmodels
    if isempty(R{mm}), continue; end
    yreal = R{mm}.fcstYrealized;  % N x H x Njumpoffs
    yhat  = R{mm}.fcstYhat;
    SE{mm}   = (yhat - yreal).^2;
    CRPS_all{mm} = R{mm}.fcstCRPS;
end

%% ---- 4. RMSE/CRPS Ratio Tables + DM Tests ----
fprintf('\n=== Computing RMSE/CRPS ratios and DM tests ===\n');

for rr = 1:Nregimes
    ridx = regime_idx{rr};
    if sum(ridx) == 0, continue; end

    fprintf('\n--- %s (%d origins) ---\n', regime_names{rr}, sum(ridx));

    for metric_type = 1:2
        if metric_type == 1
            metric_name = 'RMSE';
        else
            metric_name = 'CRPS';
        end

        fprintf('\n  %s ratios (%s / other): <1 means %s wins\n', ...
            metric_name, modelLabels{benchmark_idx}, modelLabels{benchmark_idx});
        fprintf('  %-28s', 'Variable');
        for hh = horizons_report
            for ci = compare_idx
                fprintf('%12s h=%d', modelLabels{ci}, hh);
            end
        end
        fprintf('\n');

        for j = 1:N
            fprintf('  %-28s', labels{j});
            for hh = horizons_report
                for ci = compare_idx
                    mm = ci;
                    if isempty(R{mm}), fprintf('%12s', 'N/A'); continue; end

                    if metric_type == 1
                        % RMSE ratio = DPM-OISV / other
                        se_bench = squeeze(SE{benchmark_idx}(j, hh, ridx));
                        se_model = squeeze(SE{mm}(j, hh, ridx));
                        valid = ~isnan(se_bench) & ~isnan(se_model);
                        if sum(valid) < 5
                            fprintf('%12s', 'N/A'); continue;
                        end
                        ratio = sqrt(mean(se_bench(valid))) / sqrt(mean(se_model(valid)));

                        % DM test: d = loss_other - loss_DPMOISV
                        % d > 0 means DPM-OISV has lower loss (wins)
                        d = se_model(valid) - se_bench(valid);
                    else
                        % CRPS ratio = DPM-OISV / other
                        crps_bench = squeeze(CRPS_all{benchmark_idx}(j, hh, ridx));
                        crps_model = squeeze(CRPS_all{mm}(j, hh, ridx));
                        valid = ~isnan(crps_bench) & ~isnan(crps_model);
                        if sum(valid) < 5
                            fprintf('%12s', 'N/A'); continue;
                        end
                        ratio = mean(crps_bench(valid)) / mean(crps_model(valid));

                        % DM test: d = loss_other - loss_DPMOISV
                        d = crps_model(valid) - crps_bench(valid);
                    end

                    % DM-West test with Newey-West HAC
                    [~, pval] = dm_west_test(d, hh);
                    stars = significance_stars(pval, ratio);
                    fprintf('%8.2f%-4s', ratio, stars);
                end
            end
            fprintf('\n');
        end
    end
end

%% ---- 5. Joint Log Predictive Scores ----
fprintf('\n=== Joint Log Predictive Scores (h=1) ===\n');

JLS = NaN(Njumpoffs, Nmodels);
for mm = 1:Nmodels
    if isempty(R{mm}), continue; end
    JLS(:, mm) = R{mm}.fcstJointLogscore(:);
end

% Table by regime
fprintf('\n%-20s', 'Period');
for mm = 1:Nmodels
    fprintf('%12s', modelLabels{mm});
end
fprintf('\n');

for rr = 1:Nregimes
    ridx = regime_idx{rr};
    fprintf('%-20s', regime_names{rr});
    for mm = 1:Nmodels
        vals = JLS(ridx, mm);
        fprintf('%12.2f', mean(vals, 'omitnan'));
    end
    fprintf('\n');
end

% Time series plot
figure('Name','Joint Log Scores','Position',[100 100 1200 400]);

% Plot without extreme outliers
jls_plot = JLS;
for mm = 1:Nmodels
    q01 = prctile(JLS(:,mm), 1);
    jls_plot(JLS(:,mm) < q01, mm) = NaN;
end

hold on;
colors = lines(Nmodels);
for mm = 1:Nmodels  % all 5 models: SVo, SVO, SV-t, DPM-CSV, DPM-OISV
    plot(origin_dates, jls_plot(:,mm), '-', 'Color', colors(mm,:), ...
        'LineWidth', 1, 'DisplayName', modelLabels{mm});
end
hold off;
datetick('x','yyyy','keeplimits');
ylabel('Joint Log Predictive Score');
title('One-Step-Ahead Joint Log Scores (1st/99th percentile trimmed)');
legend('Location','southwest');
grid on; box off;
saveas(gcf, [output_dir 'joint_logscores.png']);

%% ---- 6. PIT Histograms ----
fprintf('\n=== Generating PIT Histograms ===\n');

pit_models = [1, 3, 4, 5];  % SVo, SV-t, DPM-CSV, DPM-OISV

% Search for IP, PCE, Unemployment using multiple possible name patterns
pit_search = {
    {'INDPRO', 'Industrial', 'IP ('}       % Industrial Production
    {'PCEPI', 'PCE Price', 'PCE price'}    % PCE Price Index
    {'UNRATE', 'Unemployment', 'Unemp'}    % Unemployment Rate
};
pit_var_idx = NaN(1, length(pit_search));
for vv = 1:length(pit_search)
    for pp = 1:length(pit_search{vv})
        idx_found = find(contains(labels, pit_search{vv}{pp}, 'IgnoreCase', true));
        if ~isempty(idx_found)
            pit_var_idx(vv) = idx_found(1);
            break;
        end
    end
end

% Diagnostic: show what was found
fprintf('PIT variable matching:\n');
for vv = 1:length(pit_search)
    if isnan(pit_var_idx(vv))
        fprintf('  %s: NOT FOUND\n', pit_search{vv}{1});
    else
        fprintf('  %s: matched -> labels{%d} = "%s"\n', ...
            pit_search{vv}{1}, pit_var_idx(vv), labels{pit_var_idx(vv)});
    end
end

% Fallback: if nothing matched, use variables 3, 5, 10 (or whatever exists)
if all(isnan(pit_var_idx))
    fprintf('  WARNING: No matches found. Using variables 3, 5, 10 as fallback.\n');
    pit_var_idx = [min(3,N), min(5,N), min(10,N)];
end

figure('Name','PIT Histograms','Position',[100 100 1400 800]);
Nbins = 10;
panel = 0;
n_pit_vars = sum(~isnan(pit_var_idx));
for vv = 1:length(pit_var_idx)
    if isnan(pit_var_idx(vv)), continue; end
    jj = pit_var_idx(vv);
    for mi = 1:length(pit_models)
        mm = pit_models(mi);
        if isempty(R{mm}), continue; end
        panel = panel + 1;
        subplot(n_pit_vars, length(pit_models), panel);

        pit_vals = squeeze(R{mm}.fcstPIT(jj, 1, :));  % h=1
        pit_vals = pit_vals(~isnan(pit_vals));
        
        fprintf('  PIT: model=%s, var=%s, n_valid=%d\n', ...
            modelLabels{mm}, labels{jj}, length(pit_vals));

        if isempty(pit_vals)
            text(0.5, 0.5, 'No data', 'HorizontalAlignment', 'center');
            title(sprintf('%s: %s', modelLabels{mm}, labels{jj}), 'FontSize', 8);
            continue;
        end

        histogram(pit_vals, linspace(0,1,Nbins+1), 'Normalization','pdf', ...
            'FaceColor', colors(mm,:), 'FaceAlpha', 0.7);
        hold on;
        yline(1, 'r--', 'LineWidth', 1.5);  % uniform reference
        hold off;
        xlim([0 1]); ylim([0 3]);
        title(sprintf('%s: %s', modelLabels{mm}, labels{jj}), 'FontSize', 8);
        if mi == 1, ylabel('Density'); end
        if vv == n_pit_vars, xlabel('PIT'); end
    end
end
sgtitle('PIT Histograms (h=1)');
saveas(gcf, [output_dir 'pit_histograms.png']);

%% ---- 7. Export LaTeX Tables ----
fprintf('\n=== Exporting LaTeX Tables ===\n');

for rr = 1:Nregimes
    ridx = regime_idx{rr};
    if sum(ridx) < 5, continue; end

    for metric_type = 1:2
        if metric_type == 1
            metric_name = 'RMSE';
        else
            metric_name = 'CRPS';
        end

        fname_tex = sprintf('%s%s_%s.tex', output_dir, metric_name, ...
            strrep(regime_names{rr}, ' ', '_'));
        fid = fopen(fname_tex, 'w');

        % Header
        Ncompare = length(compare_idx);
        Ncols = 1 + Ncompare * Nhorizons_report;
        fprintf(fid, '\\begin{table}[H]\n\\centering\n');
        fprintf(fid, '\\caption{%s Ratios, %s (DPM-OISV in numerator)}\n', ...
            metric_name, regime_names{rr});
        fprintf(fid, '\\scalebox{0.55}{\n');
        fprintf(fid, '\\begin{tabular}{|l|');
        for ci = 1:Ncompare
            fprintf(fid, '%s', repmat('c', 1, Nhorizons_report));
            if ci < Ncompare, fprintf(fid, '|'); end
        end
        fprintf(fid, '|}\n\\hline\n');

        % Model headers
        fprintf(fid, '\\multirow{2}{*}{Variable}');
        for ci = 1:Ncompare
            mm = compare_idx(ci);
            sep = iif(ci < Ncompare, '|', '');
            fprintf(fid, ' & \\multicolumn{%d}{c%s}{%s}', ...
                Nhorizons_report, sep, modelLabels{mm});
        end
        fprintf(fid, ' \\\\ \\cline{2-%d}\n', Ncols);

        % Horizon headers
        fprintf(fid, ' ');
        for ci = 1:Ncompare
            for hh = horizons_report
                fprintf(fid, ' & %d', hh);
            end
        end
        fprintf(fid, ' \\\\ \\hline\n');

        % Data rows
        for j = 1:N
            fprintf(fid, '%s', labels{j});
            for ci = 1:Ncompare
                mm = compare_idx(ci);
                for hh = horizons_report
                    if isempty(R{mm})
                        fprintf(fid, ' & --'); continue;
                    end

                    if metric_type == 1
                        se_b = squeeze(SE{benchmark_idx}(j,hh,ridx));
                        se_m = squeeze(SE{mm}(j,hh,ridx));
                        valid = ~isnan(se_b) & ~isnan(se_m);
                        if sum(valid) < 5
                            fprintf(fid, ' & --'); continue;
                        end
                        ratio = sqrt(mean(se_b(valid)))/sqrt(mean(se_m(valid)));
                        d = se_m(valid) - se_b(valid);
                    else
                        cr_b = squeeze(CRPS_all{benchmark_idx}(j,hh,ridx));
                        cr_m = squeeze(CRPS_all{mm}(j,hh,ridx));
                        valid = ~isnan(cr_b) & ~isnan(cr_m);
                        if sum(valid) < 5
                            fprintf(fid, ' & --'); continue;
                        end
                        ratio = mean(cr_b(valid))/mean(cr_m(valid));
                        d = cr_m(valid) - cr_b(valid);
                    end

                    [~, pval] = dm_west_test(d, hh);
                    stars = significance_stars_tex(pval, ratio);
                    fprintf(fid, ' & %.2f%s', ratio, stars);
                end
            end
            fprintf(fid, ' \\\\\n');
        end

        fprintf(fid, '\\hline\n\\end{tabular}}\n');
        fprintf(fid, '\\vspace{1ex}\n\n');
        fprintf(fid, '\\raggedright \\footnotesize{DPM-OISV values are in the numerator. ');
        fprintf(fid, 'Values below 1 indicate DPM-OISV improvement. ');
        fprintf(fid, 'Significance (DM-West with Newey-West HAC): ');
        fprintf(fid, '\\textsuperscript{***} $p < 0.001$, ');
        fprintf(fid, '\\textsuperscript{**} $0.001 \\leq p < 0.01$, ');
        fprintf(fid, '\\textsuperscript{*} $0.01 \\leq p < 0.05$.}\n');
        fprintf(fid, '\\end{table}\n');
        fclose(fid);
        fprintf('  Wrote %s\n', fname_tex);
    end
end

% Joint log score table
fname_jls = [output_dir 'joint_logscores.tex'];
fid = fopen(fname_jls, 'w');
fprintf(fid, '\\begin{table}[H]\n\\centering\n');
fprintf(fid, '\\caption{Average One-Step Ahead Joint Log Predictive Scores}\n');
fprintf(fid, '\\begin{tabular}{|l|%s|}\n\\hline\n', repmat('c', 1, Nmodels));
fprintf(fid, '\\textbf{Period}');
for mm = 1:Nmodels
    fprintf(fid, ' & \\textbf{%s}', modelLabels{mm});
end
fprintf(fid, ' \\\\\n\\hline\n');
for rr = 1:Nregimes
    ridx = regime_idx{rr};
    fprintf(fid, '%s', regime_names{rr});
    for mm = 1:Nmodels
        vals = JLS(ridx, mm);
        fprintf(fid, ' & %.2f', mean(vals, 'omitnan'));
    end
    fprintf(fid, ' \\\\\n');
end
fprintf(fid, '\\hline\n\\end{tabular}\n\\end{table}\n');
fclose(fid);
fprintf('  Wrote %s\n', fname_jls);

fprintf('\n=== Post-processing complete ===\n');
fprintf('Output in: %s\n', output_dir);
fprintf('Console log saved to: %s\n', [output_dir 'console_output.txt']);
diary off;

%% ===================== Helper Functions ==================================

function [tstat, pval] = dm_west_test(d, h)
% DM-West test with Newey-West HAC standard errors
% d: vector of loss differentials (positive = model beats benchmark)
% h: forecast horizon (determines bandwidth)
    T = length(d);
    dbar = mean(d);

    % Newey-West bandwidth: max(h-1, 1) for multi-step forecasts
    bw = max(h-1, 1);

    % HAC variance
    gamma0 = var(d);
    gammasum = 0;
    for lag = 1:bw
        gam_lag = mean((d(1+lag:end) - dbar) .* (d(1:end-lag) - dbar));
        w = 1 - lag/(bw+1);  % Bartlett kernel
        gammasum = gammasum + 2 * w * gam_lag;
    end
    V_nw = gamma0 + gammasum;

    if V_nw <= 0
        tstat = 0; pval = 1; return;
    end

    tstat = dbar / sqrt(V_nw / T);
    pval = 2 * (1 - normcdf(abs(tstat)));
end

function stars = significance_stars(pval, ratio)
% Console stars: show when ratio < 1 (DPM-OISV in numerator has lower loss)
    if ratio >= 1
        stars = '';
    elseif pval < 0.001
        stars = '***';
    elseif pval < 0.01
        stars = '**';
    elseif pval < 0.05
        stars = '*';
    else
        stars = '';
    end
end

function stars = significance_stars_tex(pval, ratio)
% LaTeX stars: show when ratio < 1 (DPM-OISV in numerator has lower loss)
    if ratio >= 1
        stars = '';
    elseif pval < 0.001
        stars = '\textsuperscript{***}';
    elseif pval < 0.01
        stars = '\textsuperscript{**}';
    elseif pval < 0.05
        stars = '\textsuperscript{*}';
    else
        stars = '';
    end
end

function out = iif(cond, a, b)
    if cond, out = a; else, out = b; end
end
