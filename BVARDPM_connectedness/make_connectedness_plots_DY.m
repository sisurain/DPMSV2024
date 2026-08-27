% =========================================================================
% make_connectedness_plots_DY.m
%
% Headline figure for the DY application: total system-wide connectedness
% C^H_t time series for DPM-OISV and SV, with 68% credible bands and
% shaded crisis windows.
%
% Two views generated:
%   (1) full sample (2003-09 to 2014-02)
%   (2) crisis window zoom (2007 to 2014) for finer resolution
%
% Reads from:
%   results/results_DPM_OIMSV_DY.mat
%   results/results_SV_DY.mat
%
% Writes to:
%   dy_figures/connectedness_total_full.png/.fig
%   dy_figures/connectedness_total_zoom.png/.fig
%   dy_figures/connectedness_total_data.csv  (medians + bands for re-plotting)
% =========================================================================

clear; clc; close all;

%% ===================== Configuration ====================================
results_dir = 'results/';
output_dir  = 'dy_figures/';
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

modelNames  = {'DPM_OIMSV', 'SV'};
modelLabels = {'DPM-OISV', 'SV'};
modelFiles  = {'results_DPM_OIMSV_DY.mat', 'results_SV_DY.mat'};
Nmodels     = length(modelNames);

% Color scheme
colorDPM = [0.1216, 0.4667, 0.7059];   % blue
colorSV  = [0.8902, 0.1020, 0.1098];   % red
modelColors = [colorDPM; colorSV];

% Credible band quantiles
band_lo = 16;
band_hi = 84;

% Crisis windows for shading (drawn behind the lines)
crisis_windows = {
    {'2007-08-01', '2009-06-30', 'Subprime / GFC'};
    {'2010-04-01', '2012-07-31', 'Greek / Euro debt crisis'};
    {'2013-05-01', '2013-12-31', 'Taper Tantrum'};
};
crisis_color = [0.85 0.85 0.85];    % light grey
crisis_alpha = 0.45;

%% ===================== Load Results =====================================
fprintf('Loading results...\n');
R = cell(Nmodels, 1);
for mm = 1:Nmodels
    fname = [results_dir modelFiles{mm}];
    if ~exist(fname, 'file')
        error('Missing results file: %s', fname);
    end
    R{mm} = load(fname, 'store_C_total', 'dates_subsample');
    fprintf('  Loaded %s (%d draws x %d time points)\n', modelLabels{mm}, ...
            size(R{mm}.store_C_total, 2), size(R{mm}.store_C_total, 1));
end

% Verify the two models share the same time grid
if length(R{1}.dates_subsample) ~= length(R{2}.dates_subsample) || ...
   any(R{1}.dates_subsample ~= R{2}.dates_subsample)
    error('The two model results files have different time grids; cannot overlay.');
end
dates_sub = R{1}.dates_subsample(:);   % column vector
T_sub     = length(dates_sub);

%% ===================== Compute summaries ================================
% Posterior median + 16/84 percentile bands across draws (axis 2)
C_med = zeros(T_sub, Nmodels);
C_lo  = zeros(T_sub, Nmodels);
C_hi  = zeros(T_sub, Nmodels);
for mm = 1:Nmodels
    Cm = R{mm}.store_C_total;       % T_sub x nsims
    C_med(:, mm) = median(Cm, 2);
    C_lo(:, mm)  = prctile(Cm, band_lo, 2);
    C_hi(:, mm)  = prctile(Cm, band_hi, 2);
end

%% ===================== Two views ========================================
views = {
    struct('name','full',  'xlim', [],                                                      'suffix','full');
    struct('name','zoom',  'xlim', [datenum(2007,1,1), datenum(2014,3,1)],                  'suffix','zoom');
};

for vw = 1:length(views)
    view_cfg = views{vw};
    fprintf('\nGenerating %s-sample plot...\n', view_cfg.name);

    fig = figure('Name', sprintf('Connectedness: %s', view_cfg.name), ...
                 'Position', [100 100 1200 540], ...
                 'Color', 'white');

    hold on;

    % Determine y-axis range first so the crisis patches cover the full panel
    if isempty(view_cfg.xlim)
        in_view = true(T_sub, 1);
    else
        in_view = dates_sub >= view_cfg.xlim(1) & dates_sub <= view_cfg.xlim(2);
    end
    y_data_min = min(min(C_lo(in_view, :)));
    y_data_max = max(max(C_hi(in_view, :)));
    y_pad = 0.05 * (y_data_max - y_data_min);
    y_lim = [y_data_min - y_pad, y_data_max + y_pad];

    % --- Crisis window shading (drawn first, behind everything) ---
    for cw = 1:length(crisis_windows)
        c_start = datenum(crisis_windows{cw}{1}, 'yyyy-mm-dd');
        c_end   = datenum(crisis_windows{cw}{2}, 'yyyy-mm-dd');
        % Clip to view window if zoomed
        if ~isempty(view_cfg.xlim)
            c_start = max(c_start, view_cfg.xlim(1));
            c_end   = min(c_end,   view_cfg.xlim(2));
        end
        if c_end > c_start
            patch([c_start c_end c_end c_start], ...
                  [y_lim(1) y_lim(1) y_lim(2) y_lim(2)], ...
                  crisis_color, ...
                  'FaceAlpha', crisis_alpha, ...
                  'EdgeColor', 'none', ...
                  'HandleVisibility', 'off');
        end
    end

    % --- 68% credible bands (translucent fill) ---
    for mm = 1:Nmodels
        % Build closed polygon: forward along upper, back along lower
        x_poly = [dates_sub; flipud(dates_sub)];
        y_poly = [C_hi(:, mm); flipud(C_lo(:, mm))];
        patch(x_poly, y_poly, modelColors(mm, :), ...
              'FaceAlpha', 0.18, ...
              'EdgeColor', 'none', ...
              'HandleVisibility', 'off');
    end

    % --- Median lines (drawn on top) ---
    h_lines = gobjects(Nmodels, 1);
    for mm = 1:Nmodels
        h_lines(mm) = plot(dates_sub, C_med(:, mm), '-', ...
                           'Color', modelColors(mm, :), ...
                           'LineWidth', 1.6);
    end

    % --- Crisis labels at top of plot (only on full view to avoid clutter) ---
    if isempty(view_cfg.xlim)
        for cw = 1:length(crisis_windows)
            c_start = datenum(crisis_windows{cw}{1}, 'yyyy-mm-dd');
            c_end   = datenum(crisis_windows{cw}{2}, 'yyyy-mm-dd');
            c_mid   = (c_start + c_end) / 2;
            text(c_mid, y_lim(2) - 0.04*(y_lim(2)-y_lim(1)), ...
                 crisis_windows{cw}{3}, ...
                 'HorizontalAlignment', 'center', ...
                 'VerticalAlignment',   'top', ...
                 'FontSize', 8, ...
                 'Color', [0.35 0.35 0.35], ...
                 'FontAngle', 'italic');
        end
    end

    hold off;
    box on;
    grid on;
    set(gca, 'GridAlpha', 0.15, 'Layer', 'top', 'FontSize', 11);

    % Axes
    ylim(y_lim);
    if ~isempty(view_cfg.xlim)
        xlim(view_cfg.xlim);
    else
        xlim([dates_sub(1)-30, dates_sub(end)+30]);
    end

    if isempty(view_cfg.xlim)
        datetick('x', 'yyyy', 'keeplimits');
    else
        datetick('x', 'mmm yy', 'keeplimits');
    end

    xlabel('Date');
    ylabel('Total connectedness C^H_t  (% of system FEV)');
    title(sprintf('Total Connectedness: DPM-OISV vs SV  (%s sample, n=30 banks)', ...
                  view_cfg.name), ...
          'FontWeight', 'normal', 'FontSize', 12);

    % Legend with credible-band note in the entries
    legend_labels = cell(Nmodels, 1);
    for mm = 1:Nmodels
        legend_labels{mm} = sprintf('%s (median, %d-%d%% band)', ...
                                    modelLabels{mm}, band_lo, band_hi);
    end
    lh = legend(h_lines, legend_labels, ...
                'Location', 'southwest', ...
                'Box', 'off', ...
                'FontSize', 10);

    % Save
    fname_base = sprintf('connectedness_total_%s', view_cfg.suffix);
    saveas(fig, [output_dir fname_base '.png']);
    savefig(fig, [output_dir fname_base '.fig']);
    fprintf('  Saved %s.png and .fig\n', fname_base);
end

%% ===================== Save data as CSV for re-plotting =================
csv_path = [output_dir 'connectedness_total_data.csv'];
fid = fopen(csv_path, 'w');
fprintf(fid, 'date,date_num,DPM_median,DPM_q16,DPM_q84,SV_median,SV_q16,SV_q84\n');
for t = 1:T_sub
    fprintf(fid, '%s,%.2f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n', ...
        datestr(dates_sub(t), 'yyyy-mm-dd'), dates_sub(t), ...
        C_med(t,1), C_lo(t,1), C_hi(t,1), ...
        C_med(t,2), C_lo(t,2), C_hi(t,2));
end
fclose(fid);
fprintf('\nSaved %s\n', csv_path);

%% ===================== Print summary statistics =========================
fprintf('\n----- Headline summary -----\n');
fprintf('  %-10s %12s %12s %12s\n', 'Model', 'mean(med)', 'min(med)', 'max(med)');
for mm = 1:Nmodels
    fprintf('  %-10s %12.2f %12.2f %12.2f\n', modelLabels{mm}, ...
            mean(C_med(:,mm)), min(C_med(:,mm)), max(C_med(:,mm)));
end
fprintf('\n  Dynamic range (max - min of medians):\n');
for mm = 1:Nmodels
    fprintf('    %-10s: %.2f pp\n', modelLabels{mm}, ...
            max(C_med(:,mm)) - min(C_med(:,mm)));
end
fprintf('\nDone.\n');
