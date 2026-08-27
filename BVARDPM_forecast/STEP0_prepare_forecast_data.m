% =========================================================================
% STEP0_prepare_forecast_data.m
%
% Loads the January 2026 FRED-MD vintage, extracts 42 variables spanning
% all major macroeconomic categories, applies CCMM-style transformations,
% and saves the result for use in the forecasting exercise.
%
% CCMM transformation convention (less differencing than standard FRED-MD):
%   dlog  = 1200*(log(x_t) - log(x_{t-1}))   [annualised % growth rate]
%   log   = log(x_t)                           [log level]
%   level = x_t                                [no transformation]
%
% EXCLUDED:
%   ACOGNO  — raw data starts Feb 1992 only (would lose 33 yrs of sample)
%   UMCSENTx — quarterly before ~1978, creates interior NaN every other month
%
% OUTPUT: data/FREDMD42_forecast_data.mat
% =========================================================================

clear; clc;

%% ---- 0. Parameters ----
p         = 12;
data_file = 'data/2026-01-MD.csv';

%% ---- 1. Variable definitions: 42 variables ----
var_defs = {
    % === Group 1: Output & Income (5) ===
    'RPI',               'Real Income (*)',          'dlog',  0;
    'W875RX1',           'Real Income ex Transf',    'dlog',  0;
    'DPCERA3M086SBEA',   'Real Consumption (*)',     'dlog',  0;
    'INDPRO',            'IP (*)',                   'dlog',  0;
    'CUMFNS',            'Capacity Util (*)',        'level', 1;
    % === Group 2: Labor Market (7) ===
    'CE16OV',            'Civilian Employment',      'dlog',  0;
    'UNRATE',            'Unemployment Rate (*)',    'level', 1;
    'UEMPMEAN',          'Mean Unemp Duration',     'level', 1;
    'CLAIMSx',           'Initial Claims',          'dlog',  0;
    'PAYEMS',            'Nonfarm Payrolls (*)',     'dlog',  0;
    'CES0600000007',     'Hours (*)',               'level', 0;
    'AWHMAN',            'Hours: Manufacturing',    'level', 0;
    % === Group 3: Housing (2) ===
    'HOUST',             'Housing Starts (*)',       'log',   1;
    'PERMIT',            'Building Permits',         'log',   1;
    % === Group 4: Orders & Inventories (3) ===
    'CMRMTSPLx',         'Real Mfg & Trade Sales',  'dlog',  0;
    'AMDMNOx',           'Durable Goods Orders',    'dlog',  0;
    'BUSINVx',           'Business Inventories',    'dlog',  0;
    % === Group 5: Money & Credit (4) ===
    'M1SL',              'M1 Money Supply',          'dlog',  1;
    'M2SL',              'M2 Money Supply',          'dlog',  1;
    'BUSLOANS',          'Business Loans',           'dlog',  1;
    'NONREVSL',          'Consumer Credit (NR)',     'dlog',  1;
    % === Group 6: Interest Rates (6) ===
    'FEDFUNDS',          'Fed Funds Rate',           'level', 1;
    'TB3MS',             '3-Mo T-Bill',              'level', 1;
    'GS1',               '1-Year Treasury',          'level', 1;
    'GS5',               '5-Year Treasury (*)',      'level', 1;
    'GS10',              '10-Year Treasury (*)',     'level', 1;
    'BAA',               'Baa Corporate',            'level', 1;
    % === Group 7: Spreads (4) ===
    'BAAFFM',            'Baa Spread (*)',           'level', 1;
    'AAAFFM',            'AAA Spread',               'level', 1;
    'T5YFFM',            '5Y-FF Spread',             'level', 1;
    'T10YFFM',           '10Y-FF Spread',            'level', 1;
    % === Group 8: Prices (6) ===
    'WPSFD49207',        'PPI: Fin Goods (*)',       'dlog',  1;
    'PCEPI',             'PCE Prices (*)',           'dlog',  1;
    'CPIAUCSL',          'CPI: All Items',           'dlog',  1;
    'OILPRICEx',         'Oil Price',                'dlog',  0;
    'PPICMM',            'PPI: Metals & Metal Prod', 'dlog', 0;
    'CES0600000008',     'Hourly Earnings (*)',      'dlog',  0;
    % === Group 9: Stock Market (3) ===
    'S&P 500',           'S&P 500 (*)',              'dlog',  0;
    'S&P div yield',     'S&P Div Yield',            'level', 1;
    'VIXCLSx',           'VIX',                      'level', 0;
    % === Group 10: Exchange Rates (2) ===
    'EXUSUKx',           'USD/GBP FX (*)',           'dlog',  0;
    'EXJPUSx',           'JPY/USD FX',               'dlog',  0;
};

n_vars             = size(var_defs, 1);
codes              = var_defs(:,1);
labels             = var_defs(:,2);
tcodes_str         = var_defs(:,3);
mn_prior_raw       = cell2mat(var_defs(:,4));

fprintf('Selected %d variables for forecasting exercise\n', n_vars);

%% ---- 2. Identify CCMM 16 ----
ccmm16_codes = {'RPI','DPCERA3M086SBEA','INDPRO','CUMFNS','UNRATE','PAYEMS', ...
    'CES0600000007','CES0600000008','WPSFD49207','PCEPI','HOUST', ...
    'S&P 500','EXUSUKx','GS5','GS10','BAAFFM'};
is_ccmm16 = false(n_vars, 1);
ccmm16_ndx = NaN(16, 1);
for j = 1:16
    idx = find(strcmp(codes, ccmm16_codes{j}));
    if ~isempty(idx)
        is_ccmm16(idx) = true;
        ccmm16_ndx(j)  = idx;
    end
end
fprintf('  Of these, %d are from the original CCMM 16-variable set\n', sum(is_ccmm16));

%% ---- 3. Load raw CSV ----
fprintf('Loading %s ...\n', data_file);
opts = detectImportOptions(data_file);
opts.VariableNamingRule = 'preserve';
opts = setvartype(opts, 'sasdate', 'char');
raw_tbl = readtable(data_file, opts);

first_sasdate_val = raw_tbl{1, 'sasdate'};
if iscell(first_sasdate_val), first_sasdate_val = first_sasdate_val{1}; end
first_sasdate_str = strtrim(char(first_sasdate_val));
if isempty(strfind(first_sasdate_str, '/'))
    raw_tbl(1, :) = [];
    fprintf('  Dropped tcode row.\n');
end

date_col = raw_tbl{:, 'sasdate'};
if iscell(date_col), date_strings = date_col;
elseif isstring(date_col), date_strings = cellstr(date_col);
elseif isdatetime(date_col), dn_raw = datenum(date_col); date_strings = {};
else, date_strings = cellstr(string(date_col));
end
if exist('date_strings','var') && ~isempty(date_strings)
    try, dn_raw = datenum(date_strings, 'M/d/yyyy');
    catch
        try, dn_raw = datenum(date_strings, 'mm/dd/yyyy');
        catch, dn_raw = datenum(date_strings);
        end
    end
end

Nraw = length(dn_raw);
fprintf('  Raw data: %d obs,  %s  to  %s\n', Nraw, ...
    datestr(dn_raw(1),'yyyy:mmm'), datestr(dn_raw(end),'yyyy:mmm'));

%% ---- 4. Extract raw level columns ----
raw_levels = NaN(Nraw, n_vars);
all_col_names = raw_tbl.Properties.VariableNames;

for j = 1:n_vars
    cname = codes{j};
    col_idx = find(strcmp(all_col_names, cname));
    if isempty(col_idx)
        col_idx = find(strcmpi(all_col_names, strrep(cname,' ','')) | ...
                       strcmpi(all_col_names, strrep(cname,'&','.')));
    end
    if isempty(col_idx) && contains(cname, 'S&P')
        all_sp = find(contains(all_col_names, 'S') & contains(all_col_names, 'P'));
        for cc = 1:length(all_sp)
            if contains(cname, 'div') == contains(all_col_names{all_sp(cc)}, 'div')
                col_idx = all_sp(cc); break;
            end
        end
    end
    if isempty(col_idx)
        error('Cannot locate column for variable %d: "%s"', j, cname);
    end
    raw_levels(:, j) = raw_tbl{:, col_idx(1)};
    fprintf('  Var %2d: %-28s <- "%s"\n', j, labels{j}, all_col_names{col_idx(1)});
end

%% ---- 5. Apply CCMM-style transformations ----
transformed = NaN(Nraw, n_vars);
for j = 1:n_vars
    x = raw_levels(:, j);
    switch tcodes_str{j}
        case 'dlog'
            lx = log(x);
            transformed(2:end, j) = 1200 * diff(lx);
        case 'log'
            transformed(:, j) = log(x);
        case 'level'
            transformed(:, j) = x;
    end
end

%% ---- 6. Find usable sample ----
complete_rows = all(~isnan(transformed), 2);
first_row = find(complete_rows, 1, 'first');
last_row  = find(complete_rows, 1, 'last');

fprintf('\n  First complete row: %d  (%s)\n', first_row, datestr(dn_raw(first_row),'yyyy:mmm'));
fprintf('  Last  complete row: %d  (%s)\n', last_row,  datestr(dn_raw(last_row), 'yyyy:mmm'));

data_trans  = transformed(first_row:last_row, :);
dates_trans = dn_raw(first_row:last_row);
T_trans     = size(data_trans, 1);

% Final NaN check
assert(~any(isnan(data_trans(:))), 'ERROR: data_trans still contains NaN!');
fprintf('  Total usable observations: %d  (verified NaN-free)\n', T_trans);

%% ---- 7. Split ----
Y0           = data_trans(1:p, :);
shortY       = data_trans(p+1:end, :);
dates_Y0     = dates_trans(1:p);
dates_shortY = dates_trans(p+1:end);
T            = size(shortY, 1);

fprintf('\n  Pre-sample  Y0    : %d obs  (%s to %s)\n', p, ...
    datestr(dates_Y0(1),'yyyy:mmm'), datestr(dates_Y0(end),'yyyy:mmm'));
fprintf('  Estimation shortY : %d obs  (%s to %s)\n', T, ...
    datestr(dates_shortY(1),'yyyy:mmm'), datestr(dates_shortY(end),'yyyy:mmm'));
fprintf('  n=%d variables,  p=%d lags,  T=%d obs\n', n_vars, p, T);

%% ---- 8. Minnesota prior means ----
minnesotaPriorMean = mn_prior_raw;
n = n_vars;

%% ---- 9. Cumulative codes ----
cumcode = strcmp(tcodes_str, 'dlog');
fprintf('\n  %d variables use cumulative forecast evaluation (dlog)\n', sum(cumcode));

%% ---- 10. Forecasting parameters ----
date_preCOVID_end    = datenum(2020, 2, 1);
date_COVID_start     = datenum(2020, 3, 1);
date_COVID_end       = datenum(2021, 12, 1);
date_postCOVID_start = datenum(2022, 1, 1);

fcstNhorizons = 12;
fcstHorizons_report = [1, 6, 12];

fcst_firstOrigin_date = datenum(2004, 12, 1);
fcst_firstOrigin = find(dates_trans == fcst_firstOrigin_date);
if isempty(fcst_firstOrigin)
    [~, fcst_firstOrigin] = min(abs(dates_trans - fcst_firstOrigin_date));
end
fprintf('  First forecast origin: row %d (%s)\n', ...
    fcst_firstOrigin, datestr(dates_trans(fcst_firstOrigin),'yyyy:mmm'));
fprintf('  Effective T at first origin: %d\n', fcst_firstOrigin - p);

Njumpoffs = T_trans - fcst_firstOrigin + 1;
fprintf('  Total forecast origins: %d  (%s to %s)\n', Njumpoffs, ...
    datestr(dates_trans(fcst_firstOrigin),'yyyy:mmm'), ...
    datestr(dates_trans(end),'yyyy:mmm'));

%% ---- 11. Save ----
save('data/FREDMD42_forecast_data.mat', ...
    'data_trans', 'dates_trans', 'T_trans', ...
    'Y0', 'shortY', 'dates_Y0', 'dates_shortY', ...
    'n', 'T', 'p', ...
    'labels', 'codes', 'tcodes_str', ...
    'mn_prior_raw', 'minnesotaPriorMean', ...
    'first_row', 'last_row', ...
    'is_ccmm16', 'ccmm16_ndx', 'ccmm16_codes', ...
    'cumcode', ...
    'date_preCOVID_end', 'date_COVID_start', 'date_COVID_end', ...
    'date_postCOVID_start', ...
    'fcstNhorizons', 'fcstHorizons_report', 'fcst_firstOrigin', 'Njumpoffs');
fprintf('\nSaved data/FREDMD42_forecast_data.mat\n');

%% ---- 12. Summary ----
fprintf('\n========================================\n');
fprintf('VARIABLE SUMMARY TABLE\n');
fprintf('========================================\n');
fprintf('%3s  %-28s  %-6s  %2s  %s\n', '#', 'Label', 'Trans', 'mn', 'CCMM16');
fprintf('%s\n', repmat('-',1,60));
for j = 1:n_vars
    ccmm_flag = '';
    if is_ccmm16(j), ccmm_flag = '  *'; end
    fprintf('%3d  %-28s  %-6s  %2d%s\n', j, labels{j}, tcodes_str{j}, mn_prior_raw(j), ccmm_flag);
end
fprintf('Total: %d variables  (* = original CCMM 16)\n', n_vars);

%% ---- 13. Plot ----
figure('Name','42 Forecasting Variables','NumberTitle','off','Position',[50 50 1400 900]);
ncols = 7; nrows = ceil(n_vars/ncols);
for j = 1:n_vars
    subplot(nrows, ncols, j);
    plot(dates_shortY, shortY(:,j), 'b-', 'LineWidth', 0.5);
    datetick('x','yy','keeplimits');
    title(labels{j}, 'FontSize', 5, 'Interpreter','none');
    box off; set(gca,'FontSize',4);
end
sgtitle(sprintf('42 Variables: %s to %s', ...
    datestr(dates_shortY(1),'yyyy:mmm'), datestr(dates_shortY(end),'yyyy:mmm')));
fprintf('\nDone.\n');
