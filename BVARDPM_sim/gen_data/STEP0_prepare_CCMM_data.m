% =========================================================================
% STEP0_prepare_CCMM_data.m
%
% Loads the January 2026 FRED-MD vintage (2026-01-MD.csv), extracts the
% 16 CCMM variables, applies the exact transformations from Table S.1 of
% Carriero, Clark, Marcellino & Mertens (2022, REStat), and saves the
% result ready for SVo estimation.
%
% Raw file structure of 202601MD.csv:
%   Row 1  : variable names  (sasdate, RPI, W875RX1, ...)
%   Row 2  : tcode numbers   (1,2,2,...) -- we IGNORE these, apply our own
%   Rows 3+: monthly data, starting 1959:M01
%
% Transformations applied (CCMM Table S.1):
%   dlog  -> 1200 * (log(x_t) - log(x_{t-1}))   [annualised % growth rate]
%   log   -> log(x_t)                             [log level]
%   level -> x_t                                  [no transformation]
%
% After transformation the usable sample starts 1959:M03 (two dlog vars
% lose one observation each; the binding constraint is losing the first obs
% for dlog variables, so transformed data starts at row 2 of levels = Feb
% 1959 for dlog, but CCMM convention uses 1959:M03 as the first obs
% because the data originally starts Jan 1959 and tcode row must be dropped).
%
% OUTPUTS (saved to CCMM16_data.mat):
%   data_trans     : (T_trans x 16) full transformed dataset
%   dates_trans    : (T_trans x 1)  datenum vector
%   Y0             : (p x 16)       pre-sample initial conditions
%   shortY         : (T x 16)       estimation sample
%   dates_shortY   : (T x 1)        dates for shortY
%   n, T, p        : scalars
%   labels, codes  : variable info
%   minnesotaPriorMean : (16 x 1)   Minnesota prior means on own-lag 1
% =========================================================================

clear; clc;

%% ---- 0. Parameters ----
p         = 12;           % VAR lag order (CCMM uses p=12 for monthly)
data_file = '2026-01-MD.csv';

%% ---- 1. The 16 CCMM variables (Table S.1 / Table 1 in CCMM 2022) ----
%
%  col 1: FRED-MD column name  (exactly as in CSV header)
%  col 2: descriptive label
%  col 3: transformation  ('dlog', 'log', or 'level')
%  col 4: Minnesota prior mean on first own-lag coefficient

var_defs = {
    'RPI',               'Real Income',           'dlog',  0;
    'DPCERA3M086SBEA',   'Real Consumption',      'dlog',  0;
    'INDPRO',            'IP',                    'dlog',  0;
    'CUMFNS',            'Capacity Utilization',  'level', 1;
    'UNRATE',            'Unemployment Rate',     'level', 1;
    'PAYEMS',            'Nonfarm Payrolls',      'dlog',  0;
    'CES0600000007',     'Hours',                 'level', 0;
    'CES0600000008',     'Hourly Earnings',       'dlog',  0;
    'WPSFD49207',        'PPI (Fin. Goods)',      'dlog',  1;
    'PCEPI',             'PCE Prices',            'dlog',  1;
    'HOUST',             'Housing Starts',        'log',   1;
    'S&P 500',           'S&P 500',               'dlog',  0;
    'EXUSUKx',           'USD/GBP FX Rate',       'dlog',  0;
    'GS5',               '5-Year Yield',          'level', 1;
    'GS10',              '10-Year Yield',         'level', 1;
    'BAAFFM',            'Baa Spread',            'level', 1;
};

n_vars             = size(var_defs, 1);    % 16
codes              = var_defs(:,1);        % FRED-MD column names
labels             = var_defs(:,2);        % descriptive labels
tcodes_str         = var_defs(:,3);        % transformation strings
mn_prior_raw       = cell2mat(var_defs(:,4)); % 16x1 Minnesota prior means

%% ---- 2. Load raw CSV ----
fprintf('Loading %s ...\n', data_file);

% Step 2a: Detect options and force 'sasdate' to be read as TEXT.
% Without this, MATLAB auto-parses it as datetime, which fails
% on the tcode row (which contains numbers, not dates).
opts = detectImportOptions(data_file);
opts.VariableNamingRule = 'preserve';   % keep exact names e.g. 'S&P 500'
opts = setvartype(opts, 'sasdate', 'char');   % read date column as text
raw_tbl = readtable(data_file, opts);

% Step 2b: Drop the tcode row.
% 202601MD.csv structure:
%   Header row  : variable names
%   Data row 1  : tcode numbers (e.g. "5", "1", "2" ...) -- NOT dates
%   Data rows 2+: monthly observations starting 1959:M01
% A real date contains '/' (e.g. '1/1/1959'); tcode row has a plain number.
first_sasdate_val = raw_tbl{1, 'sasdate'};
if iscell(first_sasdate_val), first_sasdate_val = first_sasdate_val{1}; end
first_sasdate_str = strtrim(char(first_sasdate_val));

if isempty(strfind(first_sasdate_str, '/'))
    raw_tbl(1, :) = [];
    fprintf('  Dropped tcode row (sasdate was "%s").\n', first_sasdate_str);
else
    fprintf('  No tcode row (first sasdate is "%s").\n', first_sasdate_str);
end

% Step 2c: Parse sasdate strings -> datenums.
% FRED-MD format: M/D/YYYY (e.g. "1/1/1959"), no leading zeros.
date_col = raw_tbl{:, 'sasdate'};

if iscell(date_col)
    date_strings = date_col;
elseif isstring(date_col)
    date_strings = cellstr(date_col);
elseif isdatetime(date_col)
    dn_raw = datenum(date_col);
    date_strings = {};
elseif isnumeric(date_col)
    dn_raw = date_col;
    date_strings = {};
else
    date_strings = cellstr(string(date_col));
end

if exist('date_strings','var') && ~isempty(date_strings)
    try
        dn_raw = datenum(date_strings, 'M/d/yyyy');
    catch
        try
            dn_raw = datenum(date_strings, 'mm/dd/yyyy');
        catch
            dn_raw = datenum(date_strings);  % let MATLAB guess as last resort
        end
    end
end

Nraw = length(dn_raw);
fprintf('  Raw data: %d obs,  %s  to  %s\n', Nraw, ...
    datestr(dn_raw(1),'yyyy:mmm'), datestr(dn_raw(end),'yyyy:mmm'));

%% ---- 3. Extract raw level columns for the 16 variables ----
raw_levels = NaN(Nraw, n_vars);
all_col_names = raw_tbl.Properties.VariableNames;

for j = 1:n_vars
    cname = codes{j};
    
    % Direct match first
    col_idx = find(strcmp(all_col_names, cname));
    
    % Fallback for 'S&P 500': FRED-MD sometimes stores as 'S.P.500' or 'SP500'
    if isempty(col_idx)
        col_idx = find(strcmpi(all_col_names, strrep(cname,' ','')) | ...
                       strcmpi(all_col_names, strrep(cname,'&','.')));
    end
    if isempty(col_idx) && strcmp(cname,'S&P 500')
        col_idx = find(contains(all_col_names,'S') & ...
                       (contains(all_col_names,'P') | contains(all_col_names,'p')) & ...
                       contains(all_col_names,'500'));
    end
    
    if isempty(col_idx)
        error('Cannot locate column for variable %d: "%s"', j, cname);
    end
    raw_levels(:, j) = raw_tbl{:, col_idx(1)};
    fprintf('  Var %2d: %-22s <- column "%s"\n', j, labels{j}, all_col_names{col_idx(1)});
end

%% ---- 4. Apply transformations row by row ----
%
% dlog : transformed(t) = 1200*(log(raw(t)) - log(raw(t-1))), so row 1 is NaN
% log  : transformed(t) = log(raw(t)),  no NaN from transformation
% level: transformed(t) = raw(t),        no NaN from transformation

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

%% ---- 5. Find usable sample (all 16 non-NaN) ----
complete_rows = all(~isnan(transformed), 2);

% First complete row (row 2 for dlog variables = 1959:M02,
% but some series may have early NaNs in the raw data too)
first_row = find(complete_rows, 1, 'first');

% Last complete row: some Dec-2025 entries may be NaN
last_row  = find(complete_rows, 1, 'last');

fprintf('\n  First complete row: %d  (%s)\n', first_row, datestr(dn_raw(first_row),'yyyy:mmm'));
fprintf('  Last  complete row: %d  (%s)\n', last_row,  datestr(dn_raw(last_row), 'yyyy:mmm'));

data_trans  = transformed(first_row:last_row, :);   % T_trans x 16
dates_trans = dn_raw(first_row:last_row);            % T_trans x 1
T_trans     = size(data_trans, 1);
fprintf('  Total usable observations: %d\n', T_trans);

%% ---- 6. Split into pre-sample Y0 and estimation sample shortY ----
%
%  Rows 1..p        : Y0     (pre-sample, used to form lags at t=1)
%  Rows p+1..T_trans: shortY (estimation sample)

if T_trans <= p
    error('T_trans=%d <= p=%d. Not enough data.', T_trans, p);
end

Y0           = data_trans(1:p,     :);      % p   x 16
shortY       = data_trans(p+1:end, :);      % T   x 16
dates_Y0     = dates_trans(1:p);
dates_shortY = dates_trans(p+1:end);
T            = size(shortY, 1);

fprintf('\n  Pre-sample  Y0    : %d obs  (%s to %s)\n', p, ...
    datestr(dates_Y0(1),'yyyy:mmm'), datestr(dates_Y0(end),'yyyy:mmm'));
fprintf('  Estimation shortY : %d obs  (%s to %s)\n', T, ...
    datestr(dates_shortY(1),'yyyy:mmm'), datestr(dates_shortY(end),'yyyy:mmm'));
fprintf('  n=%d variables,  p=%d lags,  T=%d estimation obs\n', n_vars, p, T);

%% ---- 7. Minnesota prior means (CCMM Table S.1) ----
%
% Variables with prior mean 1 on first own-lag (level variables):
%   CUMFNS, UNRATE, WPSFD49207, PCEPI, HOUST, GS5, GS10, BAAFFM
% All others: prior mean 0

minnesotaPriorMean = mn_prior_raw;   % 16 x 1  (already set from var_defs)
n = n_vars;

%% ---- 8. Save ----
save('CCMM16_data.mat', ...
    'data_trans', 'dates_trans', 'T_trans', ...
    'Y0', 'shortY', 'dates_Y0', 'dates_shortY', ...
    'n', 'T', 'p', ...
    'labels', 'codes', 'tcodes_str', ...
    'mn_prior_raw', 'minnesotaPriorMean', ...
    'first_row', 'last_row');
fprintf('\nSaved CCMM16_data.mat\n');
fprintf('Next step: run STEP1_run_SVo_fullsample.m\n\n');

%% ---- 9. Quick diagnostic plot ----
figure('Name','CCMM 16 Variables','NumberTitle','off');
for j = 1:n_vars
    subplot(4,4,j);
    plot(dates_shortY, shortY(:,j), 'b-', 'LineWidth', 0.7);
    datetick('x','yy','keeplimits');
    title(labels{j}, 'FontSize', 6, 'Interpreter','none');
    box off;
end
sgtitle(sprintf('CCMM 16 Variables — Estimation Sample (%s to %s)', ...
    datestr(dates_shortY(1),'yyyy:mmm'), datestr(dates_shortY(end),'yyyy:mmm')));
