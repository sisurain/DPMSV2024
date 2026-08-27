% =========================================================================
% load_ddly_data.m
%
% Reads the Demirer-Diebold-Liu-Yilmaz (2018, JAE) bank volatility panel
% from ddly-data.csv, takes the top n_keep banks by assets, transforms to
% log-variance, and saves a .mat file with the same interface as the
% existing macro pipeline (Y0, shortY, dates_shortY, labels, n, T, p).
%
% =========================================================================
% NOTE ON THE LOG-VARIANCE TRANSFORM
% =========================================================================
% The raw CSV contains daily Garman-Klass (1980) range-based variance
% estimates, which are O(1e-6) to O(1e-1). We transform to LOG VARIANCES:
%
%   y_{i,t} = log( GK_variance_{i,t} )
%
% Rationale:
%   1. Log realized variances are approximately Gaussian (Andersen-Bollerslev-
%      Diebold-Labys 2003, JFE), which makes the conditionally-Gaussian VAR
%      assumption substantially more credible than working on the variance
%      scale where the data span 5+ orders of magnitude.
%   2. This is the standard transformation in the realized-vol VAR literature
%      and in Diebold-Yilmaz (2014) connectedness applications.
%   3. The DPM-OISV outlier mechanism still has plenty to absorb: even log
%      variances exhibit sharp jumps at Lehman, the European debt crisis,
%      and COVID. The DPM clusters detect these as discrete shifts rather
%      than the SV process having to track 4-orders-of-magnitude variation.
%   4. The Cholesky-SV implementation in run_SV_fullsample.m would be
%      numerically unstable on raw variances at this scale.
%
% Note that Chan-Yu (2022, ICEEE) Section 6.2 describe their input as the
% raw GK variance series, but their reported system-wide connectedness range
% (70-82) is consistent with a log-vol VAR; the literal raw-variance reading
% is unlikely to match their numerics.
% =========================================================================
%
% Inputs (set in the configuration block below):
%   csv_path  : path to ddly-data.csv
%   n_keep    : number of top banks to retain (default 30)
%   p         : number of VAR lags (default 3, Diebold-Yilmaz standard)
%   out_path  : output .mat filename
%
% Output .mat fields (matching the interface used by the existing pipeline):
%   Y0           : p x n   pre-sample (first p rows)
%   shortY       : T x n   estimation sample (T = total - p)
%   data_full    : (T+p) x n  full series in log-variance units
%   dates_full   : (T+p) x 1  MATLAB datenum
%   dates_shortY : T x 1      datenum for the estimation sample
%   labels       : n x 1 cell of bank tickers (CSV column names)
%   bank_info    : n x 1 struct with .ticker, .name, .country, .region
%   n, T, p, k   : dimensions (k = n*p)
%   minnesotaPriorMean : n x 1 vector of zeros (log-vols are stationary)
% =========================================================================

clear; clc;

%% ===================== Configuration ====================================
csv_path = 'data/ddly-data.csv';
n_keep   = 30;
p        = 3;
out_path = 'data/ddly_top30_data.mat';

%% ===================== Read raw CSV =====================================
fprintf('Reading %s ...\n', csv_path);
if ~exist(csv_path, 'file')
    error('load_ddly_data:FileNotFound', ...
        'Cannot find %s. Set csv_path to the correct location.', csv_path);
end

% The CSV uses semicolon delimiter, dd/mm/yy date format, and includes
% 96 bank columns + 10 government bond columns. We use readtable with
% explicit options to handle the dd/mm/yy date format reliably.
opts = detectImportOptions(csv_path, 'Delimiter', ';', 'NumHeaderLines', 0, ...
    'VariableNamingRule', 'preserve');
opts = setvartype(opts, opts.VariableNames(1), 'char');  % keep date as text
for j = 2:length(opts.VariableNames)
    opts = setvartype(opts, opts.VariableNames{j}, 'double');
end
T_table = readtable(csv_path, opts);

% Total: 1 date col + 96 banks + 10 bonds = 107 columns expected
all_var_names = T_table.Properties.VariableNames;
fprintf('  Read %d rows x %d cols\n', height(T_table), width(T_table));

if width(T_table) ~= 107
    warning('load_ddly_data:UnexpectedColCount', ...
        'Expected 107 columns (1 date + 96 banks + 10 bonds), got %d.', ...
        width(T_table));
end

%% ===================== Parse dates ======================================
date_strings = T_table{:, 1};
dates_full   = datenum(date_strings, 'dd/mm/yy');

% datenum interprets two-digit years using a pivot year. The DDLY sample
% spans 2003-2014, so all years should resolve to 2003-2014. Verify.
year_vals = year(datetime(dates_full, 'ConvertFrom', 'datenum'));
if any(year_vals < 2000) || any(year_vals > 2020)
    warning('load_ddly_data:DateParseProblem', ...
        'Some parsed years are outside [2000, 2020]. Check the date format.');
end

fprintf('  Date range: %s to %s\n', ...
    datestr(dates_full(1),  'yyyy-mm-dd'), ...
    datestr(dates_full(end),'yyyy-mm-dd'));

%% ===================== Extract top n_keep banks =========================
% Banks are columns 2..97 (col 1 = date, cols 98..107 = bonds).
% The CSV column order matches DDLY appendix Table A1 (descending by assets),
% so the top n_keep are simply the first n_keep bank columns.
bank_col_idx = 2 : (1 + n_keep);
bank_tickers = all_var_names(bank_col_idx)';   % n_keep x 1 cell
data_raw     = T_table{:, bank_col_idx};       % (T+p) x n_keep, GK variances

fprintf('  Selected top %d banks (columns 2..%d of CSV)\n', n_keep, 1+n_keep);

%% ===================== Data quality check ===============================
n_missing = sum(isnan(data_raw(:)));
n_zero    = sum(data_raw(:) == 0);
n_neg     = sum(data_raw(:) < 0);
fprintf('  Quality: %d missing, %d zero, %d negative out of %d total cells\n', ...
    n_missing, n_zero, n_neg, numel(data_raw));

if n_missing > 0
    error('load_ddly_data:MissingValues', ...
        'Found %d missing values in top-%d panel. Cannot proceed without imputation.', ...
        n_missing, n_keep);
end
if n_neg > 0
    error('load_ddly_data:NegativeValues', ...
        'Found %d negative values - GK variances must be non-negative.', n_neg);
end
if n_zero > 0
    warning('load_ddly_data:ZeroValues', ...
        'Found %d zero values; will floor at machine epsilon before log transform.', ...
        n_zero);
    data_raw(data_raw == 0) = eps;
end

%% ===================== Log-variance transform ===========================
fprintf('Applying log transform: y_{i,t} = log(GK_variance_{i,t})\n');
data_full = log(data_raw);

% Sanity check on transformed range
fprintf('  log-variance range: [%.2f, %.2f]\n', ...
    min(data_full(:)), max(data_full(:)));
fprintf('  cross-sectional mean of time-series std: %.3f\n', ...
    mean(std(data_full, 0, 1)));

%% ===================== Bank metadata (top 30 by assets) =================
% Hardcoded mapping from CSV ticker -> bank name, country, region.
% Source: DDLY (2018) appendix Table A1, descending by assets at end-2013.
% This is locked to n_keep = 30; if you change n_keep, extend this table.
bank_info_full = {
    'hsba.ln',    'HSBC Holdings',                      'UK',          'Europe';
    'X8306.to',   'Mitsubishi UFJ Financial Group',     'Japan',       'Asia';
    'bnp.fr',     'BNP Paribas',                        'France',      'Europe';
    'jpm',        'JPMorgan Chase & Co',                'US',          'N.America';
    'dbk.xe',     'Deutsche Bank',                      'Germany',     'Europe';
    'barc.ln',    'Barclays',                           'UK',          'Europe';
    'aca.fr',     'Credit Agricole',                    'France',      'Europe';
    'bac',        'Bank of America',                    'US',          'N.America';
    'c',          'Citigroup',                          'US',          'N.America';
    'X8411.to',   'Mizuho Financial Group',             'Japan',       'Asia';
    'gle.fr',     'Societe Generale',                   'France',      'Europe';
    'rbs.ln',     'Royal Bank of Scotland Group',       'UK',          'Europe';
    'X8316.to',   'Sumitomo Mitsui Financial Group',    'Japan',       'Asia';
    'san.mc',     'Banco Santander',                    'Spain',       'Europe';
    'wfc',        'Wells Fargo',                        'US',          'N.America';
    'inga.ae',    'ING Groep',                          'Netherlands', 'Europe';
    'lloy.ln',    'Lloyds Banking Group',               'UK',          'Europe';
    'ucg.mi',     'Unicredit',                          'Italy',       'Europe';
    'ubsn.vx',    'UBS',                                'Switzerland', 'Europe';
    'csgn.vx',    'Credit Suisse Group',                'Switzerland', 'Europe';
    'gs',         'Goldman Sachs Group',                'US',          'N.America';
    'ndasek.sk',  'Nordea Bank',                        'Sweden',      'Europe';
    'isp.mi',     'Intesa Sanpaolo',                    'Italy',       'Europe';
    'MS',         'Morgan Stanley',                     'US',          'N.America';
    'td.t',       'Toronto-Dominion Bank',              'Canada',      'N.America';
    'ry.t',       'Royal Bank of Canada',               'Canada',      'N.America';
    'bbva.mc',    'Banco Bilbao Vizcaya Argentaria',    'Spain',       'Europe';
    'cbk.xe',     'Commerzbank',                        'Germany',     'Europe';
    'nab.au',     'National Australia Bank',            'Australia',   'Oceania';
    'bns.t',      'Bank of Nova Scotia',                'Canada',      'N.America';
};

if n_keep ~= 30
    error('load_ddly_data:MetadataMismatch', ...
        'bank_info_full is hardcoded for n_keep = 30. Extend it for n_keep = %d.', ...
        n_keep);
end

% Verify ordering: ticker in metadata table must match CSV column ticker
for i = 1:n_keep
    csv_tk  = bank_tickers{i};
    meta_tk = bank_info_full{i, 1};
    if ~strcmp(csv_tk, meta_tk)
        error('load_ddly_data:TickerMismatch', ...
            'Position %d: CSV ticker "%s" does not match metadata "%s". ', ...
            i, csv_tk, meta_tk);
    end
end

% Build struct array
bank_info = struct('ticker', bank_info_full(:,1), ...
                   'name',   bank_info_full(:,2), ...
                   'country',bank_info_full(:,3), ...
                   'region', bank_info_full(:,4));

%% ===================== Construct Y0 / shortY interface ==================
% Match the convention used by the existing DPM-OISV and SV pipelines:
%   Y0     : p x n   pre-sample (initial conditions for lag construction)
%   shortY : T x n   estimation sample
n = n_keep;
Y0           = data_full(1:p, :);
shortY       = data_full(p+1:end, :);
T            = size(shortY, 1);
k            = n * p;
dates_shortY = dates_full(p+1:end);

labels = bank_tickers;             % n x 1 cell of CSV ticker strings

% Bank vol panel is stationary -> Minnesota prior mean is zero.
minnesotaPriorMean = zeros(n, 1);

%% ===================== Summary printout =================================
fprintf('\n----- Constructed dataset -----\n');
fprintf('  n = %d banks\n', n);
fprintf('  p = %d lags\n', p);
fprintf('  T = %d (estimation), T+p = %d (total)\n', T, T+p);
fprintf('  k = n*p = %d regressors per equation (no intercept by default)\n', k);
fprintf('  Estimation sample: %s to %s\n', ...
    datestr(dates_shortY(1),  'yyyy-mm-dd'), ...
    datestr(dates_shortY(end),'yyyy-mm-dd'));

% Regional breakdown
regions = {bank_info.region};
[uniq_reg, ~, ridx] = unique(regions);
fprintf('  Regional composition:\n');
for r = 1:length(uniq_reg)
    fprintf('    %-12s: %d banks\n', uniq_reg{r}, sum(ridx == r));
end

%% ===================== Save =============================================
save(out_path, ...
    'Y0', 'shortY', 'data_full', 'dates_full', 'dates_shortY', ...
    'labels', 'bank_info', 'n', 'T', 'p', 'k', 'minnesotaPriorMean', ...
    '-v7.3');
fprintf('\nSaved to %s\n', out_path);
