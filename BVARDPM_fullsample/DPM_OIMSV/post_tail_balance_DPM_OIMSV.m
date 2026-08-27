% =========================================================================
% post_tail_balance_DPM_OIMSV.m
%
% Post-estimation tail-balance diagnostic for the DPM-OIMSV model.
% Implements the procedure outlined in Section 5.1 (pages 22-23 of draft).
%
% Computes the asymmetry index A_i = log(p_i^+ / p_i^-) for each variable,
% where p_i^+ and p_i^- are the time-averaged upper and lower tail
% exceedance probabilities of the standardized predictive innovations.
%
% Required: results_DPM_OIMSV.mat  (from BVAR_DPM_OIMSV_FE_01.m)
%           CCMM16_data.mat         (or the original data file)
%
% =========================================================================

clear; clc; close all;

%% ===================== Configuration ====================================
q = 2;          % Tail threshold (number of std devs)
thin = 1;       % Thinning factor (use every thin-th draw; set 1 for all)

%% ===================== Load MCMC Results ================================
if exist('results_DPM_OIMSV.mat','file')
    load('results_DPM_OIMSV.mat');
    fprintf('Loaded results_DPM_OIMSV.mat: n=%d, T=%d, nsims=%d\n', n, T, nsims);
else
    error('results_DPM_OIMSV.mat not found.');
end

%% ===================== Reload Data & Reconstruct X ======================
% NOTE: shortY and X are NOT saved in results_DPM_OIMSV.mat.
%       We reconstruct them here from the original data.
%       Consider adding shortY to the save block in the main script.

if exist('CCMM16_data.mat','file')
    tmp = load('CCMM16_data.mat', 'Y0', 'shortY');
    Y0 = tmp.Y0; shortY = tmp.shortY;
    fprintf('Loaded data from CCMM16_data.mat\n');
else
    % Fallback: load from text file
    data_file = '2022-12-ccmm16-r.txt';
    if ~exist(data_file, 'file')
        error('Cannot find data. Place CCMM16_data.mat or %s in the working directory.', data_file);
    end
    data = readmatrix(data_file);
    Y0 = data(1:p,:);
    shortY = data(p+1:end,:);
    fprintf('Loaded data from %s\n', data_file);
end

tmpY = [Y0; shortY];
X = zeros(T, k);
for i = 1:p
    X(:, (i-1)*n+1 : i*n) = tmpY(p-i+1:end-i, :);
end

%% ===================== Select Draws (with optional thinning) ============
draw_idx = 1:thin:nsims;
S = length(draw_idx);
fprintf('Using S = %d posterior draws (thin = %d)\n', S, thin);

%% ===================== Step 1: Compute Innovations ======================
% For each draw s:  epsilon_t^(s) = y_t - X_t * A^(s) - mu_t^(s)
%
% Store as T x n x S array

fprintf('Computing posterior predictive innovations...\n');
eps_draws = zeros(T, n, S);

for ss = 1:S
    s = draw_idx(ss);
    A_s  = store_A(:,:,s);      % k x n
    mu_s = store_mu(:,:,s);     % T x n
    eps_draws(:,:,ss) = shortY - X*A_s - mu_s;
end

%% ===================== Step 2: Compute sigma_{i,t} ======================
% sigma_{i,t} = sqrt( E_s[ lambda_t^(s) * Sigma_{ii,t}^(s) ] )
% where Sigma_t = B0^{-1} * D_t * (B0^{-1})'
%       D_t = diag(exp(h_{1,t}), ..., exp(h_{n,t}))
%
% The (i,i) element: Sigma_{ii,t} = sum_j (iB0(i,j))^2 * exp(h_{j,t})

fprintf('Computing scaling factors sigma_{i,t}...\n');
lam_Sig_ii = zeros(T, n);   % Will accumulate E[lambda_t * Sigma_{ii,t}]

for ss = 1:S
    s = draw_idx(ss);
    
    B0_s  = store_B0(:,:,s);     % n x n
    H_s   = store_H(:,:,s);      % T x n
    lam_s = store_lam(:,s);      % T x 1
    
    % Invert B0
    iB0_s = B0_s \ eye(n);       % n x n
    
    % Compute Sigma_{ii,t} for all i and t simultaneously
    % Sigma_{ii,t} = sum_j iB0(i,j)^2 * exp(h_{j,t})
    % = (iB0.^2) * exp(H)' evaluated per time
    iB0_sq = iB0_s.^2;           % n x n: element (i,j) = iB0(i,j)^2
    expH   = exp(H_s);           % T x n
    
    % Sig_ii: T x n matrix where entry (t,i) = Sigma_{ii,t}
    Sig_ii = expH * iB0_sq';     % (T x n) * (n x n)' = T x n
    
    % Accumulate lambda_t * Sigma_{ii,t}
    lam_Sig_ii = lam_Sig_ii + lam_s .* Sig_ii;
end

% Posterior mean
lam_Sig_ii = lam_Sig_ii / S;

% sigma_{i,t} = sqrt(E[lambda_t * Sigma_{ii,t}])
sigma_it = sqrt(lam_Sig_ii);

% Safety: floor at machine epsilon to avoid division by zero
sigma_it = max(sigma_it, 1e-15);

%% ===================== Step 3: Standardize & Compute Exceedances ========
% tilde_eps_{i,t}^(s) = eps_{i,t}^(s) / sigma_{i,t}
% p_{i,t}^+ = (1/S) sum_s 1( tilde_eps > q )
% p_{i,t}^- = (1/S) sum_s 1( tilde_eps < -q )

fprintf('Computing tail exceedance probabilities (q = %.1f)...\n', q);

count_upper = zeros(T, n);   % count of tilde_eps > q
count_lower = zeros(T, n);   % count of tilde_eps < -q

for ss = 1:S
    eps_std = eps_draws(:,:,ss) ./ sigma_it;   % T x n
    count_upper = count_upper + (eps_std > q);
    count_lower = count_lower + (eps_std < -q);
end

p_upper_it = count_upper / S;   % p_{i,t}^+  (T x n)
p_lower_it = count_lower / S;   % p_{i,t}^-  (T x n)

%% ===================== Step 4: Time-Average & Asymmetry Index ===========
% p_i^+ = (1/T) sum_t p_{i,t}^+
% p_i^- = (1/T) sum_t p_{i,t}^-
% A_i   = log(p_i^+ / p_i^-)

p_upper = mean(p_upper_it, 1)';   % n x 1
p_lower = mean(p_lower_it, 1)';   % n x 1

% Robustness: if either tail probability is exactly zero, 
% set a floor at 1/(S*T) to avoid log(0) or log(0/0)
floor_prob = 1 / (S * T);
p_upper_safe = max(p_upper, floor_prob);
p_lower_safe = max(p_lower, floor_prob);

A_index = log(p_upper_safe ./ p_lower_safe);

%% ===================== Display Results & Save to Text File ===============

% Build variable name list
vnames = cell(n,1);
for i = 1:n
    if exist('labels','var') && length(labels) >= i
        vnames{i} = labels{i};
    else
        vnames{i} = sprintf('Var_%d', i);
    end
end

% Print to console AND write to text file simultaneously
fid = fopen('tail_balance_table.txt', 'w');
targets = {1, fid};   % 1 = stdout, fid = file

for tt = 1:2
    f = targets{tt};
    fprintf(f, '\n');
    fprintf(f, '=================================================================\n');
    fprintf(f, '  Tail-Balance Diagnostic (q = %.1f, S = %d draws)\n', q, S);
    fprintf(f, '=================================================================\n');
    fprintf(f, '%-25s  %8s  %8s  %8s\n', 'Variable', 'p_i^+', 'p_i^-', 'A_i');
    fprintf(f, '-----------------------------------------------------------------\n');
    for i = 1:n
        fprintf(f, '%-25s  %8.4f  %8.4f  %8.4f\n', vnames{i}, p_upper(i), p_lower(i), A_index(i));
    end
    fprintf(f, '-----------------------------------------------------------------\n');
    fprintf(f, 'A_i > 0: heavier right tail;  A_i < 0: heavier left tail\n');
    fprintf(f, 'A_i ~ 0: approximately symmetric tails\n');
    fprintf(f, '=================================================================\n');
end
fclose(fid);
fprintf('Table saved to tail_balance_table.txt\n\n');

%% ===================== Bar Chart of A_i =================================
figure('Position', [100 100 900 500]);

bar_colors = zeros(n, 3);
bar_colors(A_index > 0, :) = repmat([0.2 0.4 0.8], sum(A_index > 0), 1);  % blue for right-heavy
bar_colors(A_index < 0, :) = repmat([0.8 0.2 0.2], sum(A_index < 0), 1);  % red for left-heavy
bar_colors(A_index == 0, :) = repmat([0.5 0.5 0.5], sum(A_index == 0), 1); % grey

b = bar(A_index, 'FaceColor', 'flat');
b.CData = bar_colors;

hold on;
yline(0, 'k--', 'LineWidth', 1);
hold off;

if exist('labels','var') && length(labels) >= n
    set(gca, 'XTick', 1:n, 'XTickLabel', labels, 'XTickLabelRotation', 45);
else
    set(gca, 'XTick', 1:n);
end

ylabel('Asymmetry Index A_i');
title(sprintf('Tail-Balance Diagnostic (q = %.0f)', q));
set(gca, 'FontSize', 11);
grid on;

% Save figure
print('-dpdf', 'tail_balance_Ai.pdf');
print('-djpeg', '-r300', 'tail_balance_Ai.jpg');
fprintf('Figure saved to tail_balance_Ai.pdf and tail_balance_Ai.jpg\n');

%% ===================== Time-Varying Tail Probabilities (all variables) ===
% Plot p_{i,t}^+ and p_{i,t}^- for ALL n variables, 4 per page.

vars_per_page = 4;
n_pages = ceil(n / vars_per_page);

for pg = 1:n_pages
    figure('Position', [100 100 1000 700]);
    
    idx_start = (pg-1)*vars_per_page + 1;
    idx_end   = min(pg*vars_per_page, n);
    n_this    = idx_end - idx_start + 1;
    
    for jj = 1:n_this
        ii = idx_start + jj - 1;
        subplot(vars_per_page, 1, jj);
        
        if exist('dates_shortY','var') && length(dates_shortY) == T
            plot(dates_shortY, p_upper_it(:,ii), 'b-', 'LineWidth', 1); hold on;
            plot(dates_shortY, p_lower_it(:,ii), 'r-', 'LineWidth', 1);
            datetick('x', 'yyyy', 'keeplimits');
        else
            plot(1:T, p_upper_it(:,ii), 'b-', 'LineWidth', 1); hold on;
            plot(1:T, p_lower_it(:,ii), 'r-', 'LineWidth', 1);
        end
        hold off;
        
        title(sprintf('%s  (A_i = %.3f)', vnames{ii}, A_index(ii)), 'FontSize', 11);
        
        if jj == 1
            legend('p_{i,t}^+ (upper)', 'p_{i,t}^- (lower)', ...
                   'Location', 'northeast', 'FontSize', 9);
        end
        ylabel('Prob');
        grid on;
    end
    
    sgtitle(sprintf('Time-Varying Tail Probabilities (page %d/%d)', pg, n_pages), ...
            'FontSize', 13, 'FontWeight', 'bold');
    
    fname = sprintf('tail_balance_time_varying_p%d', pg);
    print('-dpdf',  [fname '.pdf']);
    print('-djpeg', '-r300', [fname '.jpg']);
    fprintf('Saved %s.pdf\n', fname);
end

%% ===================== Save Time-Varying Data to CSV ====================
% Two CSV files: one for p_{i,t}^+ and one for p_{i,t}^-
% Each has T rows and n+1 columns (time index or date + 16 variables)

% --- p_upper CSV ---
fid = fopen('tail_prob_upper.csv', 'w');
% Header
fprintf(fid, 't');
for i = 1:n; fprintf(fid, ',%s', vnames{i}); end
fprintf(fid, '\n');
% Data
for t = 1:T
    if exist('dates_shortY','var') && length(dates_shortY) == T
        fprintf(fid, '%.2f', dates_shortY(t));  % datenum or year-fraction
    else
        fprintf(fid, '%d', t);
    end
    for i = 1:n
        fprintf(fid, ',%.6f', p_upper_it(t,i));
    end
    fprintf(fid, '\n');
end
fclose(fid);

% --- p_lower CSV ---
fid = fopen('tail_prob_lower.csv', 'w');
fprintf(fid, 't');
for i = 1:n; fprintf(fid, ',%s', vnames{i}); end
fprintf(fid, '\n');
for t = 1:T
    if exist('dates_shortY','var') && length(dates_shortY) == T
        fprintf(fid, '%.2f', dates_shortY(t));
    else
        fprintf(fid, '%d', t);
    end
    for i = 1:n
        fprintf(fid, ',%.6f', p_lower_it(t,i));
    end
    fprintf(fid, '\n');
end
fclose(fid);

fprintf('Time-varying data saved to tail_prob_upper.csv and tail_prob_lower.csv\n');

%% ===================== Save All Diagnostics to .mat =====================
save('tail_balance_diagnostics.mat', ...
     'A_index', 'p_upper', 'p_lower', ...
     'p_upper_it', 'p_lower_it', ...
     'sigma_it', 'q', 'S', ...
     'vnames', 'n', 'T');
fprintf('Full diagnostics saved to tail_balance_diagnostics.mat\n');
fprintf('  - p_upper_it, p_lower_it: T x n matrices (all 16 variables)\n');
fprintf('  - A_index, p_upper, p_lower: n x 1 summary vectors\n');
fprintf('  - sigma_it: T x n scaling factors\n');
