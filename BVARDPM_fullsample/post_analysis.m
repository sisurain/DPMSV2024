% =========================================================================
% post_analysis.m
%
% Post-estimation analysis for the DPM-OIMSV model full-sample results.
% Produces:
%   - Posterior means of lambda_t and h_t (saved to Excel)
%   - Posterior average number of clusters
%   - Asymmetry diagnostics for implied error distribution
%   - Convergence diagnostics (Geweke, inefficiency factors)
%   - Plots (volatility paths, lambda, cluster trace)
%
% Run AFTER the estimation script has saved results_DPM_OIMSV.mat.
% =========================================================================

clear; clc; close all;

%% ===================== User Configuration ===============================
run_OIMSV = true;    % Analyse DPM-OIMSV results?

%% ========================================================================
%  DPM-OIMSV posterior summaries
%% ========================================================================
if run_OIMSV
    oimsv_file = fullfile('DPM_OIMSV', 'results_DPM_OIMSV.mat');
    if ~exist(oimsv_file, 'file')
        warning('results_DPM_OIMSV.mat not found. Skipping OIMSV analysis.');
        run_OIMSV = false;
    end
end

if run_OIMSV
    fprintf('\n============================================================\n');
    fprintf('  POST-ANALYSIS: DPM-OIMSV\n');
    fprintf('============================================================\n');

    R = load(oimsv_file);

    % --- Posterior means ---
    H_mean_oimsv   = mean(R.store_H, 3);          % T x n
    lam_mean_oimsv = mean(R.store_lam, 2);         % T x 1
    B0_mean        = mean(R.store_B0, 3);          % n x n

    alpha_mean_oimsv = mean(R.store_DPM_params(:,1));
    K_mean_oimsv     = mean(R.store_DPM_params(:,2));

    SV_params_mean = mean(R.store_SV_params, 1);
    phi_mean   = SV_params_mean(1:R.n);
    omega2_mean = SV_params_mean(R.n+1:end);

    % --- Time axis ---
    if isfield(R, 'dates_shortY')
        T_id = datetime(datestr(R.dates_shortY));
    else
        T_id = (1:R.T)';
    end

    % --- Print summary ---
    fprintf('\nPosterior means:\n');
    fprintf('  alpha (concentration) = %.4f\n', alpha_mean_oimsv);
    fprintf('  K (clusters)          = %.2f\n', K_mean_oimsv);
    fprintf('  phi  (SV persistence) = %s\n', mat2str(round(phi_mean, 4)));
    fprintf('  omega2 (SV innov.)    = %s\n', mat2str(round(omega2_mean, 6)));

    fprintf('\nLambda statistics:\n');
    fprintf('  Min    = %.4f\n', min(lam_mean_oimsv));
    fprintf('  Median = %.4f\n', median(lam_mean_oimsv));
    fprintf('  Max    = %.4f\n', max(lam_mean_oimsv));

    % Outlier periods
    idx_large = find(lam_mean_oimsv > 20);
    if ~isempty(idx_large)
        fprintf('\nPeriods with lambda_mean > 20:\n');
        for ii = 1:length(idx_large)
            fprintf('  t=%d  (%s)  lambda=%.2f\n', idx_large(ii), ...
                    datestr(R.dates_shortY(idx_large(ii)),'yyyy-mmm'), ...
                    lam_mean_oimsv(idx_large(ii)));
        end
    end

    % --- Asymmetry diagnostic ---
    mu_mean_oimsv = mean(R.store_mu, 3);
    mu_overall    = mean(mu_mean_oimsv, 1);
    fprintf('\nAsymmetry diagnostic (posterior mean of mu across time):\n');
    fprintf('  Overall mean of mu: %s\n', mat2str(round(mu_overall, 4)));
    skew_mu = skewness(mu_mean_oimsv, 0, 1);
    fprintf('  Skewness of mu_t:   %s\n', mat2str(round(skew_mu, 4)));
    if any(abs(mu_overall) > 0.05)
        fprintf('  => Evidence of ASYMMETRY in the implied error distribution.\n');
    else
        fprintf('  => Implied error distribution appears roughly SYMMETRIC.\n');
    end

    % --- Save to Excel ---
    oimsv_excel = 'posterior_means_OIMSV.xlsx';
    % H sheet (individual SV paths)
    H_table_data = [R.dates_shortY, H_mean_oimsv];
    H_varnames = ['Date', arrayfun(@(i) sprintf('h_%d_t', i), 1:R.n, 'UniformOutput', false)];
    H_table = array2table(H_table_data, 'VariableNames', H_varnames);
    writetable(H_table, oimsv_excel, 'Sheet', 'h_it');

    % Lambda sheet
    L_table = table(R.dates_shortY, lam_mean_oimsv, ...
        'VariableNames', {'Date','lambda_t'});
    writetable(L_table, oimsv_excel, 'Sheet', 'lambda_t');
    fprintf('\nSaved posterior means to %s\n', oimsv_excel);

    % --- Convergence diagnostics ---
    fprintf('\n--- Convergence Diagnostics (DPM-OIMSV) ---\n');

    % Alpha and K
    run_convergence_diagnostics(R.store_DPM_params, {'alpha','K'}, 'OIMSV DPM');

    % SV parameters
    sv_labels = [arrayfun(@(i) sprintf('phi_%d',i), 1:R.n, 'UniformOutput', false), ...
                 arrayfun(@(i) sprintf('omega2_%d',i), 1:R.n, 'UniformOutput', false)];
    run_convergence_diagnostics(R.store_SV_params, sv_labels, 'OIMSV SV params');

    % Lambda diagnostics
    lam_rep_idx = round(linspace(1, R.T, min(5, R.T)));
    lam_rep_draws = R.store_lam(lam_rep_idx, :)';
    lam_labels = arrayfun(@(i) sprintf('lam_t=%d', i), lam_rep_idx, 'UniformOutput', false);
    run_convergence_diagnostics(lam_rep_draws, lam_labels, 'OIMSV lambda');

    % --- Plots ---
    % Individual SV paths (first 4 variables)
    figure('Name','DPM-OIMSV: Individual SV paths','Position',[100 400 900 600]);
    n_plot = min(4, R.n);
    for j = 1:n_plot
        subplot(2, 2, j);
        plot(T_id, exp(H_mean_oimsv(:,j)/2), 'b-', 'LineWidth', 1);
        if isfield(R,'labels')
            title(sprintf('$\\sqrt{\\exp(h_{%d,t})}$ (%s)', j, R.labels{j}), 'Interpreter','latex');
        else
            title(sprintf('$\\sqrt{\\exp(h_{%d,t})}$', j), 'Interpreter','latex');
        end
        xlim([T_id(1) T_id(end)]); grid on; box off;
    end

    figure('Name','DPM-OIMSV: Lambda','Position',[100 50 900 300]);
    lam_plot = min(lam_mean_oimsv, 30);
    plot(T_id, lam_plot, 'r-', 'LineWidth', 1.2);
    title('DPM-OIMSV: Posterior Mean of $\lambda_t$ (capped at 30)','Interpreter','latex');
    xlim([T_id(1) T_id(end)]); grid on; box off;

    figure('Name','DPM-OIMSV: K trace');
    plot(R.store_DPM_params(:,2));
    title('DPM-OIMSV: Trace of K (number of clusters)');
    xlabel('MCMC iteration'); ylabel('K');

    clear R;
end

%% ========================================================================
%  Figure 1: posterior mean of lambda_t with pointwise 90% credible band
%% ========================================================================
if run_OIMSV
    R = load(oimsv_file, 'store_lam', 'store_DPM_params', 'dates_shortY');

    lam_mean = mean(R.store_lam, 2);              % T x 1
    lam_band = prctile(R.store_lam, [5 95], 2);   % T x 2 (pointwise)
    dts      = R.dates_shortY(:);

    figure('Position', [100 100 900 400]);
    fill([dts; flipud(dts)], [lam_band(:,1); flipud(lam_band(:,2))], ...
         [0.75 0.75 0.75], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    hold on;
    plot(dts, lam_mean, 'k-', 'LineWidth', 1);
    hold off;
    datetick('x', 'yyyy', 'keeplimits');
    ylabel('\lambda_t', 'Interpreter', 'tex');
    box off;
    saveas(gcf, 'fig_lambda_fullsample.png');
    fprintf('Saved fig_lambda_fullsample.png\n');

    % Posterior statistics for the DPM concentration and cluster count
    alpha_draws = R.store_DPM_params(:,1);
    K_draws     = R.store_DPM_params(:,2);
    alpha_ci    = prctile(alpha_draws, [5 95]);
    fprintf('\nDPM concentration alpha: mean %.2f, 90%% CI (%.2f, %.2f)\n', ...
            mean(alpha_draws), alpha_ci(1), alpha_ci(2));
    fprintf('Active clusters K: min %d, max %d, mode %d\n', ...
            min(K_draws), max(K_draws), mode(K_draws));
    clear R;
end

fprintf('\nPost-analysis complete.\n');

%% ========================================================================
%  LOCAL FUNCTIONS
%% ========================================================================

function run_convergence_diagnostics(draws, param_names, group_name)
% DRAWS: nsims x nparams matrix
% PARAM_NAMES: cell array of parameter name strings
% GROUP_NAME: label for the group

    [nsims, nparams] = size(draws);
    fprintf('\n  %s  (nsims=%d)\n', group_name, nsims);
    fprintf('  %-16s %10s %10s %10s %10s\n', 'Parameter', 'Mean', 'IF', 'Geweke_z', 'Geweke_p');
    fprintf('  %s\n', repmat('-', 1, 60));

    for j = 1:nparams
        chain = draws(:, j);
        mn = mean(chain);

        % Inefficiency Factor (batch means)
        IF_val = compute_IF(chain);

        % Geweke diagnostic (first 10% vs last 50%)
        [z_stat, p_val] = geweke_test(chain, 0.1, 0.5);

        fprintf('  %-16s %10.4f %10.2f %10.3f %10.4f', ...
                param_names{j}, mn, IF_val, z_stat, p_val);
        if p_val < 0.05
            fprintf('  *');
        end
        fprintf('\n');
    end
end

function IF = compute_IF(chain)
% Inefficiency factor via batch means estimator
    N = length(chain);
    batch_size = max(10, floor(sqrt(N)));
    n_batches  = floor(N / batch_size);
    chain_trunc = chain(1:n_batches*batch_size);
    batches = reshape(chain_trunc, batch_size, n_batches);
    batch_means = mean(batches, 1);

    var_overall = var(chain);
    if var_overall < 1e-20
        IF = 1;
        return;
    end
    var_bm = var(batch_means) * batch_size;
    IF = var_bm / var_overall;
    IF = max(IF, 1);  % IF >= 1 by definition
end

function [z_stat, p_val] = geweke_test(chain, frac_a, frac_b)
% Geweke (1992) convergence diagnostic
    N = length(chain);
    na = floor(frac_a * N);
    nb = floor(frac_b * N);

    chain_a = chain(1:na);
    chain_b = chain(end-nb+1:end);

    mean_a = mean(chain_a);
    mean_b = mean(chain_b);

    % Spectral density at frequency zero (simple AR(1) approximation)
    var_a = spectral_var(chain_a);
    var_b = spectral_var(chain_b);

    z_stat = (mean_a - mean_b) / sqrt(var_a/na + var_b/nb);
    p_val  = 2 * (1 - normcdf(abs(z_stat)));
end

function sv = spectral_var(x)
% Simple spectral density at zero using AR(1) approximation
    x = x - mean(x);
    N = length(x);
    if N < 3
        sv = var(x);
        return;
    end
    rho1 = x(1:end-1)' * x(2:end) / (x' * x);
    rho1 = max(min(rho1, 0.999), -0.999);  % bound
    sv   = var(x) / (1 - rho1)^2;
end
