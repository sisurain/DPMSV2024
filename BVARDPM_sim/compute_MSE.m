% =========================================================================
% compute_MSE.m
%
% PURPOSE: Load simulation results from all three model subfolders and
%          compute MSE comparison tables for lambda_t, h_t, and tau_t
%          recovery across SVo, SVt, and DPM-CSV models.
%
% PREREQUISITES:
%   Results files in:
%     model_SVo/results/     SVo_DGP{d}_T{T}_rep{rr}.mat
%     model_SVt/results/     SVt_DGP{d}_T{T}_rep{rr}.mat
%     model_DPMCSV/results/  DPMCSV_DGP{d}_T{T}_rep{rr}.mat
%
%   Each .mat file must contain:
%     lam_hat   T x 1   posterior mean of lambda_t
%     h_hat     T x 1   posterior mean of h_t
%     lam_true  T x 1   true lambda_t from DGP
%     h_true    T x 1   true h_t from DGP
%
% IDENTIFICATION CORRECTION:
%   lam: divided by non-outlier-regime mean (DGP2,3) or overall median (DGP1)
%   h:   subtracted non-outlier-regime mean
%   tau: divided by non-outlier-regime mean of sqrt(lam)*exp(h/2)
%
% MSE DEFINITION:
%   Pooled MSE = (1/RT) sum_r sum_t (estimate - truth)^2
%   Regime-split MSE computed per replication before averaging.
%
% OUTPUT:
%   Tables 1-3: MSE (with SE) for lambda, h, tau
%   Tables 4-6: Relative MSE (normalised by DPM-CSV)
%   Time-profile regime-split summaries
%   MSE_results.mat saved in this folder
% =========================================================================

clear; clc;

%% ---- 0. Configuration ---------------------------------------------------
models   = {'SVo', 'SVt', 'DPMCSV'};
DGP_list = [1, 2, 3, 4];   % R2: DGP4 = pure Gaussian-SV. DGPs with no
                            % results present are dropped automatically below.
T_list   = [300, 800];
R        = 20;

DGP_labels = { ...
    'DGP1: fat-tail continuous mixture (favours SVt)', ...
    'DGP2: discrete outliers q=0.05   (favours SVo)', ...
    'DGP3: pandemic block t0=0.7T, L=6, lam=25 (favours DPM-CSV)', ...
    'DGP4: pure Gaussian-SV, lambda_t = 1 (no outliers; cost of over-generality)'};

% --- R2: auto-detect which DGPs have result files, drop the rest ----------
have_dgp = false(size(DGP_list));
for di = 1:numel(DGP_list)
    for m = 1:numel(models)
        probe = sprintf('model_%s/results/%s_DGP%d_T%d_rep01.mat', ...
                        models{m}, models{m}, DGP_list(di), T_list(1));
        if exist(probe, 'file'); have_dgp(di) = true; end
    end
end
if ~all(have_dgp)
    fprintf('NOTE: no results found for DGP(s) %s - excluded from tables.\n', ...
            num2str(DGP_list(~have_dgp)));
end
DGP_labels = DGP_labels(have_dgp);
DGP_list   = DGP_list(have_dgp);

n_models = numel(models);
n_dgp    = numel(DGP_list);
n_T      = numel(T_list);
T_max    = max(T_list);

%% ---- 1. Storage ---------------------------------------------------------
% Scalar MSE (1/RT pooled)
MSE_lam      = nan(n_models, n_dgp, n_T);
MSE_h        = nan(n_models, n_dgp, n_T);
MSE_tau      = nan(n_models, n_dgp, n_T);
MSE_prod     = nan(n_models, n_dgp, n_T);       % R2 (R5-S2): lam*exp(h) product
MSE_lam_reps = nan(n_models, n_dgp, n_T, R);
MSE_h_reps   = nan(n_models, n_dgp, n_T, R);
MSE_tau_reps = nan(n_models, n_dgp, n_T, R);
MSE_prod_reps = nan(n_models, n_dgp, n_T, R);   % R2 (R5-S2)
n_found      = zeros(n_models, n_dgp, n_T);

% Per-t squared errors for time-profile plots
sqerr_lam_t = nan(n_models, n_dgp, n_T, T_max, R);
sqerr_h_t   = nan(n_models, n_dgp, n_T, T_max, R);
sqerr_tau_t = nan(n_models, n_dgp, n_T, T_max, R);

% Regime-split MSE per replication
% Shape: (model, dgp, T_idx, rep, regime) where regime 1=non-outlier, 2=outlier
regime_mse_lam = nan(n_models, n_dgp, n_T, R, 2);
regime_mse_h   = nan(n_models, n_dgp, n_T, R, 2);
regime_mse_tau = nan(n_models, n_dgp, n_T, R, 2);

%% ---- 2. Load and compute ------------------------------------------------
fprintf('Loading and computing...\n\n');

for m = 1:n_models
    model_name = models{m};
    res_dir    = sprintf('model_%s/results/', model_name);

    for di = 1:n_dgp
        dgp = DGP_list(di);

        for ti = 1:n_T
            T_sim = T_list(ti);

            for r = 1:R
                fname = sprintf('%s%s_DGP%d_T%d_rep%02d.mat', ...
                    res_dir, model_name, dgp, T_sim, r);
                if ~exist(fname, 'file')
                    warning('MISSING: %s', fname); continue;
                end

                d = load(fname, 'lam_hat','h_hat','lam_true','h_true');
                if numel(d.lam_hat) ~= T_sim
                    warning('Dimension mismatch: %s', fname); continue;
                end

                % Aggregate n-column paths to scalar if needed
                if size(d.lam_hat,  2)>1, d.lam_hat  = mean(d.lam_hat,  2); end
                if size(d.h_hat,    2)>1, d.h_hat    = mean(d.h_hat,    2); end
                if size(d.lam_true, 2)>1, d.lam_true = mean(d.lam_true, 2); end
                if size(d.h_true,   2)>1, d.h_true   = mean(d.h_true,   2); end

                %% Identify non-outlier indices and outlier mask (per rep)
                switch dgp
                    case 1
                        % Normalise by overall median
                        % "Outlier" defined as top 25% of lam_true
                        pct75   = prctile(d.lam_true, 75);
                        is_out  = d.lam_true >= pct75;
                        no_idx  = ~is_out;
                        norm_fn = @(x) median(x);
                    case 2
                        % Baseline = lam_true < 1.5
                        is_out  = d.lam_true >= 1.5;
                        no_idx  = ~is_out;
                        if sum(no_idx) < 5, no_idx = true(T_sim,1); end
                        norm_fn = @(x) mean(x(no_idx));
                    case 3
                        % 6-observation block starting at t0 = floor(0.7*T_sim)
                        L       = 6;
                        t0      = floor(0.7 * T_sim);
                        t_end   = min(t0 + L - 1, T_sim);
                        is_out  = false(T_sim, 1);
                        is_out(t0:t_end) = true;
                        no_idx  = ~is_out;
                        norm_fn = @(x) mean(x(no_idx));
                    case 4
                        % R2: pure Gaussian-SV - no outliers by construction
                        is_out  = false(T_sim, 1);
                        no_idx  = true(T_sim, 1);
                        norm_fn = @(x) mean(x);
                end

                %% Normalise lambda (identification correction)
                lam_hat_lv  = norm_fn(d.lam_hat);
                lam_true_lv = norm_fn(d.lam_true);
                if lam_hat_lv  <= 0, lam_hat_lv  = 1; end
                if lam_true_lv <= 0, lam_true_lv = 1; end
                lam_hat_n  = d.lam_hat  / lam_hat_lv;
                lam_true_n = d.lam_true / lam_true_lv;

                %% Centre h on non-outlier regime mean
                h_hat_c  = d.h_hat  - mean(d.h_hat(no_idx));
                h_true_c = d.h_true - mean(d.h_true(no_idx));

                %% Combined scale tau = sqrt(lam)*exp(h/2)
                tau_hat_r  = sqrt(d.lam_hat)  .* exp(d.h_hat  / 2);
                tau_true_r = sqrt(d.lam_true) .* exp(d.h_true / 2);
                tau_hat_lv  = mean(tau_hat_r(no_idx));
                tau_true_lv = mean(tau_true_r(no_idx));
                if tau_hat_lv  <= 0, tau_hat_lv  = 1; end
                if tau_true_lv <= 0, tau_true_lv = 1; end
                tau_hat_n  = tau_hat_r  / tau_hat_lv;
                tau_true_n = tau_true_r / tau_true_lv;

                %% R2 (R5-S2): tracking metric for the product lam_t*exp(h_t)
                % (the total variance multiplier; normalised on the
                % non-outlier regime like tau)
                prod_hat_r  = d.lam_hat  .* exp(d.h_hat);
                prod_true_r = d.lam_true .* exp(d.h_true);
                prod_hat_lv  = mean(prod_hat_r(no_idx));
                prod_true_lv = mean(prod_true_r(no_idx));
                if prod_hat_lv  <= 0, prod_hat_lv  = 1; end
                if prod_true_lv <= 0, prod_true_lv = 1; end
                prod_hat_n  = prod_hat_r  / prod_hat_lv;
                prod_true_n = prod_true_r / prod_true_lv;

                %% Per-t squared errors
                se_lam = (lam_hat_n - lam_true_n).^2;
                se_h   = (h_hat_c   - h_true_c  ).^2;
                se_tau = (tau_hat_n - tau_true_n).^2;
                se_prod = (prod_hat_n - prod_true_n).^2;   % R2 (R5-S2)

                %% Scalar MSE for this rep
                MSE_lam_reps(m,di,ti,r) = mean(se_lam);
                MSE_h_reps(m,di,ti,r)   = mean(se_h);
                MSE_tau_reps(m,di,ti,r) = mean(se_tau);
                MSE_prod_reps(m,di,ti,r) = mean(se_prod);  % R2 (R5-S2)
                n_found(m,di,ti)        = n_found(m,di,ti) + 1;

                %% Per-t storage (for time-axis plots)
                sqerr_lam_t(m,di,ti,1:T_sim,r) = se_lam;
                sqerr_h_t(m,di,ti,1:T_sim,r)   = se_h;
                sqerr_tau_t(m,di,ti,1:T_sim,r) = se_tau;

                %% Regime-split MSE (computed per replication)
                % Non-outlier regime (index 1)
                regime_mse_lam(m,di,ti,r,1) = mean(se_lam(no_idx),  'omitnan');
                regime_mse_h(m,di,ti,r,1)   = mean(se_h(no_idx),    'omitnan');
                regime_mse_tau(m,di,ti,r,1) = mean(se_tau(no_idx),  'omitnan');
                % Outlier regime (index 2)
                if any(is_out)
                    regime_mse_lam(m,di,ti,r,2) = mean(se_lam(is_out), 'omitnan');
                    regime_mse_h(m,di,ti,r,2)   = mean(se_h(is_out),   'omitnan');
                    regime_mse_tau(m,di,ti,r,2) = mean(se_tau(is_out), 'omitnan');
                end
            end

            %% Pooled MSE = (1/RT) sum_r sum_t
            MSE_lam(m,di,ti) = mean(MSE_lam_reps(m,di,ti,:), 4, 'omitnan');
            MSE_h(m,di,ti)   = mean(MSE_h_reps(m,di,ti,:),   4, 'omitnan');
            MSE_tau(m,di,ti) = mean(MSE_tau_reps(m,di,ti,:), 4, 'omitnan');
            MSE_prod(m,di,ti) = mean(MSE_prod_reps(m,di,ti,:), 4, 'omitnan');

            fprintf('  %-8s DGP%d T=%3d n=%2d  MSE_lam=%7.4f  MSE_h=%7.4f  MSE_tau=%7.4f\n', ...
                model_name, dgp, T_sim, n_found(m,di,ti), ...
                MSE_lam(m,di,ti), MSE_h(m,di,ti), MSE_tau(m,di,ti));
        end
    end
    fprintf('\n');
end

%% ---- 3. Standard errors -------------------------------------------------
SE_lam = nan(n_models,n_dgp,n_T);
SE_h   = nan(n_models,n_dgp,n_T);
SE_tau = nan(n_models,n_dgp,n_T);
SE_prod = nan(n_models,n_dgp,n_T);
for m=1:n_models; for di=1:n_dgp; for ti=1:n_T
    Ra = n_found(m,di,ti);
    if Ra > 1
        SE_lam(m,di,ti) = std(MSE_lam_reps(m,di,ti,:),0,4,'omitnan')/sqrt(Ra);
        SE_h(m,di,ti)   = std(MSE_h_reps(m,di,ti,:),  0,4,'omitnan')/sqrt(Ra);
        SE_tau(m,di,ti) = std(MSE_tau_reps(m,di,ti,:),0,4,'omitnan')/sqrt(Ra);
        SE_prod(m,di,ti) = std(MSE_prod_reps(m,di,ti,:),0,4,'omitnan')/sqrt(Ra);
    end
end; end; end

%% ---- 4. Time-axis MSEt (averaged over reps) -----------------------------
MSEt_lam = mean(sqerr_lam_t, 5, 'omitnan');   % (model, dgp, T_idx, t)
MSEt_h   = mean(sqerr_h_t,   5, 'omitnan');
MSEt_tau = mean(sqerr_tau_t, 5, 'omitnan');

%% ---- 5. Build flat matrices (DGP-outer / T-inner) -----------------------
n_col      = n_dgp * n_T;
col_labels = cell(1, n_col);
c = 0;
for di = 1:n_dgp
    for ti = 1:n_T
        c = c+1;
        col_labels{c} = sprintf('DGP%d/T=%d', DGP_list(di), T_list(ti));
    end
end

MSE_lam_flat = nan(n_models,n_col); MSE_h_flat = nan(n_models,n_col);
MSE_tau_flat = nan(n_models,n_col); MSE_prod_flat = nan(n_models,n_col);
SE_lam_flat  = nan(n_models,n_col); SE_h_flat  = nan(n_models,n_col);
SE_tau_flat  = nan(n_models,n_col); SE_prod_flat  = nan(n_models,n_col);
c = 0;
for di=1:n_dgp; for ti=1:n_T
    c = c+1;
    MSE_lam_flat(:,c) = MSE_lam(:,di,ti);
    MSE_h_flat(:,c)   = MSE_h(:,di,ti);
    MSE_tau_flat(:,c) = MSE_tau(:,di,ti);
    MSE_prod_flat(:,c) = MSE_prod(:,di,ti);
    SE_lam_flat(:,c)  = SE_lam(:,di,ti);
    SE_h_flat(:,c)    = SE_h(:,di,ti);
    SE_tau_flat(:,c)  = SE_tau(:,di,ti);
    SE_prod_flat(:,c)  = SE_prod(:,di,ti);
end; end

hdr_w = 10; col_w = 24;
sep   = repmat('=', 1, hdr_w + col_w*n_col);

%% ---- 6. Print MSE tables ------------------------------------------------
tbl_specs = { ...
    MSE_lam_flat, SE_lam_flat, ...
        'TABLE 1: MSE_lambda (normalised to non-outlier baseline)'; ...
    MSE_h_flat,   SE_h_flat, ...
        'TABLE 2: MSE_h (non-outlier regime centred)'; ...
    MSE_tau_flat, SE_tau_flat, ...
        'TABLE 3: MSE_tau (tau = sqrt(lam)*exp(h/2), combined scale)'; ...
    MSE_prod_flat, SE_prod_flat, ...
        'TABLE 3b (R5-S2): MSE of the product lam_t*exp(h_t) (total variance multiplier)'};

for tbl = 1:4
    dat=tbl_specs{tbl,1}; se=tbl_specs{tbl,2}; ttl=tbl_specs{tbl,3};
    fprintf('\n%s\n%s\n%s\n', sep, ttl, sep);
    fprintf('%-*s', hdr_w, 'Model');
    for c=1:n_col, fprintf('  %-*s', col_w-2, col_labels{c}); end
    fprintf('\n%s\n', repmat('-',1,hdr_w+col_w*n_col));
    for m=1:n_models
        fprintf('%-*s', hdr_w, models{m});
        for c=1:n_col
            fprintf('  %7.4f (±%6.4f)  ', dat(m,c), se(m,c));
        end
        fprintf('\n');
    end
end

%% ---- 7. Relative MSE tables (vs DPM-CSV) --------------------------------
dpmcsv_idx = find(strcmp(models,'DPMCSV'));
rel_specs = { ...
    MSE_lam_flat, 'TABLE 4: Relative MSE_lambda (vs DPM-CSV; >1 = worse than DPM-CSV)'; ...
    MSE_h_flat,   'TABLE 5: Relative MSE_h      (vs DPM-CSV; >1 = worse than DPM-CSV)'; ...
    MSE_tau_flat, 'TABLE 6: Relative MSE_tau     (vs DPM-CSV; >1 = worse than DPM-CSV)'; ...
    MSE_prod_flat, 'TABLE 6b (R5-S2): Relative MSE of lam*exp(h) (vs DPM-CSV)'};

for tbl = 1:4
    raw = rel_specs{tbl,1};
    dat = raw ./ raw(dpmcsv_idx,:);
    ttl = rel_specs{tbl,2};
    fprintf('\n%s\n%s\n%s\n', sep, ttl, sep);
    fprintf('%-*s', hdr_w, 'Model');
    for c=1:n_col, fprintf('  %-*s', col_w-2, col_labels{c}); end
    fprintf('\n%s\n', repmat('-',1,hdr_w+col_w*n_col));
    for m=1:n_models
        fprintf('%-*s', hdr_w, models{m});
        for c=1:n_col, fprintf('  %11.3f           ', dat(m,c)); end
        fprintf('\n');
    end
end

%% ---- 8. Time-profile regime-split table ---------------------------------
fprintf('\n%s\n', repmat('=',1,74));
fprintf('TIME-PROFILE: MSEt in non-outlier vs outlier regime\n');
fprintf('(averaged over replications; outlier mask applied per replication)\n');
fprintf('%s\n', repmat('=',1,74));

for di = 1:n_dgp
    dgp = DGP_list(di);
    for ti = 1:n_T
        T_sim = T_list(ti);

        switch dgp
            case 1
                regime_note = 'non-outlier: lower 75% lam_true | high-lam: top 25% lam_true';
            case 2
                regime_note = 'non-outlier: lam_true<1.5 | outlier: lam_true>=1.5 (~5% obs each rep)';
            case 3
                L_blk = 6;
                t0    = floor(0.7*T_sim);
                t_end = min(t0 + L_blk - 1, T_sim);
                regime_note = sprintf('non-outlier: t=1..%d & t=%d..%d | block: t=%d..%d (%d obs)', ...
                    t0-1, t_end+1, T_sim, t0, t_end, L_blk);
        end

        fprintf('\n  DGP%d / T=%d  [%s]\n', dgp, T_sim, regime_note);
        fprintf('  %-8s  %-5s  %-22s  %-22s\n', ...
            'Model','Var','MSEt (non-outlier)','MSEt (outlier/block)');
        fprintf('  %s\n', repmat('-',1,62));

        for m = 1:n_models
            no_lam  = mean(regime_mse_lam(m,di,ti,:,1), 4, 'omitnan');
            out_lam = mean(regime_mse_lam(m,di,ti,:,2), 4, 'omitnan');
            no_h    = mean(regime_mse_h(m,di,ti,:,1),   4, 'omitnan');
            out_h   = mean(regime_mse_h(m,di,ti,:,2),   4, 'omitnan');
            no_tau  = mean(regime_mse_tau(m,di,ti,:,1), 4, 'omitnan');
            out_tau = mean(regime_mse_tau(m,di,ti,:,2), 4, 'omitnan');

            fprintf('  %-8s  lam  %12.4f          %12.4f\n', models{m}, no_lam, out_lam);
            fprintf('  %-8s  h    %12.4f          %12.4f\n', '', no_h,   out_h);
            fprintf('  %-8s  tau  %12.4f          %12.4f\n', '', no_tau, out_tau);
        end
    end
end

%% ---- 9. Notes -----------------------------------------------------------
fprintf('\n--- DGP key ---\n');
for di=1:n_dgp, fprintf('  DGP%d: %s\n', di, DGP_labels{di}); end
fprintf('\n--- Normalisation ---\n');
fprintf('  lam: divided by non-outlier-regime mean (DGP2,3) or overall median (DGP1)\n');
fprintf('  h:   subtracted non-outlier-regime mean\n');
fprintf('  tau: divided by non-outlier-regime mean of sqrt(lam)*exp(h/2)\n');
fprintf('\n--- Regime-split approach ---\n');
fprintf('  DGP1: "outlier" = top 25%% of lam_true in each replication\n');
fprintf('  DGP2: "outlier" = t where lam_true>=1.5 in each replication (~5%% of obs)\n');
fprintf('  DGP3: "outlier" = 6-obs block starting at t = floor(0.7*T)\n');
fprintf('  Regime split applied per replication before averaging\n\n');

%% ---- 10. Save ----------------------------------------------------------
save('MSE_results.mat', ...
    'MSE_lam_flat','MSE_h_flat','MSE_tau_flat','MSE_prod_flat', ...
    'SE_lam_flat','SE_h_flat','SE_tau_flat','SE_prod_flat', ...
    'MSE_lam_reps','MSE_h_reps','MSE_tau_reps','MSE_prod_reps', ...
    'regime_mse_lam','regime_mse_h','regime_mse_tau', ...
    'MSEt_lam','MSEt_h','MSEt_tau', ...
    'n_found','models','DGP_list','T_list','col_labels');

fprintf('Results saved to MSE_results.mat\n');
