% =========================================================================
% run_DPMCSV_sim.m
%
% Runs the DPM-CSV model on all 6 simulation datasets and saves
% per-replication results to model_DPMCSV/results/.
%
% MODEL: BVAR-DPM-CSV
%   y_t = X_t * A + eps_t
%   eps_t ~ N(mu_t,  lambda_t * exp(h_t) * Sigma)
%   (mu_t, lambda_t) | z_t ~ DPM  (Dirichlet Process Mixture)
%   h_t = rho * h_{t-1} + eta_t,   eta_t ~ N(0, sigh2)  [scalar AR(1) CSV]
%   Sigma ~ IW(S0, nu0)
%   A ~ Minnesota prior (no intercept, k = n*p)
%
% OUTPUT per .mat file:
%   lam_hat    T x 1   posterior mean of lambda_t  (variance scale)
%   h_hat      T x 1   posterior mean of h_t       (mean-centred)
%   mu_hat     T x n   posterior mean of mu_t      (cluster mean)
%   lam_true   T x 1   true lambda_t from DGP
%   h_true     T x 1   true h_t from DGP
%   K_hat      scalar  posterior mean number of clusters
%   alpha_hat  scalar  posterior mean concentration parameter
%   rho_hat    scalar  posterior mean of rho
%   sigh2_hat  scalar  posterior mean of sigh2
%   dgp, T_sim, r : bookkeeping scalars
%
% FILES REQUIRED in model_DPMCSV/:
%   run_DPMCSV_sim.m          (this file)
%   construct_prior_A_func.m  (Minnesota prior construction)
%   nigSV.m                   (NIG base measure class)
%   oneDPM.m                  (collapsed Gibbs initialisation)
%   sample_h.m                (AR-MH sampler for stochastic volatility)
%   randmn.m                  (categorical distribution sampler)
% =========================================================================

clear; clc;

%% ---- 0. Paths ---------------------------------------------------------------
data_dir    = '../sim_data/';
results_dir = 'results/';
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
    fprintf('Created results directory: %s\n', results_dir);
end

%% ---- 1. MCMC settings -------------------------------------------------------
p      = 12;    % VAR lag order (monthly data)
nsims  = 2000;  % posterior draws to keep
burnin = 500;   % burn-in draws

include_intercept = false;  % no intercept (mu_t from DPM absorbs it)

%% ---- 2. Prior hyperparameters ------------------------------------------------

% IW prior on Sigma:
%   S0 = eye(n), nu0 = n+3   [set inside loop, n is data-dependent]

% SV process priors:
nuh0 = 1e4;               % IG shape for sigh2
Sh0  = 0.001 * (nuh0-1);  % IG scale for sigh2 => prior mean = 0.001
rho0 = 0.9;               % prior mean for rho
Vrho = 0.2^2;             % prior variance for rho (truncated Normal MH)

% DPM concentration alpha:
a0 = 2;  b0 = 4;   % Gamma sampler via auxiliary variable

% NIG base measure for DPM (nigSV parameters):
nu0_dp = 10;        % NIG shape
S0_dp  = 8;         % NIG scale

%% ---- 3. Main loop -----------------------------------------------------------
DGP_list = [1, 2, 3];
T_list   = [300, 800];

for dgp = DGP_list
    for T_sim = T_list

        %% Load simulation dataset
        fname = sprintf('%ssimdata_DGP%d_T%d_R20.mat', data_dir, dgp, T_sim);
        if ~exist(fname, 'file')
            warning('File not found, skipping: %s', fname);
            continue;
        end
        load(fname, 'all_Y', 'all_lam', 'all_h', 'Y0_fixed', 'R');

        n = size(all_Y{1}, 2);   % 16 variables
        k = n * p;               % no intercept

        % NIG base measure (size depends on n)
        mu0_dp = zeros(n, 1);
        V0_dp  = 100 * eye(n);

        fprintf('\n=== DGP%d  T=%d  n=%d  k=%d  (R=%d reps) ===\n', ...
                dgp, T_sim, n, k, R);

        for r = 1 : R
            fprintf('  Rep %2d/%d ... ', r, R);
            t_rep = tic;

            rng(1000 * dgp + T_sim * 100 + r);   % reproducible seed

            %% ----------------------------------------------------------------
            %% Data setup
            %% ----------------------------------------------------------------
            Y0     = Y0_fixed;        % p x n   pre-sample
            shortY = all_Y{r};        % T x n   estimation sample
            T      = T_sim;
            tmpY   = [Y0; shortY];    % (p+T) x n  for regressor construction

            %% Construct Minnesota prior
            [A0, VA0, ~] = construct_prior_A_func(Y0, shortY, p, include_intercept);

            %% IW prior on Sigma
            S0  = eye(n);
            nu0 = n + 3;

            %% Construct regressor matrix X (no intercept)
            X = zeros(T, k);
            for i = 1:p
                X(:, (i-1)*n+1 : i*n) = tmpY(p-i+1 : end-i, :);
            end

            %% Pre-compute prior terms
            iVA0_sparse  = spdiags(1./VA0, 0, k, k);
            prior_term_A = iVA0_sparse * A0;
            prior_term_S = S0 + A0' * iVA0_sparse * A0;

            %% Initialise DPM base measure
            prior_dp = nigSV(mu0_dp, V0_dp, nu0_dp, S0_dp);

            %% Storage arrays
            store_A_sum   = zeros(k, n);
            store_Sig_sum = zeros(n, n);
            store_h       = zeros(T, nsims);
            store_lam     = zeros(T, nsims);
            store_mu      = zeros(T, n, nsims);
            store_z       = zeros(T, nsims);
            store_theta   = zeros(nsims, 4);    % [alpha, K, rho, sigh2]

            counth = 0; countrho = 0;

            %% Chain initialisation
            h     = zeros(T, 1);
            rho   = 0.8;
            sigh2 = 0.01;

            lam = 1 ./ gamrnd(nu0_dp/2, 2/S0_dp, T, 1);
            mu  = zeros(T, n);
            alpha = 1;

            % Initial cluster assignments via one-pass collapsed Gibbs
            [z, Theta, nk] = oneDPM(shortY, h, S0, alpha, prior_dp);

            %% ================================================================
            %% MCMC loop
            %% ================================================================
            for isim = 1 : nsims + burnin

                %% Step 1: Sample Sigma and A
                iOm_diag = exp(-h) ./ lam;
                iOm      = spdiags(iOm_diag, 0, T, T);
                XiOm     = X' * iOm;
                Y_adj    = shortY - mu;

                KA   = iVA0_sparse + XiOm * X;
                Ahat = KA \ (prior_term_A + XiOm * Y_adj);
                Shat = prior_term_S + Y_adj' * iOm * Y_adj - Ahat' * KA * Ahat;
                Shat = (Shat + Shat') / 2;

                Sig  = iwishrnd(Shat, nu0 + T);
                CSig = chol(Sig, 'lower');
                A    = Ahat + (chol(KA, 'lower')' \ randn(k, n)) * CSig';

                %% Step 2: Sample h_t (scalar AR(1) CSV)
                U   = shortY - X * A - mu;
                tmp = U / CSig';
                s2  = sum(tmp.^2, 2) ./ lam;
                [h, flag] = sample_h(s2, rho, sigh2, h, n);
                counth = counth + flag;

                %% Step 3: Sample (lambda_t, mu_t) jointly via DPM Collapsed Gibbs
                X_dp = shortY - X * A;

                % Part 3a: Rebuild Theta from current z, h, Sig
                K = max(z);
                if isempty(K), K = 0; end

                Theta = cell(1, K);
                nk    = zeros(1, K);
                for j = 1:K
                    idxT  = find(z == j);
                    nk(j) = length(idxT);
                    if nk(j) > 0
                        Theta{j} = prior_dp.clone();
                        Theta{j}.addData(X_dp(idxT, :), h(idxT), Sig);
                    end
                end

                % Part 3b: Sequential Collapsed Gibbs sweep (random order)
                for i = randperm(T)
                    x  = X_dp(i, :);
                    hi = h(i);
                    kk = z(i);

                    % Remove obs i from its current cluster
                    Theta{kk}.delSample(x, hi, Sig);
                    nk(kk) = nk(kk) - 1;

                    % Log probabilities for existing clusters + new cluster
                    if ~isempty(nk)
                        Pk = log(nk) + cellfun(@(t) t.logPredPdf(x, hi, Sig), Theta);
                    else
                        Pk = [];
                    end
                    P0 = log(alpha) + prior_dp.logPredPdf(x, hi, Sig);
                    pp = [Pk, P0];

                    % Normalize (log-sum-exp) and sample
                    pp = pp - max(pp);
                    pp = exp(pp);
                    pp = pp ./ sum(pp);
                    kk = randmn(pp);

                    % Add obs i to selected cluster
                    if kk == numel(Theta) + 1    % new cluster
                        Theta{kk} = prior_dp.clone();
                        Theta{kk}.addSample(x, hi, Sig);
                        nk = [nk, 1]; %#ok<AGROW>
                    else                         % existing cluster
                        Theta{kk}.addSample(x, hi, Sig);
                        nk(kk) = nk(kk) + 1;
                    end
                    z(i) = kk;

                    % Remove empty clusters and relabel
                    idx0 = find(nk == 0);
                    if ~isempty(idx0)
                        Theta(idx0) = [];
                        nk(idx0)    = [];
                        for ii = idx0
                            z(z > ii) = z(z > ii) - 1;
                        end
                    end
                end

                % Part 3c: Sample cluster parameters (lambda_k, mu_k)
                K    = length(Theta);
                lamK = zeros(K, 1);
                muK  = zeros(K, n);

                for i = 1:K
                    Ui  = Theta{i}.U_;
                    iVi = Theta{i}.iV_;
                    mui = Theta{i}.mu_;
                    nui = Theta{i}.nu_;

                    % Marginal scale S = U - 0.5 * mu' * iV * mu
                    S_k = Ui - mui' * iVi * mui / 2;
                    if S_k <= 0
                        S_k = S0_dp; nui = nu0_dp;
                    end

                    % lambda_k ~ IG(nui, S_k)
                    lamK(i) = 1 / gamrnd(nui, 1/S_k);

                    % mu_k | lambda_k ~ N(mui, lambda_k * iVi^{-1})
                    [L, flag_chol] = chol(iVi / lamK(i), 'lower');
                    if flag_chol == 0
                        muK(i, :) = (mui + L' \ randn(n, 1))';
                    else
                        muK(i, :) = mui';
                    end
                end

                % Map cluster parameters back to observation dimension
                if iscolumn(z), z = z'; end
                zK  = sparse(1:T, z, ones(T,1), T, K);
                lam = full(zK * lamK);
                mu  = full(zK * muK);

                %% Step 4: Sample sigh2 | h, rho
                eh    = [h(1) * sqrt(1 - rho^2); h(2:end) - rho * h(1:end-1)];
                sigh2 = 1 / gamrnd(nuh0 + T/2, 1 / (Sh0 + sum(eh.^2)/2));

                %% Step 5: Sample rho via MH (truncated Normal proposal)
                Krho   = 1/Vrho + sum(h(1:T-1).^2) / sigh2;
                rhohat = Krho \ (rho0/Vrho + h(1:T-1)' * h(2:T) / sigh2);
                rhoc   = rhohat + sqrt(1/Krho) * randn;
                grho   = @(x) -0.5*log(sigh2./(1-x.^2)) - 0.5*(1-x.^2)/sigh2 * h(1)^2;
                if abs(rhoc) < 0.9999
                    alpMH = grho(rhoc) - grho(rho);
                    if alpMH > log(rand)
                        rho = rhoc;
                        countrho = countrho + 1;
                    end
                end

                %% Step 6: Sample DPM concentration alpha
                xi        = betarnd(alpha + 1, T);
                pi_alpha  = (a0 + K - 1) / ((a0 + K - 1) + T * (b0 - log(xi)));
                if rand < pi_alpha
                    alpha = gamrnd(a0 + K,     1 / (b0 - log(xi)));
                else
                    alpha = gamrnd(a0 + K - 1, 1 / (b0 - log(xi)));
                end

                %% Store post-burnin draws
                if isim > burnin
                    isave = isim - burnin;
                    store_A_sum   = store_A_sum + A;
                    store_Sig_sum = store_Sig_sum + Sig;
                    store_h(:, isave)        = h;
                    store_lam(:, isave)      = lam;
                    store_mu(:, :, isave)    = mu;
                    store_z(:, isave)        = z';
                    store_theta(isave, :)    = [alpha, K, rho, sigh2];
                end

                if mod(isim, 200) == 0
                    fprintf('\n    iter %d/%d (K=%d, alpha=%.1f)', ...
                            isim, nsims+burnin, K, alpha);
                end

            end  % MCMC loop

            %% ----------------------------------------------------------------
            %% Posterior summaries
            %% ----------------------------------------------------------------
            lam_hat   = mean(store_lam, 2);
            h_hat_raw = mean(store_h,   2);
            h_hat     = h_hat_raw - mean(h_hat_raw);
            mu_hat    = mean(store_mu, 3);

            theta_mean = mean(store_theta, 1);
            alpha_hat  = theta_mean(1);
            K_hat      = theta_mean(2);
            rho_hat    = theta_mean(3);
            sigh2_hat  = theta_mean(4);

            %% True values (mean-centre h to match h_hat scale)
            lam_true   = all_lam{r};
            h_true_raw = all_h{r};
            h_true     = h_true_raw - mean(h_true_raw);

            %% Save
            out_file = sprintf('%sDPMCSV_DGP%d_T%d_rep%02d.mat', ...
                                results_dir, dgp, T_sim, r);
            save(out_file, ...
                'lam_hat', 'h_hat', 'mu_hat', ...
                'lam_true', 'h_true', ...
                'K_hat', 'alpha_hat', 'rho_hat', 'sigh2_hat', ...
                'dgp', 'T_sim', 'r');

            fprintf('\n  done (%.1f sec) -> %s\n', toc(t_rep), out_file);

        end  % rep loop
    end  % T_sim loop
end  % dgp loop

fprintf('\n=== DPM-CSV simulation complete. Results in %s ===\n', results_dir);
