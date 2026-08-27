% =========================================================================
% run_SVt_sim.m
%
% Runs the BVAR-t-CSV model (Chan, 2020) on all 6 simulation datasets
% and saves per-replication results to model_SVt/results/.
%
% MODEL: BVAR-t-CSV
%   y_t = X_t * A + eps_t
%   eps_t ~ N(0, lambda_t * exp(h_t) * Sigma)
%   lambda_t ~ IG(nu/2, nu/2)     i.i.d. fat-tail scaling
%   h_t = rho * h_{t-1} + eta_t,  eta_t ~ N(0, sigh2)   AR(1) log-vol
%   Sigma ~ IW(S0, nu0)
%   A ~ Minnesota prior (with intercept, k = n*p + 1)
%   nu ~ MH sampler (sample_nu.m)
%
% OUTPUT per .mat file:
%   lam_hat    T x 1   posterior mean of lambda_t  (variance scale)
%   h_hat      T x 1   posterior mean of h_t       (mean-centred)
%   lam_true   T x 1   true lambda_t from DGP
%   h_true     T x 1   true h_t from DGP          (mean-centred)
%   nu_hat     scalar  posterior mean of nu
%   rho_hat    scalar  posterior mean of rho
%   sigh2_hat  scalar  posterior mean of sigh2
%   dgp, T_sim, r : bookkeeping scalars
%
% FILES REQUIRED in model_SVt/:
%   run_SVt_sim.m    (this file)
%   sample_h.m       (AR-MH sampler for stochastic volatility)
%   sample_nu.m      (MH sampler for degrees of freedom nu)
% =========================================================================

clear; clc;

%% ---- 0. Paths ----------------------------------------------------------------
data_dir    = '../sim_data/';
results_dir = 'results/';
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
    fprintf('Created results directory: %s\n', results_dir);
end

%% ---- 1. MCMC settings -------------------------------------------------------
p      = 12;    % VAR lag order (monthly data)
nsims  = 500;   % posterior draws to keep
burnin = 100;   % burn-in draws
nuub   = 100;   % upper bound for nu

%% ---- 2. Prior hyperparameters ------------------------------------------------
% Minnesota prior on VAR coefficients A:
c1 = 0.2^2;        % shrinkage for slope coefficients
c2 = 100;           % prior variance for intercept

% SV process priors:
nuh0 = 5;                    % IG shape for sigh2
Sh0  = 0.01 * (nuh0 - 1);   % IG scale for sigh2  => prior mean = 0.01
rho0 = 0.9;                  % prior mean for rho
Vrho = 0.2^2;                % prior variance for rho

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
        k = n * p + 1;           % regressors (with intercept)

        fprintf('\n=== DGP%d  T=%d  n=%d  k=%d  (R=%d reps) ===\n', ...
                dgp, T_sim, n, k, R);

        for r = 1 : R
            fprintf('  Rep %2d/%d ... ', r, R);
            t_rep = tic;

            rng(1000 * dgp + T_sim * 100 + r);   % reproducible seed

            %% Data setup
            Y0     = Y0_fixed;
            shortY = all_Y{r};
            T      = T_sim;
            tmpY   = [Y0; shortY];

            %% Construct Minnesota prior
            % Fit AR(4) per variable to get residual variance sig2,
            % then build diagonal prior variances VA0.
            A0   = zeros(k, n);
            VA0  = zeros(k, 1);
            sig2 = zeros(n, 1);
            p_ar = min(4, size(tmpY, 1) - 2);
            for i = 1:n
                T_ar = size(tmpY, 1) - p_ar - 1;
                Z_ar = ones(T_ar, 1);
                for lag = 1:p_ar
                    Z_ar = [Z_ar, tmpY(p_ar+1-lag : end-lag-1, i)]; %#ok<AGROW>
                end
                y_ar = tmpY(p_ar+2 : end, i);
                tmpb = (Z_ar' * Z_ar) \ (Z_ar' * y_ar);
                sig2(i) = mean((y_ar - Z_ar * tmpb).^2);
            end
            VA0(1) = c2;
            for i = 2:k
                l   = ceil((i-1) / n);
                idx = mod(i-1, n);
                if idx == 0, idx = n; end
                VA0(i) = c1 / (l^2 * sig2(idx));
            end

            %% Construct regressor matrix X (with intercept)
            X_lags = zeros(T, n * p);
            for lag = 1:p
                X_lags(:, (lag-1)*n+1 : lag*n) = tmpY(p-lag+1 : end-lag, :);
            end
            X = [ones(T, 1), X_lags];

            S0  = eye(n);
            nu0 = n + 3;

            %% Storage arrays
            store_Sig   = zeros(n, n);
            store_A     = zeros(k, n);
            store_h     = zeros(nsims, T);
            store_lam   = zeros(nsims, T);
            store_theta = zeros(nsims, 3);    % [nu, rho, sigh2]

            counth = 0; countrho = 0; countnu = 0;

            %% Chain initialisation
            h     = zeros(T, 1);
            nu    = 5;
            rho   = 0.8;
            sigh2 = 0.1;
            lam   = 1 ./ gamrnd(nu/2, 2/nu, T, 1);

            %% MCMC loop
            for isim = 1 : nsims + burnin

                %% Step 1: Sample Sigma and A
                iOm  = spdiags(exp(-h) ./ lam, 0, T, T);
                XiOm = X' * iOm;
                KA   = spdiags(1./VA0, 0, k, k) + XiOm * X;
                Ahat = KA \ (spdiags(1./VA0, 0, k, k) * A0 + XiOm * shortY);
                Shat = S0 + A0' * spdiags(1./VA0, 0, k, k) * A0 ...
                         + shortY' * iOm * shortY - Ahat' * KA * Ahat;
                Shat = (Shat + Shat') / 2;
                Sig  = iwishrnd(Shat, nu0 + T);
                CSig = chol(Sig, 'lower');
                A    = Ahat + (chol(KA, 'lower')' \ randn(k, n)) * CSig';

                %% Step 2: Sample h_t
                U   = shortY - X * A;
                tmp = U / CSig';
                s2  = sum(tmp.^2, 2) ./ lam;
                [h, flag] = sample_h(s2, rho, sigh2, h, n);
                counth = counth + flag;

                %% Step 3: Sample lambda_t
                s2_lam = sum(tmp.^2, 2) ./ exp(h);
                lam    = 1 ./ gamrnd((n + nu)/2, 2 ./ (s2_lam + nu));

                %% Step 4: Sample nu
                [nu, flag, ~] = sample_nu(lam, nu, nuub);
                countnu = countnu + flag;

                %% Step 5: Sample sigh2
                eh    = [h(1) * sqrt(1 - rho^2); h(2:end) - rho * h(1:end-1)];
                sigh2 = 1 / gamrnd(nuh0 + T/2, 1 / (Sh0 + sum(eh.^2)/2));

                %% Step 6: Sample rho via MH
                Krho   = 1/Vrho + sum(h(1:T-1).^2) / sigh2;
                rhohat = Krho \ (rho0/Vrho + h(1:T-1)' * h(2:T) / sigh2);
                rhoc   = rhohat + sqrt(Krho)' \ randn;
                grho   = @(x) -0.5*log(sigh2./(1-x.^2)) ...
                               - 0.5*(1-x.^2)/sigh2 * h(1)^2;
                if abs(rhoc) < 0.9999
                    alpMH = exp(grho(rhoc) - grho(rho));
                    if alpMH > rand
                        rho = rhoc;
                        countrho = countrho + 1;
                    end
                end

                %% Store post-burnin
                if isim > burnin
                    isave = isim - burnin;
                    store_A   = store_A + A;
                    store_Sig = store_Sig + Sig;
                    store_h(isave, :)     = h';
                    store_lam(isave, :)   = lam';
                    store_theta(isave, :) = [nu, rho, sigh2];
                end

            end  % MCMC loop

            %% Posterior summaries
            lam_hat   = mean(store_lam, 1)';
            h_hat_raw = mean(store_h, 1)';
            h_hat     = h_hat_raw - mean(h_hat_raw);

            nu_hat    = mean(store_theta(:, 1));
            rho_hat   = mean(store_theta(:, 2));
            sigh2_hat = mean(store_theta(:, 3));

            %% True values (mean-centre h)
            lam_true   = all_lam{r};
            h_true_raw = all_h{r};
            h_true     = h_true_raw - mean(h_true_raw);

            %% Save
            out_file = sprintf('%sSVt_DGP%d_T%d_rep%02d.mat', ...
                                results_dir, dgp, T_sim, r);
            save(out_file, ...
                'lam_hat', 'h_hat', ...
                'lam_true', 'h_true', ...
                'nu_hat', 'rho_hat', 'sigh2_hat', ...
                'dgp', 'T_sim', 'r');

            fprintf('done (%.1f sec) -> %s\n', toc(t_rep), out_file);

        end  % rep loop
    end  % T_sim loop
end  % dgp loop

fprintf('\n=== BVAR-t-CSV simulation complete. Results in %s ===\n', results_dir);
