=========================================================================
README — DPM_OIMSV/  (BVARDPM_fullsample)
=========================================================================

This folder contains the full-sample estimation of the DPM-OIMSV model:
the Dirichlet process mixture innovation specification combined with
the order-invariant multivariate stochastic volatility structure of
Chan, Koop and Yu (2024).

Model specification:

  y_t = X_t * A + mu_{z_t} + u_t
  u_t ~ N(0, lambda_{z_t} * Sigma_t)
  Sigma_t = B0^{-1} * D_t * (B0^{-1})'
  D_t = diag(exp(h_{1,t}), ..., exp(h_{n,t}))
  h_{i,t} = phi_i * h_{i,t-1} + eta_{i,t},  eta_{i,t} ~ N(0, omega2_i)
  (mu_k, lambda_k) ~ DPM with normal-inverse-gamma base measure
  B0 sampled row by row with the order-invariant construction of
  Chan, Koop and Yu (2024)
  A ~ Minnesota-type horseshoe prior (no intercept; the DPM cluster
  means mu_k absorb the location)

=========================================================================
Execution
=========================================================================

  1.  BVAR_DPM_OIMSV_FE.m

        Input:   CCMM16_data.mat  (produced by
                 ../STEP0_prepare_CCMM_data.m and copied here by
                 ../RUN_fullsample.m)
        Output:  results_DPM_OIMSV.mat

        Run from within the DPM_OIMSV/ directory, either directly or
        via ../RUN_fullsample.m. The sampler uses 10,000 posterior
        draws after 5,000 burn-in. One MCMC iteration consists of:
        horseshoe shrinkage hyperparameters; VAR coefficients A
        equation by equation; the contemporaneous matrix B0 row by
        row (order-invariant step); individual log-volatilities via
        the KSC mixture sampler; SV parameters (phi_i, omega2_i);
        DPM cluster assignments via collapsed Gibbs; cluster
        parameters (mu_k, lambda_k) from their normal-inverse-gamma
        posteriors; and the DPM concentration parameter.

  2.  post_tail_balance_DPM_OIMSV.m

        Input:   results_DPM_OIMSV.mat, CCMM16_data.mat
        Output:  tail-asymmetry figure and summary table

        Run from within the DPM_OIMSV/ directory after estimation.
        Computes the tail-balance diagnostic: posterior draws of the
        innovations eps_t = y_t - X_t*A - mu_t; the fixed scales
        sigma_{i,t} = sqrt(E[lambda_t * Sigma_{ii,t}]); standardized
        tail exceedance probabilities at two standard deviations; and
        the asymmetry index A_i = log(p_i^+ / p_i^-) for each
        variable.

=========================================================================
File list
=========================================================================

  BVAR_DPM_OIMSV_FE.m
    Main estimation script (see above).

  post_tail_balance_DPM_OIMSV.m
    Tail-balance diagnostic (see above).

  nigMSV.m
    Handle class implementing the normal-inverse-gamma sufficient
    statistics for the DPM base measure under the time-varying
    covariance Sigma_t = B0^{-1} D_t (B0^{-1})'. Supports sequential
    addSample/delSample updates and computes the multivariate
    Student-t log predictive density.

  oneDPM_MSV.m
    One-pass collapsed Gibbs sampler used to initialise the DPM
    cluster assignments.

  randmn.m
    Categorical distribution sampler (inverse transform).

  sample_SV_OISV.m
    Kim-Shephard-Chib (1998) mixture sampler for the individual
    log-volatilities with a stationary AR(1) prior. Adapted from the
    Chan, Koop and Yu (2024) replication code.

  sample_SV0para_OISV.m
    Samples the SV persistence parameters phi_i and innovation
    variances omega2_i. Adapted from the Chan, Koop and Yu (2024)
    replication code.

  get_C_OISV_FE.m
    Scaling matrix for the Minnesota-type horseshoe prior
    (no-intercept version).

  get_resid_var_OISV_FE.m
    Univariate AR(4) residual variances (no intercept) used to scale
    the Minnesota prior.

=========================================================================
Output format
=========================================================================

results_DPM_OIMSV.mat contains, with T observations, n variables,
k = n*p regressors, and nsims kept draws:

  store_A           k x n x nsims    VAR coefficient draws
  store_B0          n x n x nsims    contemporaneous matrix draws
  store_H           T x n x nsims    log-volatility draws
  store_lam         T x nsims        lambda_t draws
  store_mu          T x n x nsims    DPM location draws
  store_z           T x nsims        cluster assignment draws
  store_DPM_params  nsims x 2        [alpha, K]
  store_SV_params   nsims x 2n       [phi_1..phi_n, omega2_1..omega2_n]

  plus the sample information (dates, labels, dimensions) and all
  prior hyperparameters used in the run.

=========================================================================
