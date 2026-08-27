=========================================================================
README — model_DPMOISV/  (BVARDPM_forecast)
=========================================================================

This folder runs the expanding-window forecasting exercise for the
DPM-OISV model: the Dirichlet process mixture innovation specification
combined with the order-invariant stochastic volatility structure of
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

  run_DPMOISV_forecast.m

    Input:   ../data/FREDMD42_forecast_data.mat
    Output:  ../results/results_DPMOISV_forecast.mat

  Run from within the model_DPMOISV/ directory. The script loops over
  all forecast origins in parallel (parfor). Each origin re-estimates
  the model on the available data (10,000 posterior draws after 5,000
  burn-in) and simulates one predictive path per kept draw. Per-origin
  random seeds make origins independent of the parallel configuration.

=========================================================================
File list
=========================================================================

  run_DPMOISV_forecast.m
    Driver script: loads the data, allocates storage, runs the
    per-origin estimation and forecasting in parallel, and saves the
    stacked results.

  fcst_DPMOISV_oneorigin.m
    Estimation and forecasting at a single origin. One MCMC iteration
    consists of: horseshoe shrinkage hyperparameters; VAR coefficients
    A equation by equation; the contemporaneous matrix B0 row by row
    (order-invariant step); individual log-volatilities via the KSC
    mixture sampler; SV parameters (phi_i, omega2_i); DPM cluster
    assignments via collapsed Gibbs; cluster parameters
    (mu_k, lambda_k) from their normal-inverse-gamma posteriors; and
    the DPM concentration parameter. After burn-in, each iteration
    simulates a forecast path and evaluates the one-step-ahead joint
    predictive density.

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
    log-volatilities with a stationary AR(1) prior, implemented with
    a sparse precision sampler. Adapted from the Chan, Koop and Yu
    (2024) replication code.

  sample_SV0para_OISV.m
    Samples the SV persistence parameters phi_i (Metropolis-Hastings
    with a Gaussian proposal and stationary-initial-condition
    correction) and innovation variances omega2_i (Gibbs). Adapted
    from the Chan, Koop and Yu (2024) replication code.

  get_C_OISV_FE.m
    Scaling matrix for the Minnesota-type horseshoe prior on the VAR
    coefficients (no-intercept version; k = n*p regressors).

  get_resid_var_OISV_FE.m
    Univariate AR(4) residual variances (no intercept) used to scale
    the Minnesota prior.

  crpsDraws.m
    CRPS from predictive draws (empirical estimator of Krueger,
    Lerch, Thorarinsdottir and Gneiting, 2021). From the Carriero,
    Clark, Marcellino and Mertens (2024) replication codebase.

=========================================================================
Output format
=========================================================================

results_DPMOISV_forecast.mat contains, with N variables, H forecast
horizons, and J forecast origins:

  fcstYrealized      N x H x J     realized values
  fcstYhat           N x H x J     predictive means
  fcstYmedian        N x H x J     predictive medians
  fcstCRPS           N x H x J     CRPS per variable and horizon
  fcstLogscore       N x H x J     marginal log predictive densities
  fcstJointLogscore  1 x J         one-step-ahead joint log scores
  fcstYquantiles     N x H x 5 x J predictive quantiles [5,16,50,84,95]
  fcstPIT            N x H x J     probability integral transforms
  fcstYdraws_h1      N x 1000 x J  subsampled h=1 predictive draws
  Tjumpoffs, Njumpoffs, dates_trans, labels, codes, and the MCMC
  settings used.

=========================================================================
