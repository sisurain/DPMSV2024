=========================================================================
README — model_DPMCSV/  (BVARDPM_forecast)
=========================================================================

This folder runs the expanding-window forecasting exercise for the
DPM-CSV model: the Dirichlet process mixture innovation specification
combined with a common scalar stochastic volatility factor.

Model specification:

  y_t = X_t * A + mu_{z_t} + u_t
  u_t ~ N(0, lambda_{z_t} * exp(h_t) * Sigma)
  h_t = rho * h_{t-1} + eta_t,  eta_t ~ N(0, sigh2)   [scalar AR(1) CSV]
  Sigma ~ IW(S0, nu0)
  (mu_k, lambda_k) ~ DPM with normal-inverse-gamma base measure
  A ~ Minnesota prior (no intercept; the DPM cluster means mu_k
  absorb the location)

=========================================================================
Execution
=========================================================================

  run_DPMCSV_forecast.m

    Input:   ../data/FREDMD42_forecast_data.mat
    Output:  ../results/results_DPMCSV_forecast.mat

  Run from within the model_DPMCSV/ directory. The script loops over
  all forecast origins in parallel (parfor). Each origin re-estimates
  the model on the available data (10,000 posterior draws after 5,000
  burn-in) and simulates one predictive path per kept draw. Per-origin
  random seeds make origins independent of the parallel configuration.

=========================================================================
File list
=========================================================================

  run_DPMCSV_forecast.m
    Driver script: loads the data, allocates storage, runs the
    per-origin estimation and forecasting in parallel, and saves the
    stacked results.

  fcst_DPMCSV_oneorigin.m
    Estimation and forecasting at a single origin. One MCMC iteration
    consists of: a joint draw of the VAR coefficients A and covariance
    Sigma from their natural-conjugate conditional (matric-normal /
    inverse-Wishart); the common log-volatility h via an
    acceptance-rejection Metropolis-Hastings step; the SV parameters
    rho and sigh2; DPM cluster assignments via collapsed Gibbs;
    cluster parameters (mu_k, lambda_k) from their normal-inverse-
    gamma posteriors; and the DPM concentration parameter. After
    burn-in, each iteration simulates a forecast path and evaluates
    the one-step-ahead joint predictive density.

  nigSV.m
    Handle class implementing the normal-inverse-gamma sufficient
    statistics for the DPM base measure under the common-volatility
    covariance exp(h_t)*Sigma. Supports sequential addSample/delSample
    updates and computes the multivariate Student-t log predictive
    density.

  oneDPM.m
    One-pass collapsed Gibbs sampler used to initialise the DPM
    cluster assignments.

  randmn.m
    Categorical distribution sampler (inverse transform).

  sample_h.m
    Samples the common scalar stochastic volatility h using an
    acceptance-rejection Metropolis-Hastings algorithm with a
    Gaussian proposal constructed via Newton-Raphson.

  construct_prior_A_func.m
    Constructs the Minnesota prior for the VAR slope coefficients
    (no-intercept version): prior means, prior variances, and
    AR-residual scaling factors.

  crpsDraws.m
    CRPS from predictive draws (empirical estimator of Krueger,
    Lerch, Thorarinsdottir and Gneiting, 2021). From the Carriero,
    Clark, Marcellino and Mertens (2024) replication codebase.

=========================================================================
Output format
=========================================================================

results_DPMCSV_forecast.mat contains, with N variables, H forecast
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
