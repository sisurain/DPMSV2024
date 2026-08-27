=========================================================================
README — model_SVt/  (BVARDPM_forecast)
=========================================================================

This folder runs the expanding-window forecasting exercise for the
SV-t benchmark: the BVAR-t-CSV model of Chan (2020), combining a
common scalar stochastic volatility factor with i.i.d. Student-t
fat-tail scaling.

Model specification:

  y_t = c + X_t * A + u_t
  u_t ~ N(0, lambda_t * exp(h_t) * Sigma)
  lambda_t ~ IG(nu/2, nu/2)         i.i.d. fat-tail scaling
  h_t = rho * h_{t-1} + eta_t,  eta_t ~ N(0, sigh2)
  Sigma ~ IW(S0, nu0)
  A ~ Minnesota prior (with intercept)
  nu ~ Metropolis-Hastings with Newton-Raphson proposal

=========================================================================
Execution
=========================================================================

  run_SVt_forecast.m

    Input:   ../data/FREDMD42_forecast_data.mat
    Output:  ../results/results_SVt_forecast.mat

  Run from within the model_SVt/ directory. The script loops over all
  forecast origins in parallel (parfor). Each origin re-estimates the
  model on the available data (10,000 posterior draws after 5,000
  burn-in) and simulates one predictive path per kept draw. Per-origin
  random seeds make origins independent of the parallel configuration.

=========================================================================
File list
=========================================================================

  run_SVt_forecast.m
    Driver script: loads the data, allocates storage, runs the
    per-origin estimation and forecasting in parallel, and saves the
    stacked results.

  fcst_SVt_oneorigin.m
    Estimation and forecasting at a single origin. One MCMC iteration
    consists of: VAR coefficients A and covariance Sigma jointly;
    the common scalar log-volatility h_t via acceptance-rejection
    Metropolis-Hastings; lambda_t from its inverse-gamma full
    conditional; the degrees of freedom nu via Metropolis-Hastings;
    and the SV parameters (rho, sigh2). After burn-in, each iteration
    simulates a forecast path and evaluates the one-step-ahead joint
    predictive density.

--- Utility files from the Chan (2020) codebase ---

  sample_h.m
    Samples the common scalar stochastic volatility h_t using an
    acceptance-rejection Metropolis-Hastings algorithm with a
    Gaussian proposal constructed via Newton-Raphson.

  sample_nu.m
    Samples the degrees-of-freedom parameter nu using a Metropolis-
    Hastings step with a Gaussian proposal centred at the Newton-
    Raphson mode of the full conditional.

--- Other utilities ---

  crpsDraws.m
    CRPS from predictive draws (empirical estimator of Krueger,
    Lerch, Thorarinsdottir and Gneiting, 2021), following the
    implementation in the Carriero, Clark, Marcellino and Mertens
    replication codebase.

=========================================================================
Output format
=========================================================================

results_SVt_forecast.mat contains, with N variables, H forecast
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
