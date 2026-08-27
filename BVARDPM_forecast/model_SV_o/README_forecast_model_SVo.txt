=========================================================================
README — model_SVo/  (BVARDPM_forecast)
=========================================================================

This folder runs the expanding-window forecasting exercise for the SVo
benchmark: the BVAR with stochastic volatility and a common discrete
outlier scale of Carriero, Clark, Marcellino and Mertens (2024).

Model specification:

  y_t = PAI(L) * y_{t-1} + c + v_t
  v_t = obar_t * A^{-1} * Lambda_t^{1/2} * e_t,   e_t ~ N(0, I)
  obar_t in {1, 2, ..., 20}  (scalar common outlier scale)
  P(obar_t = 1) = 1 - pi,  P(obar_t = k > 1) = pi / 19
  pi ~ Beta(2.5, 117.5)  (prior outlier frequency of one per four years)
  log volatilities follow correlated random walks (KSC mixture sampler)
  PAI ~ Minnesota prior (CCMM settings)

=========================================================================
Execution
=========================================================================

  run_SVo_forecast.m

    Input:   ../data/FREDMD42_forecast_data.mat
    Output:  ../results/results_SV_o_forecast.mat

  Run from within the model_SVo/ directory. The script loops over all
  forecast origins in parallel (parfor). Each origin re-estimates the
  model on the available data (1,000 posterior draws after 200 burn-in)
  and simulates 10 predictive paths per kept draw, for 10,000
  predictive draws per origin. Per-origin random seeds make origins
  independent of the parallel configuration.

=========================================================================
File list
=========================================================================

  run_SVo_forecast.m
    Driver script: loads the data, allocates storage, runs the
    per-origin estimation and forecasting in parallel, and saves the
    stacked results.

  fcst_SVo_oneorigin.m
    Estimation and forecasting at a single origin, following the Gibbs
    ordering of the CCMM sampler (mcmcVARSVobar.m): VAR coefficients
    via the triangular algorithm; impact matrix A equation by
    equation; the common outlier scale obar_t by Gibbs on the discrete
    grid; log volatilities via the KSC mixture sampler; and the
    outlier probability from its Beta posterior. After burn-in, each
    iteration simulates forecast paths and evaluates the one-step-
    ahead joint predictive density.

--- Utility files from the CCMM replication codebase ---

The following files are from the replication code of Carriero, Clark,
Marcellino and Mertens and are included unmodified.

  CTA.m
    Triangular algorithm for drawing the VAR coefficients conditional
    on the impact matrix and stochastic volatilities.

  StochVolKSCcorrsqrt.m
    Kim-Shephard-Chib (1998) mixture sampler for correlated random-
    walk stochastic volatilities, using the precision-based approach.

  getKSC7values.m
    Seven-component normal mixture approximation constants from Kim,
    Shephard and Chib (1998).

  rwnoisePrecisionBasedSampler.m
    Precision-based sampler for the random-walk-plus-noise state
    space representation (used by StochVolKSCcorrsqrt).

  vec.m
    Vectorises a matrix: vec(X) = X(:) (used by CTA).

  crpsDraws.m
    CRPS from predictive draws (empirical estimator of Krueger,
    Lerch, Thorarinsdottir and Gneiting, 2021).

=========================================================================
Output format
=========================================================================

results_SV_o_forecast.mat contains, with N variables, H forecast
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
