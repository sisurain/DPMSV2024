=========================================================================
README — gen_data/
=========================================================================

This folder generates the simulated datasets used in the Monte Carlo
study (Section X of the paper). The procedure has three steps:

  STEP 0:  Prepare the CCMM 16-variable dataset from the raw FRED-MD
           vintage.
  STEP 1:  Estimate the SVo model (Carriero, Clark, Marcellino & Mertens,
           2022) on the full sample to obtain fixed DGP parameters
           (VAR coefficients, base covariance, SV persistence).
  STEP 2:  Generate R=20 Monte Carlo replications for each combination
           of DGP {1,2,3} and sample size T in {300, 800}.

The simulated datasets are saved to the sim_data/ folder. They are not
included in this replication package because of their size, but can be
fully reproduced by running the three steps in order.

=========================================================================
Execution order
=========================================================================

  1.  STEP0_prepare_CCMM_data.m
        Input:   2026-01-MD.csv
        Output:  CCMM16_data.mat

  2.  STEP1_run_SVo_fullsample.m
        Input:   CCMM16_data.mat
        Output:  SVo_fullsample_params.mat

  3.  STEP2_generate_sim_datasets.m
        Input:   SVo_fullsample_params.mat, CCMM16_data.mat
        Output:  sim_data/simdata_DGP{d}_T{T}_R20.mat  (6 files)

=========================================================================
DGP specifications
=========================================================================

All three DGPs share the same fixed VAR dynamics (A_hat, Sigma_hat) and
the same AR(1) common stochastic volatility process (phi_hat, sigh2_hat),
calibrated from the SVo full-sample posterior means. Only the transitory
scale lambda_t differs across DGPs:

  DGP 1 — Fat-tail continuous mixture (designed to favour SV-t):
    lambda_t ~ (1-q)*IG(nu1/2, nu1/2) + q*kappa*IG(nu2/2, nu2/2)
    q = 0.2, nu1 = 10, nu2 = 3, kappa = 5

  DGP 2 — Discrete outliers (designed to favour SVo):
    lambda_t = 1 with prob 1-q, lambda_t = o_t^2 with prob q
    q = 0.05, o_t ~ U{2, ..., 20}

  DGP 3 — Clustered pandemic block (designed to favour DPM-CSV):
    lambda_t = 1 for all t except a 6-observation block
    Block at t0 = floor(0.7*T), lambda_star = 25

=========================================================================
File list
=========================================================================

--- Data ---

  2026-01-MD.csv
    Raw FRED-MD vintage (January 2026). Source: Federal Reserve Bank of
    St. Louis. Input to STEP0.

--- Main scripts (run in order) ---

  STEP0_prepare_CCMM_data.m
    Loads the raw FRED-MD CSV, extracts the 16 CCMM variables, applies
    transformations (annualised growth rates, logs, or levels) following
    Table S.1 of CCMM (2022), and saves the prepared dataset.

  STEP1_run_SVo_fullsample.m
    Runs the SVo Gibbs sampler (mcmcVARSVobar.m) on the full CCMM
    16-variable sample to obtain posterior mean estimates of the VAR
    coefficients, impact matrix, SV persistence, and SV innovation
    variance. These serve as the fixed DGP parameters for simulation.

  STEP2_generate_sim_datasets.m
    Generates R=20 replications per (DGP, T) cell. For each replication,
    draws a common log-volatility path h_t from the calibrated AR(1)
    process, draws lambda_t from the DGP-specific distribution, and
    simulates VAR observations using the fixed coefficients and real
    initial conditions from the data.

--- Utility files from CCMM (2022) codebase ---

The following files are from the replication code of:

  Carriero, A., Clark, T.E., Marcellino, M. and Mertens, E. (2022).
  Addressing COVID-19 Outliers in BVARs with Stochastic Volatility.
  Review of Economics and Statistics, 104(3), 517-527.

These files are used by STEP1 and should not be modified.

  mcmcVARSVobar.m
    Main Gibbs sampler for the BVAR-SV model with common outlier scaling.

  CTA.m
    Triangular algorithm for drawing VAR coefficients conditional on
    the Cholesky factor and stochastic volatility (Carriero, Clark &
    Marcellino, 2015).

  obarGibbsdraw.m
    Draws the common outlier scale o_t and outlier probability from
    their conditional posteriors.

  StochVolKSCcorrsqrt.m
    Kim-Shephard-Chib (1998) mixture sampler for stochastic volatility
    with correlated innovations, using the precision-based approach.

  KSC.m
    KSC mixture smoother: draws mixture indicator states and then
    log-volatility states via Carter-Kohn (1994) backward simulation.

  getKSC7values.m
    Returns the 7-component normal mixture approximation constants
    from Kim, Shephard and Chib (1998).

  rwnoisePrecisionBasedSampler.m
    Precision-based sampler for the random-walk-plus-noise state space
    representation of the SV process (Chan and Jeliazkov, 2009).

  betadraw.m
    Draws from a Beta distribution via the ratio-of-chi-squares method.

  checkdiff.m
    Numerical comparison utility for debugging.

  parid.m
    Returns the current parallel worker ID (or 1 if not in a parpool).

  progressbar.m
    Text/graphical progress bar for long MCMC runs.

  vec.m
    Vectorises a matrix: vec(X) = X(:).

=========================================================================
