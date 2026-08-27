=========================================================================
README — BVARDPM_sim/
=========================================================================

This folder contains the Monte Carlo simulation study comparing three
Bayesian VAR volatility specifications in their ability to recover the
transitory variance scale lambda_t, the persistent log-volatility h_t,
and the combined volatility scale tau_t = sqrt(lambda_t)*exp(h_t/2):

  SVo      BVAR-SV with common discrete outlier scaling
           (Carriero, Clark, Marcellino & Mertens, 2022)
  SV-t     BVAR-t-CSV with Student-t fat-tail scaling
           (Chan, 2020)
  DPM-CSV  BVAR with Dirichlet Process Mixture and common stochastic
           volatility (this paper)

Simulated data are generated from a VAR calibrated to the CCMM
16-variable monthly FRED-MD dataset, with three data-generating
processes for lambda_t (fat-tail continuous mixture, discrete
outliers, clustered pandemic block), two sample sizes (T = 300, 800),
and R = 20 replications per cell. The DGP specifications are described
in detail in gen_data/README.txt.

=========================================================================
Folder structure
=========================================================================

  BVARDPM_sim/
    README.txt              this file
    MASTER_run_all.m        runs all three model runners + compute_MSE
    compute_MSE.m           computes the final MSE comparison tables
    gen_data/               data preparation and DGP generation
    sim_data/               simulated datasets (populated at run time)
    model_SVo/              SVo estimation on the simulated data
    model_SVt/              BVAR-t-CSV estimation on the simulated data
    model_DPMCSV/           DPM-CSV estimation on the simulated data

Each subfolder contains its own README describing its files in detail.

=========================================================================
Execution order (full reproduction)
=========================================================================

Step 1 — Generate the simulation datasets (gen_data/):

  1a.  STEP0_prepare_CCMM_data.m
  1b.  STEP1_run_SVo_fullsample.m
  1c.  STEP2_generate_sim_datasets.m

  This produces the six dataset files

    simdata_DGP{1,2,3}_T{300,800}_R20.mat

  which must reside in sim_data/ for the model runners to find them.

Step 2 — Estimate all three models on every dataset:

  Either run MASTER_run_all.m from this folder (runs the three model
  runners sequentially and then computes the MSE tables), or run each
  model independently from within its own subfolder:

    model_SVo/run_SVo_sim.m
    model_SVt/run_SVt_sim.m
    model_DPMCSV/run_DPMCSV_sim.m

  Each runner writes 120 result files (3 DGPs x 2 sample sizes x 20
  replications) to its own results/ subfolder. The three runners are
  independent of each other and can be executed in separate MATLAB
  sessions to parallelise across models. The DPM-CSV runner is the
  most computationally intensive of the three.

Step 3 — Compute the MSE comparison tables:

  compute_MSE.m (run from this folder). Loads all result files,
  applies the identification correction (lambda normalised to the
  non-outlier baseline, h centred on the non-outlier regime mean,
  tau normalised by its non-outlier mean), and prints:

    Tables 1-3:  MSE (with standard errors) for lambda, h, tau
    Tables 4-6:  Relative MSE (normalised by DPM-CSV)
    Time-profile MSE split by non-outlier vs outlier regime

  All quantities are also saved to MSE_results.mat.

All scripts use fixed random seeds, so results are exactly
reproducible.

=========================================================================
Files not included
=========================================================================

The simulated datasets (contents of sim_data/) and the per-replication
estimation output (contents of model_*/results/) are not included in
this package because of their size. Both are fully reproducible by
running the steps above.

=========================================================================
