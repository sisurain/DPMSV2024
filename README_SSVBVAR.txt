=========================================================================
README — SSVBVAR replication package
=========================================================================

This package contains the replication code for:

  "Semiparametric Stochastic Volatility for Large Bayesian VARs"
  Frank C. Z. Wu, Joshua C. C. Chan, and Yong Song

The paper develops Bayesian VARs whose innovations follow a Dirichlet
process mixture over a location vector and a transitory variance scale
(mu_t, lambda_t), combined either with a common scalar stochastic
volatility factor (DPM-CSV) or with the order-invariant multivariate
stochastic volatility structure of Chan, Koop and Yu (2024)
(DPM-OISV).

=========================================================================
Folder structure
=========================================================================

  SSVBVAR/
    README.txt                 this file
    THIRD_PARTY_CODE.txt       provenance of included third-party code
    BVARDPM_sim/               Monte Carlo simulation study
    BVARDPM_fullsample/        full-sample DPM-OISV estimation
    BVARDPM_forecast/          42-variable forecasting exercise
    BVARDPM_connectedness/     bank-network connectedness application

Each folder is self-contained: it carries its own README, data, and
utility files, and can be run independently of the others. See the
folder READMEs for detailed workflows and per-file descriptions.

=========================================================================
Mapping from paper outputs to code
=========================================================================

  Tables 1-3    BVARDPM_sim/compute_MSE.m
                (after running the gen_data pipeline and the three
                model runners; see BVARDPM_sim/README.txt)

  Figure 1      BVARDPM_fullsample/post_analysis.m
  Figure 2      BVARDPM_fullsample/DPM_OIMSV/
                post_tail_balance_DPM_OIMSV.m

  Table 4,      BVARDPM_forecast/postprocess_forecast_results.m
  Table 5,
  Tables 7-11
  Figures 3-4   BVARDPM_forecast/make_spaghetti_plots.m
  Figure 5      BVARDPM_forecast/postprocess_forecast_results.m

  Figure 6      BVARDPM_connectedness/make_connectedness_plots_DY.m
  Table 6       BVARDPM_connectedness/extract_DY_summary.m

=========================================================================
Data
=========================================================================

  FRED-MD (January 2026 vintage), Federal Reserve Bank of St. Louis:
    included as csv in BVARDPM_sim/gen_data/, BVARDPM_forecast/data/,
    and BVARDPM_fullsample/.

  Daily bank realized-volatility panel from Demirer, Diebold, Liu and
  Yilmaz (2018, Journal of Applied Econometrics): included with
  attribution in BVARDPM_connectedness/data/; see that folder's README.

=========================================================================
Software
=========================================================================

  MATLAB R2020a or later with the Statistics and Machine Learning
  Toolbox. The forecasting exercise uses the Parallel Computing
  Toolbox (parfor); without it the loops run serially. Some CCMM
  utility files in the simulation folder additionally use the Control
  System Toolbox (ltitr).

=========================================================================
Results files
=========================================================================

Large intermediate results (simulated datasets, per-replication and
per-origin estimation output, full-sample posterior draws) are not
included in this package because of their size. Every table and
figure in the paper is reproducible from the raw data by running the
workflows described in the folder READMEs. All scripts use fixed
random seeds where applicable.

=========================================================================
