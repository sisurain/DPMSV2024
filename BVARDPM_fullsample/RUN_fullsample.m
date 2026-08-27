% =========================================================================
% RUN_fullsample.m
%
% Master script to run the full-sample estimation of the DPM-OIMSV
% model on the CCMM 16-variable dataset.
%
% Steps:
%   1. Run STEP0 to prepare the data (if CCMM16_data.mat does not exist)
%   2. Copy the data file into the model subfolder
%   3. Run the DPM-OIMSV estimation
%   4. Run the post-estimation analysis
%
% Folder structure:
%   BVARDPM_fullsample/
%   |-- STEP0_prepare_CCMM_data.m   (data preparation)
%   |-- 202601MD.csv                (raw FRED-MD data)
%   |-- CCMM16_data.mat             (prepared data, output of STEP0)
%   |-- RUN_fullsample.m            (this file)
%   |-- post_analysis.m             (post-estimation analysis)
%   +-- DPM_OIMSV/
%       |-- BVAR_DPM_OIMSV_FE.m     (main estimation)
%       |-- get_resid_var_OISV_FE.m, get_C_OISV_FE.m
%       |-- nigMSV.m, oneDPM_MSV.m, randmn.m
%       |-- sample_SV_OISV.m, sample_SV0para_OISV.m
%       |-- construct_Sigt_DPM_OISV.m
%       +-- results_DPM_OIMSV.mat   (output)
%
% =========================================================================

clear; clc; close all;

fprintf('=========================================================\n');
fprintf('  Full-Sample Estimation: DPM-OIMSV\n');
fprintf('=========================================================\n\n');

%% Step 0: Prepare Data
if ~exist('CCMM16_data.mat', 'file')
    fprintf('CCMM16_data.mat not found. Running STEP0_prepare_CCMM_data.m...\n');
    if exist('STEP0_prepare_CCMM_data.m', 'file')
        run('STEP0_prepare_CCMM_data.m');
    else
        error('STEP0_prepare_CCMM_data.m not found in current directory.');
    end
else
    fprintf('CCMM16_data.mat found. Skipping data preparation.\n');
    load('CCMM16_data.mat', 'n', 'T', 'p');
    fprintf('  n=%d variables, T=%d observations, p=%d lags\n', n, T, p);
end

%% Step 1: Copy data to the model subfolder
copyfile('CCMM16_data.mat', 'DPM_OIMSV/CCMM16_data.mat');
fprintf('Data file copied to the model subfolder.\n\n');

%% Step 2: Run DPM-OIMSV
fprintf('=========================================================\n');
fprintf('  Running DPM-OIMSV estimation...\n');
fprintf('=========================================================\n');
cd('DPM_OIMSV');
try
    run('BVAR_DPM_OIMSV_FE.m');
catch ME
    fprintf('ERROR in DPM-OIMSV: %s\n', ME.message);
end
cd('..');

%% Step 3: Post-analysis
fprintf('\n=========================================================\n');
fprintf('  Running post-estimation analysis...\n');
fprintf('=========================================================\n');
addpath('DPM_OIMSV');
run('post_analysis.m');

fprintf('\n=========================================================\n');
fprintf('  All done.\n');
fprintf('=========================================================\n');
