% =========================================================================
% MASTER_run_all.m
%
% PURPOSE: Run all three simulation estimation scripts and produce
%          the final MSE comparison tables.
%
%   If results already exist in model_X/results/ folders, go directly
%   to compute_MSE.m — no need to run this file.
%
% FOLDER STRUCTURE:
%   BVARDPM_sim/
%     MASTER_run_all.m      <- this file
%     compute_MSE.m
%     sim_data/             <- DGP datasets (from gen_data scripts)
%     model_SVo/
%       run_SVo_sim.m
%       results/
%     model_SVt/
%       run_SVt_sim.m
%       results/
%     model_DPMCSV/
%       run_DPMCSV_sim.m
%       results/
%
% NOTE: Each model subfolder can be run independently in separate
%       MATLAB sessions to parallelise across models.
% =========================================================================

clear; clc;

main_dir = pwd;
fprintf('MASTER_run_all.m: starting from %s\n\n', main_dir);
master_start = tic;

%% ---- 1. SVo -------------------------------------------------------------
fprintf('======================================================\n');
fprintf('Running model_SVo/run_SVo_sim.m\n');
fprintf('======================================================\n');
t0 = tic;
cd(fullfile(main_dir, 'model_SVo'));
run('run_SVo_sim.m');
cd(main_dir);
fprintf('model_SVo done in %.1f minutes.\n\n', toc(t0)/60);

%% ---- 2. SVt -------------------------------------------------------------
fprintf('======================================================\n');
fprintf('Running model_SVt/run_SVt_sim.m\n');
fprintf('======================================================\n');
t0 = tic;
cd(fullfile(main_dir, 'model_SVt'));
run('run_SVt_sim.m');
cd(main_dir);
fprintf('model_SVt done in %.1f minutes.\n\n', toc(t0)/60);

%% ---- 3. DPM-CSV ---------------------------------------------------------
fprintf('======================================================\n');
fprintf('Running model_DPMCSV/run_DPMCSV_sim.m\n');
fprintf('======================================================\n');
t0 = tic;
cd(fullfile(main_dir, 'model_DPMCSV'));
run('run_DPMCSV_sim.m');
cd(main_dir);
fprintf('model_DPMCSV done in %.1f minutes.\n\n', toc(t0)/60);

%% ---- 4. Compute MSE tables ----------------------------------------------
fprintf('======================================================\n');
fprintf('Computing MSE tables via compute_MSE.m\n');
fprintf('======================================================\n');
run('compute_MSE.m');

fprintf('\nMASTER_run_all.m complete. Total time: %.1f minutes.\n', ...
        toc(master_start)/60);
