% Netflux_opt_PESTO.m
% Optimization of w using PESTO and RNA-seq data

clear; clc; close all;

% === Add PESTO to path ===
addpath(genpath('C:\Users\disma\Documents\2024-2025 School Year\CMRG and WTHPA\Coding\Netflux - CaoS\PESTO Optimization\PESTO-master'));


%% Load model parameters
[paramCells, y0] = NetfluxODE_CaoS_loadParams();
[rpar, tau, ymax, speciesNames, speciesTypes] = paramCells{:};

% Define gene indices
geneIdx = find(contains(speciesTypes, 'Gene'));
%% 

% Define validatedIdx (example: overwrite as needed)
validatedIdx = 85:654;
% optIdx = 12:20;
% optIdx = [24, 44, 45, 46, 47, 48, 55, 56, 57, 64, 65, 66, 67, 68, 71, 72, 73, 79, 80];
% optIdx = [47, 71, 72, 73, 79];
% optIdx = [47, 73];
optIdx = [12:20, 47, 73, 79];

% Make optIdx globally accessible for the simulation
setappdata(0, 'optIdx', optIdx);


%% Define optimization variables: w
w_init = 0.9 * ones(length(optIdx), 1);  % Initial guess for 9 parameters

mu = 0.9;           % Cluster center
sigma = 0.05;         % Small std dev to keep values near mu
lower = 0.7;
upper = 1.5;
n_rows = length(optIdx);
n_cols = length(optIdx);

% Preallocate matrix
w_guesses = zeros(n_rows, n_cols);

% Generate values
for i = 1:n_rows
    for j = 1:n_cols
        val = mu + sigma * randn;
        while val < lower || val > upper
            val = mu + sigma * randn;  % Reject and resample if out of bounds
        end
        w_guesses(i, j) = val;
    end
end

% Optional: round slightly for cleaner display
w_guesses = round(w_guesses, 4);

% === Define parameters ===
parameters.name   = arrayfun(@(i) sprintf('w%d', i), 1:length(optIdx), 'UniformOutput', false);
parameters.number = length(optIdx);
parameters.min    = zeros(length(optIdx), 1);
parameters.max    = 1 * ones(length(optIdx), 1);
% parameters.guess  = w_init;
parameters.guess = w_guesses;


disp('Initial guess:');
disp(parameters.guess');

% === Define options ===
options = PestoOptions();
options.obj_type        = 'negative log-posterior';
options.n_starts        = 10;
options.comp_type       = 'sequential';
options.localOptimizer  = 'fmincon';
options.objOutNumber    = 2;  % Expect [ll, grad] outputs
% options.n_starts = 10;

% === Wrap objective ===
llfun = @(w) logLikelihoodFun(w);

%% Run optimization
parameters = getMultiStarts(parameters, llfun, options);

% === Extract best result ===
if isfield(parameters.MS, 'logPost')
    fvals = parameters.MS.logPost;
elseif isfield(parameters.MS, 'fval')
    fvals = [parameters.MS.fval];
else
    error('Cannot extract log-posterior values from optimization results.');
end

validIdx = ~isnan(fvals);
if ~any(validIdx)
    error('All optimization runs failed.');
end

[~, bestIdx] = min(fvals);
disp(['Best log-posterior: ', num2str(fvals(bestIdx))]);
w_best = parameters.MS.par(:, bestIdx);

disp('--- Optimized w ---');
disp(w_best');

figure;
rmse_vals = parameters.MS.logPost;  % logPost = RMSE
plot(1:length(rmse_vals), rmse_vals, 'o-');
xlabel('Multi-Start Index'); ylabel('Final RMSE');
title('Error  per PESTO Start (10 Runs): One TF Rxn');