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
%% Defining Optimized Indices

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

%% Define optimization variables: n/EC50
% % Parameters to optimize (same optIdx)
% n_params = length(optIdx);
% 
% % Define EC50 and n bounds
% n_mu = 1.4;      % Example clustering around 1.4
% n_sigma = 0.2;
% n_lower = 0.5;
% n_upper = 4.0;
% 
% EC50_lower = 0.01;  % Avoiding zero
% EC50_upper = 0.7;   % Upper limit, will be adjusted dynamically
% 
% % Preallocate guess matrices
% EC50_guesses = zeros(n_params, n_params);
% n_guesses = zeros(n_params, n_params);
% 
% % Generate guesses with constraints
% for i = 1:n_params
%     for j = 1:n_params
%         % Sample n within bounds
%         n_val = n_mu + n_sigma * randn;
%         while n_val < n_lower || n_val > n_upper
%             n_val = n_mu + n_sigma * randn;
%         end
% 
%         % Max allowable EC50 for this n_val
%         EC50_max = (0.5)^(1 / n_val);  
%         EC50_max = min(EC50_max, EC50_upper);
% 
%         % Sample EC50 within bounds and respecting constraint
%         ec_val = EC50_max * rand();  % Uniform within 0 to EC50_max
%         if ec_val < EC50_lower
%             ec_val = EC50_lower;
%         end
% 
%         % Store values
%         n_guesses(i, j) = round(n_val, 4);
%         EC50_guesses(i, j) = round(ec_val, 4);
%     end
% end
% 
% % Stack guesses into a single parameter.guess matrix [EC50; n]
% parameters.name = [arrayfun(@(i) sprintf('EC50_%d', i), 1:n_params, 'UniformOutput', false), ...
%                    arrayfun(@(i) sprintf('n_%d', i), 1:n_params, 'UniformOutput', false)];
% parameters.number = 2 * n_params;
% parameters.min = [EC50_lower * ones(n_params, 1); n_lower * ones(n_params, 1)];
% parameters.max = [0.5 * ones(n_params, 1); n_upper * ones(n_params, 1)];
% 
% % Flatten the guesses
% parameters.guess = [EC50_guesses; n_guesses];
% 
% disp('Initial EC50 guesses:');
% disp(EC50_guesses);
% 
% disp('Initial n guesses:');
% disp(n_guesses);

%% Options
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
xlabel('Multi-Start Index'); ylabel('Final Minimization');
title('Error per PESTO Start (10 Runs): One TF Rxn');