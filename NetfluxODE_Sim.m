function [stimIdx, geneIdx, validatedIdx, ylong_4h, y4h_norm] = NetfluxODE_Sim(w_opt)
%% Loading Parameters
% Input Stimuli
S_in = [0.400, 0.7]; % Stretch Inputs (Baseline, Stimulated), change as necessary

% Load parameters 
[params, y0] = NetfluxODE_CaoS_loadParams();
[rpar, ~, ~, speciesNames, speciesTypes] = params{:};

% Identify Stimulus and Gene indices
stimIdx = find(strcmp(speciesNames, 'Stretch')); 
geneIdx = find(contains(speciesTypes, 'Gene')); 
% validatedIdx = 85:654;
% optIdx = 12:20;
[validatedIdx, ~, ~] = ParseGeneTFInputs_CaoS(geneIdx);
% optIdx = [47, 71, 72, 73, 79];
% optIdx = [24, 44, 45, 46, 47, 48, 55, 56, 57, 64, 65, 66, 67, 68, 71, 72, 73, 79, 80];
% optIdx = [47, 71, 72, 73, 79];
% optIdx = [47, 73];
optIdx = [12:20, 47, 73, 79];

% Defining Time-Scale
tspan = 0:0.1:(24*60); % used for experimental simulation timecourse
tspan1 = [0 30000]; % used for baseline simulation, should be sufficiently long to allow for 'steady-state' baseline

%% RNA Data Import
[y0_MA, ~, ~, y4h_MA, ~, ~] = ImportRNA_CaoS(speciesNames, speciesTypes, y0);
y4h_norm = y4h_MA ./ y0_MA;

%% Baseline Simulation
rpar(1,1) = 1; % Set w1 to intial stimulus as necessary
rpar(1, optIdx) = w_opt; %inputting optimized/guessed weights into the main vector to run the simulations with
params{1} = rpar;
[~, y1] = ode15s(@NetfluxODE_CaoS, tspan1, y0, [], params, S_in(1));
ss = y1(end, :)'; % End ouputs from SS model

%% Longitudinal Stretch Simulation
% Activated Simulation at S = 0.7
y0_new3 = ss; % sets initial conditions for activated simulation as the outputs from the baseline
y0_new3(stimIdx) = S_in(1); % sets stimulus output to baseline value
[t2, y2] = ode15s(@NetfluxODE_CaoS, tspan, y0_new3, [], params, S_in(2));


% Normalize longitudinal outputs by steady state
[ylong_4h, ~, ~] = ODENorm(t2, y2, ss, geneIdx);
end