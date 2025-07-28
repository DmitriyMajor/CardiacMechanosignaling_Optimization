function sensitivityResults = SensitivityAnalysis_CaoS(w_best, perturb_frac)
% SensitivityAnalysis_CaoS
% Performs local sensitivity analysis on optimized w parameters
%
% Inputs:
%   w_best        - [n x 1] vector of optimized w parameters (optIdx only)
%   perturb_frac  - Fraction for ± perturbation (e.g., 0.05 for ±5%)
%
% Output:
%   sensitivityResults - Struct array with fields:
%       .param_index
%       .y_plus   (output vector after + perturbation)
%       .y_minus  (output vector after - perturbation)
%       .delta_y  (sensitivity per species)

    n_params = length(w_best);
    sensitivityResults = struct();

    for i = 1:n_params
        % Perturb parameter up and down
        w_plus  = w_best;
        w_minus = w_best;

        w_plus(i)  = w_best(i)  * (1 + perturb_frac);
        w_minus(i) = w_best(i)  * (1 - perturb_frac);

        % Simulate with perturbed parameters
        [~, ~, ~, y_long_plus, ~]  = NetfluxODE_Sim(w_plus);
        [~, ~, ~, y_long_minus, ~] = NetfluxODE_Sim(w_minus);

        % Store results
        sensitivityResults(i).param_index = i;
        sensitivityResults(i).y_plus  = y_long_plus;
        sensitivityResults(i).y_minus = y_long_minus;

        % Compute normalized sensitivity (change per unit perturbation)
        sensitivityResults(i).delta_y = ...
            (y_long_plus - y_long_minus) / (2 * perturb_frac * w_best(i));
    end

end