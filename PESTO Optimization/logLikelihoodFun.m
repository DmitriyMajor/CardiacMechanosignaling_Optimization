%% W OPT VERSION
function [ll, grad] = logLikelihoodFun(w_opt)

    try
        fprintf('--- logLikelihoodFun called ---\n');
        disp('Input w_opt:'); disp(w_opt');

        % Run simulation with current guess of weights
        [~, ~, validatedIdx, ylong_4h, y4h_norm] = NetfluxODE_Sim(w_opt);
        disp('Simulation completed');

        % Log lengths
        fprintf('Length validatedIdx: %d\n', length(validatedIdx));
        fprintf('Length ylong_4h: %d\n', length(ylong_4h));
        fprintf('Length y4h_norm: %d\n', length(y4h_norm));
        fprintf('Max(validatedIdx): %d\n', max(validatedIdx));

        % Check valid index range
        if any(validatedIdx > length(ylong_4h)) || any(validatedIdx > length(y4h_norm))
            error('validatedIdx exceeds vector bounds');
        end

       %% Data Check
        % Extract log2 fold-changes for validated genes
        % log2_data = log2(y4h_norm);
        log2_data = (y4h_norm);

        % Identify which validated genes are significantly changing in the data
        % significantIdx = validatedIdx(abs(log2_data(validatedIdx)) > 0.5);
        significantIdx = validatedIdx;

        % Now subset both model and data outputs to those significant genes
        % modelLog = log2(ylong_4h(significantIdx));
        % dataLog = log2_data(significantIdx);
        modelLog = (ylong_4h(significantIdx));
        dataLog = log2_data(significantIdx);

        % Apply filter
        n_sig = sum(significantIdx);
        if n_sig < 100  % You can adjust this threshold
            warning('Only %d validated genes passed the |log2FC| > 0.5 filter.', n_sig);
        end


        %% Objective Function Construction
        alpha = 1.5;  % prioritize mean
        beta  = 0.5;  % lesser priority on spread
        gamma = 1.0;

        rmse = sqrt(mean((modelLog - dataLog).^2));
        fprintf('RMSE (log2FC): %.5f\n', rmse);
        % mean_diff = (mean(modelLog) - mean(dataLog))^2;
        mean_diff = (abs(mean(modelLog) - mean(dataLog)));

        fprintf('Mean Diff^2 (log2FC): %.5f\n', mean_diff);
        % std_diff  = (std(modelLog) - std(dataLog))^2;
        std_diff  = (abs(std(modelLog) - std(dataLog)));
        fprintf('Stdev Diff: %.5f\n', std_diff);

        %AUC
        % [~, ~, AUC_down] = ComputeROC(log2(ylong_4h), log2(y4h_norm), significantIdx, 'down');
        % [~, ~, AUC_up] = ComputeROC(log2(ylong_4h), log2(y4h_norm), significantIdx, 'up');
        % 
        % fprintf('AUC UP: %.5f\n', AUC_up);
        % fprintf('AUC DOWN: %.5f\n', AUC_down);
        % 
        % threshold = 0.7;
        % activeIdxLocal = validatedIdx;
        % activeIdxLocal = validatedIdx((abs(log2(y4h_norm(validatedIdx))) > threshold) | (abs(log2(ylong_4h(validatedIdx))) > threshold) );
        % [accuracy, ~, ~, ~, ~, ~, ~, ~] = ... EvaluateGenePrediction(log2(y4h_norm), log2(ylong_4h), activeIdxLocal, threshold);
        % delta = 0.5;
        % epsilon = 0.5;

        % ll = accuracy;
        ll = - (alpha * mean_diff + beta * std_diff + gamma * rmse);
        % ll = -(alpha * mean_diff + beta * std_diff + gamma * rmse - delta * AUC_up - epsilon * AUC_down); %AUC
        grad = zeros(size(w_opt));      % still usable with PESTO



    catch ME
        warning(['logLikelihoodFun crashed: ', ME.message]);
        ll = 1e6;
        grad = zeros(size(theta));
    end
end
% 
% %% FACT OPT VERSION
% function [ll, grad] = logLikelihoodFun(theta)
%     try
%         fprintf('--- logLikelihoodFun called ---\n');
%         disp('Input theta:'); disp(theta');
% 
%         % Number of parameters you're optimizing (half EC50, half n)
%         n_params = length(theta) / 2;
%         EC50_opt = theta(1:n_params);
%         n_opt    = theta(n_params+1:end);
% 
%         % Constraint check: EC50 < (0.5)^(1/n)
%         for k = 1:n_params
%             if EC50_opt(k) >= (0.5)^(1 / n_opt(k))
%                 warning('Parameter constraint violated at index %d — EC50 = %.4f, n = %.4f', ...
%                         k, EC50_opt(k), n_opt(k));
%                 ll = 1e6;
%                 grad = zeros(size(theta));
%                 return;
%             end
%         end
% 
%         % Run simulation with current EC50 and n (pass them into simulation)
%         [~, ~, validatedIdx, ylong_4h, y4h_norm] = NetfluxODE_Sim(EC50_opt, n_opt);
%         validatedIdx = 85:654;
% 
%         disp('Simulation completed');
% 
%         % Log lengths
%         fprintf('Length validatedIdx: %d\n', length(validatedIdx));
%         fprintf('Length ylong_4h: %d\n', length(ylong_4h));
%         fprintf('Length y4h_norm: %d\n', length(y4h_norm));
%         fprintf('Max(validatedIdx): %d\n', max(validatedIdx));
% 
%         % Check valid index range
%         if any(validatedIdx > length(ylong_4h)) || any(validatedIdx > length(y4h_norm))
%             error('validatedIdx exceeds vector bounds');
%         end
% 
%         % Extract log2 fold-changes for validated genes
%         log2_data = log2(y4h_norm);
% 
%         % Significant genes filter (change if needed)
%         significantIdx = validatedIdx;
%         n_sig = length(significantIdx);
%         if n_sig < 100
%             warning('Only %d validated genes used.', n_sig);
%         end
% 
%         % Extract model and data values
%         modelLog = log2(ylong_4h(significantIdx));
%         dataLog = log2_data(significantIdx);
% 
%         % Error metrics
%         rmse = sqrt(mean((modelLog - dataLog).^2));
%         mean_diff = abs(mean(modelLog) - mean(dataLog));
%         std_diff  = abs(std(modelLog) - std(dataLog));
% 
%         % KS test statistic (two-sample)
%         [~, ~, ks_stat] = kstest2(modelLog, dataLog);
% 
%         threshold = 0.1;
%         penalty_zero = sum(abs(modelLog) < threshold) / length(modelLog);
%         fprintf('Zero-region penalty: %.5f\n', penalty_zero);   
%         fprintf('KS Statistic: %.5f\n', ks_stat);
%         fprintf('RMSE: %.5f, Mean Diff: %.5f, Stdev Diff: %.5f\n', rmse, mean_diff, std_diff);
% 
%         % Combined objective
%         alpha = 1.5;  % Weight for mean difference
%         beta  = 1.0;  % Weight for stdev difference
%         gamma = 1.0;  % Weight for RMSE
%         delta = 1.5;  % Weight for KS Test
%         epsilon = 3.0;  % Weight for near-zero penalty
% 
% ll = -(alpha * mean_diff + beta * std_diff + gamma * rmse + delta * ks_stat + epsilon * penalty_zero);
%         % Placeholder zero gradient
%         grad = zeros(size(theta));
% 
%     catch ME
%         warning(['logLikelihoodFun crashed: ', ME.message]);
%         ll = 1e6;
%         grad = zeros(size(theta));
%     end
% end

function [TPR, FPR, AUC] = ComputeROC(modelLog2, dataLog2, idx, direction)
    thresholds = linspace(0, 2, 50);
    TPR = zeros(size(thresholds));
    FPR = zeros(size(thresholds));
    mVals = modelLog2(idx);
    dVals = dataLog2(idx);
    for i = 1:length(thresholds)
        th = thresholds(i);
        if strcmp(direction, 'up')
            model_pred = mVals >= th;
            data_true = dVals >= 0.5;
        elseif strcmp(direction, 'down')
            model_pred = mVals <= -th;
            data_true = dVals <= -0.5;
        end
        TP = sum(model_pred & data_true);
        FP = sum(model_pred & ~data_true);
        FN = sum(~model_pred & data_true);
        TN = sum(~model_pred & ~data_true);
        TPR(i) = TP / (TP + FN + eps);
        FPR(i) = FP / (FP + TN + eps);
    end
% Calculate AUC using trapezoidal rule
[sortedFPR, sortIdx] = sort(FPR);
sortedTPR = TPR(sortIdx);
AUC = trapz(sortedFPR, sortedTPR);
end

function [accuracy, TP, TN, FP, FN, precision, recall, f1_score] = ...
    EvaluateGenePrediction(log2data, log2model, validatedIdx, threshold)

    % Extract validated values
    data_vals = log2data(validatedIdx);
    model_vals = log2model(validatedIdx);

    % Binarize based on threshold
    dataClass = abs(data_vals) > threshold;
    modelClass = abs(model_vals) > threshold;

    % Compute confusion matrix
    [confMat, order] = confusionmat(dataClass, modelClass);

    % Extract values safely
    % Order: [0, 1] = [non-sig, sig]
    TP = 0; TN = 0; FP = 0; FN = 0;
    if isequal(order, [0; 1])
        TN = confMat(1,1);
        FP = confMat(1,2);
        FN = confMat(2,1);
        TP = confMat(2,2);
    elseif isequal(order, [1; 0])
        TP = confMat(1,1);
        FN = confMat(1,2);
        FP = confMat(2,1);
        TN = confMat(2,2);
    else
        warning('Unexpected class order in confusion matrix.');
    end

    % Compute metrics
    total = TP + TN + FP + FN;
    accuracy = (TP + TN) / total;
    precision = TP / (TP + FP + eps);
    recall = TP / (TP + FN + eps);
    f1_score = 2 * (precision * recall) / (precision + recall + eps);
end