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

        % Extract log2 fold-changes for validated genes
        log2_data = log2(y4h_norm);

        % Identify which validated genes are significantly changing in the data
        % significantIdx = validatedIdx(abs(log2_data(validatedIdx)) > 0.5);
        significantIdx = validatedIdx;

        % Now subset both model and data outputs to those significant genes
        modelLog = log2(ylong_4h(significantIdx));
        dataLog = log2_data(significantIdx);

        % Apply filter
        n_sig = sum(significantIdx);
        if n_sig < 100  % You can adjust this threshold
            warning('Only %d validated genes passed the |log2FC| > 0.5 filter.', n_sig);
        end


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

        ll = - (alpha * mean_diff + beta * std_diff + gamma * rmse);
        grad = zeros(size(w_opt));      % still usable with PESTO



    catch ME
        warning('logLikelihoodFun crashed: %s', ME.message);
        ll = 1e6;
        grad = zeros(size(w_opt));  % Also handle crash case
    end
end