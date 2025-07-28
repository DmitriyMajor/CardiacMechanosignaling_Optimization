%% Functions
% Normalization of genes by SS value (FC)

function [yout_4h, yout_30m, ynorm] = ODENorm(t_in, y_in, ss, geneIdx)
ynorm = y_in;
ynorm(:, geneIdx) = y_in(:, geneIdx) ./ ss(geneIdx)';
yout_4h = ynorm(find(t_in >= 240, 1), :)';
yout_30m = ynorm(find(t_in >= 30, 1), :)';
end
