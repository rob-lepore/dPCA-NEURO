function [nmsc_3x3] = stimuli_grid_dpca(data, W, timeBegin, timeEnd, marg)
   
epochData = data(:,:,:,timeBegin:timeEnd);
[num_neurons, num_stimuli, num_decisions, num_time_points] = size(epochData);
num_comp = 15;
X = epochData(:,:)';
Xcen = bsxfun(@minus, X, mean(X));
dPCs = Xcen * W;
dPCs = reshape(dPCs', num_comp, num_stimuli, num_decisions, num_time_points);

dPCmarg = squeeze(dPCs(marg, :, :, :));
dPCmarg_near = squeeze(dPCmarg(:, 1, :));

meanSpikeCount = mean(dPCmarg_near, 2); % mean for each stimulus wrt time 
mean_val = mean(meanSpikeCount);
std_val = std(meanSpikeCount);
normalizedMeanSpikeCount = (meanSpikeCount - mean_val) / std_val;

% Arrange stimuli by spatial disposition
nmsc_3x3 = flipud(reshape(normalizedMeanSpikeCount, [3, 3])');
