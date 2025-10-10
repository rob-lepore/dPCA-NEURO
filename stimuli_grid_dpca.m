function [nmsc_3x3] = stimuli_grid_dpca(data, W, timeBegin, timeEnd, compIdx, global_neuron_mean, g_means, g_stds)
%STIMULI_GRID_DPCA Projects data for a given epoch and component, 
%normalizes using global dPCA statistics, and reshapes into a 3×3 grid.
%
% Inputs:
%   data               - neurons × stimuli × time matrix (e.g., firingRatesAverage)
%   W                  - dPCA projection matrix (neurons × nComps)
%   timeBegin, timeEnd - start/end time indices for the epoch
%   compIdx            - component index to visualize (1..nComps)
%   global_neuron_mean - mean firing rate per neuron (1×neurons)
%   g_means, g_stds    - global per-component mean and std from full dPCA projection
%
% Output:
%   nmsc_3x3           - 3×3 matrix of normalized mean spike counts

    epochData = data(:,:,timeBegin:timeEnd); 
    [~, num_stimuli, num_time_points] = size(epochData);
    nComps = size(W,2);

    X = epochData(:,:)';       
    Xcen = bsxfun(@minus, X, global_neuron_mean);   

    dPCs = Xcen * W;                                  
    dPCs = reshape(dPCs', nComps, num_stimuli, num_time_points);

    dPC_comp = squeeze(dPCs(compIdx,:,:));          
    meanSpikeCount = mean(dPC_comp, 2);             

    normalizedMeanSpikeCount = ...
        (meanSpikeCount - g_means(compIdx)) / g_stds(compIdx);

    nmsc_3x3 = flipud(reshape(normalizedMeanSpikeCount, [3, 3])');

end
