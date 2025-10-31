function [normalized] = normalizedEpochSpikes(data, W, epochId, margId, global_neuron_mean, g_means, g_stds)
    epochData = data{epochId};     % neurons × stimuli × time
    [num_neurons, num_stimuli, num_time_points] = size(epochData);

    X = reshape(epochData, num_neurons, [])';   % (stimuli*time) × neurons

    % Center by global neuron mean
    Xcen = bsxfun(@minus, X, global_neuron_mean);

    % Project onto selected dPC
    dPC_comp = Xcen * W(:, margId);           % (stimuli*time) × 1

    % Reshape back to (time × stimuli)
    dPC_comp = reshape(dPC_comp, num_stimuli, num_time_points);
    meanSpikeCount = mean(dPC_comp, 2);

    % Normalize using global means and stds
    normalized = (meanSpikeCount - g_means(margId)) ./ g_stds(margId);

end
