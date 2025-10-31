function [angles_deg] = full_dpca_analysis(data_Path, hand_position)

% File to load/save dPCA results
data_save_file = "my_dpca_data";
load_from_file = false; % If true, skip dPCA and load from file above

epochNames = {'Fixation', 'Plan', 'Reach', 'Hold'};


if ~load_from_file
    cells_in_Directory = dir(data_Path);
    cells_in_Directory ([1,2],:) = [];

    epochEvent = {'Saccade-Off', 'GO', 'KeyUp', 'TOUCH1'};
    epochTimes = [0 700; 500 200; 200 500; 0 700];   
    sDF_bin_Size = 100;              
    
    numEpochs = numel(epochEvent);
    firingRates_all = cell(1, numEpochs);

    for e = 1:numEpochs
        event_Name = epochEvent{e};
        time_Window = epochTimes(e, :);
    
        [firingRates, trialNum] = A_general_calculate_firing_rates_dpca( ...
            data_Path, cells_in_Directory, time_Window, sDF_bin_Size, event_Name);
    
        % Keep only one hand position
        firingRates = squeeze(firingRates(:,:,hand_position,:,:));   % neurons × stimuli × time × trials
        % Average over trials
        firingRatesAverage = mean(firingRates, 4);        % neurons × stimuli × time
        firingRates_all{e} = firingRatesAverage;
    end

    firingRates_dpca = cat(4, firingRates_all{:});
    firingRates_dpca = permute(firingRates_dpca, [1 2 4 3]);
        
    combinedParams = {{1, [1 3]}, {2, [2 3]}, {[1 2], [1 2 3]}};

    num_comp = 35;

    [W,V,whichMarg] = dpca(firingRates_dpca, num_comp, ...
        'combinedParams', combinedParams);

    explVar = dpca_explainedVariance(firingRates_dpca, W, V, ...
    'combinedParams', combinedParams);

    % Decode into dPCs
    X = firingRates_dpca(:,:)';
    global_neuron_mean = mean(X, 1);
    Xcen = bsxfun(@minus, X, global_neuron_mean);
    dPCs = Xcen * W;
    g_means = mean(dPCs, 1);                               
    g_stds  = std(dPCs, 0, 1);

    if data_save_file~=""
        save(data_save_file, "W", "V", "whichMarg","firingRates_all",...
            "explVar", "g_means", "g_stds", "global_neuron_mean");
    end
else
    if data_save_file~=""
        load(data_save_file);
    else
        disp(['Error: choose a file to load data from']);
        return;
    end
end


timeMarg = find(whichMarg == 2);
timeMarg = timeMarg(1:3);
timeMargExplVar = explVar.componentVar(timeMarg);

[num_neurons, num_stimuli, num_time_points] = size(firingRates_all{1});
num_comp = size(W,2);

%% Gradient maps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

marginalizations = [find(whichMarg == 1, 2, 'first'), find(whichMarg == 3, 2, 'first')];
marginalizationsExplVar = explVar.componentVar(marginalizations);

nRows = length(marginalizations);  % number of components
nCols = length(epochNames); % number of epochs

figure;
nInterp = 100;           
pad_frac = 0.5;          
colorRange = [-2 2];  

for j = 1:nRows
    marg = marginalizations(j);

    for i = 1:nCols

        normalized = normalizedEpochSpikes(firingRates_all,W,i,marg,global_neuron_mean, g_means, g_stds);

        % Arrange into 3×3 LED grid
        nmsc_3x3 = reshape(normalized, [3, 3])';
        nmsc_3x3 = flipud(nmsc_3x3);


        [X, Y] = meshgrid(1:3, 1:3);
        [XX, YY] = meshgrid(linspace(1,3,nInterp), linspace(1,3,nInterp));
        interp_values = interp2(X, Y, nmsc_3x3, XX, YY, 'linear');

        subplot(nRows, nCols, (j-1)*nCols + i);
        hold on;

        imagesc([1 3], [1 3], interp_values, colorRange); 
        axis square;
        set(gca,'XTick',1:3,'YTick',1:3,'XTickLabel',{},'YTickLabel',{});
        set(gca,'TickLength',[0 0]);

        % Add equal padding around gradient
        xlim([1-pad_frac 3+pad_frac]);
        ylim([1-pad_frac 3+pad_frac]);

        plot(X(:), Y(:), 'k+', 'MarkerSize', 7, 'LineWidth', 1.5);

        if j == 1
            title(epochNames{i}, 'FontSize', 11);
        end
        if i == 1
            ylabel(sprintf('dPC #%d [%.1f%%]', marg, marginalizationsExplVar(j)), 'FontSize', 11, 'FontWeight', 'bold');
        end

        colormap(parula);
        clim(colorRange);
        hold off;
    end
end

cb = colorbar('Position', [0.93, 0.1, 0.015, 0.8]);
ylabel(cb, 'Normalized mean spike count', 'FontSize', 11);


%% Exaplained variance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bar plot with projected variances
figure;
axBar = subplot(1,1,1);
hold on
num_comp_bars = 15;
axis([0 num_comp_bars+1 0 Inf])
ylabel('Component variance (%)')
xlabel('Number of components')
b = bar(explVar.margVar(:,1:num_comp_bars)' , 'stacked', 'BarWidth', 0.75);

margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

for idx = 1:numel(b)
    b(idx).FaceColor = margColours(idx,:);
end 

clim([1 length(margColours)+256]);
set(gca, 'FontSize', 14);

% Pie chart with total explained variance
figure;
axPie = subplot(1,1,1);

d = explVar.totalMarginalizedVar / explVar.totalVar * 100;
roundedD = floor(d);
while sum(roundedD) < 100
    [~, ind] = max(d-roundedD);
    roundedD(ind) = roundedD(ind) + 1;
end

margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;
margNames = {'Eye', 'Context', 'Interaction'};

for i=1:length(d)
    margNamesPerc{i} = [margNames{i} ' ' num2str(roundedD(i)) '%'];
end

p = pie(d, ones(size(d)));

for k = 1:length(p)
    if strcmp(get(p(k), 'Type'), 'text')
        set(p(k), 'String', ''); 
    end
end

numSegments = length(d);
numColours = size(margColours, 1);
for k = 1:numSegments
    if strcmp(get(p(2*k-1), 'Type'), 'patch') 
        colorIndex = mod(k-1, numColours) + 1; 
        set(p(2*k-1), 'FaceColor', margColours(colorIndex, :));
    end
end

lgd = legend(margNamesPerc, 'Location', 'bestoutside');
set(lgd, 'FontSize', 14);


%% Scatter Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scatterdata = zeros(length(timeMarg), length(epochNames), num_stimuli);

for j = 1:length(timeMarg)
    compIdx = timeMarg(j);   % which dPC to project onto

    for ei = 1:length(epochNames)
        normalized = normalizedEpochSpikes(firingRates_all,W,ei,compIdx,global_neuron_mean, g_means, g_stds);
        scatterdata(j, ei, :) = normalized(:)';
    end
end


colors = [
    0.8500, 0.3250, 0.0980;  
    0.9290, 0.6940, 0.1250;  
    0.4940, 0.1840, 0.5560;  
    0.4660, 0.6740, 0.1880;  
    0.3010, 0.7450, 0.9330;  
    0.6350, 0.0780, 0.1840;  
    0.0000, 0.4470, 0.7410;  
    0.8500, 0.3250, 0.7980;  
    0.8, 0.6, 0.4;           
];


% fig = figure("Position",[0 0 1660 468]);
figure;

nPlots = length(timeMarg);
for plotIdx = 1:nPlots
    subplot(1, nPlots, plotIdx);
    hold on;
    for catIdx = 1:4
        for valIdx = 1:9
            value = scatterdata(plotIdx, catIdx, valIdx);
            scatter(catIdx, value, 80, colors(valIdx, :), 'filled');
        end
    end
    xticks(1:4);  
    xticklabels(epochNames);
    yticks(floor(min(scatterdata(:))):ceil(max(scatterdata(:))));
    title(sprintf("dPC #%d [%.1f%%]", timeMarg(plotIdx),  timeMargExplVar(plotIdx)))
    ylim([-3, 3]);
    if plotIdx == 1
        ylabel("Norm. spk count");
    end

    axis square;
    set(gca, 'FontSize', 14);
    hold off;
end

% hold off;


% Scatter legend

x = [1, 2, 3]; % 1: Left, 2: Center, 3: Right
y = [1, 2, 3]; % 1: Near, 2: Intermediate, 3: Far
figure;
hold on;

index = 1;
for i = 1:3
    for j = 1:3
        plot(x(j), y(i), 'o', 'MarkerSize', 10, 'MarkerEdgeColor', colors(index, :), 'MarkerFaceColor', colors(index, :));
        index = index + 1;
    end
end

xlim([0.5 3.5]);
ylim([0.5 3.5]);

set(gca, 'XTick', [1, 2, 3], 'XTickLabel', {'Left', 'Center', 'Right'});
set(gca, 'YTick', [1, 2, 3], 'YTickLabel', {'Near', 'Intermediate', 'Far'});

grid on;
hold off;


%% Subspace analysis
eyeIdx     = find(whichMarg == 1);
eyeIdx = eyeIdx(1:8);
contextIdx = find(whichMarg == 2);
contextIdx = contextIdx(1:3);

W_eye     = W(:, eyeIdx);      % neurons × 8
W_context = W(:, contextIdx);  % neurons × 3

angles = subspacea(W_eye, W_context);
angles_deg = angles * 180 / pi;
