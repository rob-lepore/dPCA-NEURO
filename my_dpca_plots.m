clear;
clc;
close all;

%% From D_dpca_Reaching.m *************************************************

% Get data
data_save_file = "my_dpca_data"; % file to load/save dPCA results
load_from_file = true; % true if data is already in saved in a file

if ~load_from_file
    data_Path = 'C:\Users\rober\OneDrive\Desktop\Uni\Bioinformatics\Neurosciences\project\V6a_mef24\V6A_mef24\';
    cells_in_Directory = dir(data_Path);
    cells_in_Directory ([1,2],:) = [];
    
    event_Name = 'Saccade-Off';    
    dataname = 'V6A_mef24 - Saccade-Off';
    time_Window = [200,3200];
    sDF_bin_Size = 100;              
    ifSimultaneousRecording = true;  
   
    % Compute firing rates
    [firingRates, trialNum] = A_general_calculate_firing_rates_dpca(data_Path, cells_in_Directory, time_Window, sDF_bin_Size, event_Name);
    firingRatesAverage = mean(firingRates, 5); %nanmean

    % dPCA
    disp(['dPCA with regularization'])
    
    % Time is marginalization #3
    combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
    num_comp = 15;
    
    optimalLambda = dpca_optimizeLambda(firingRatesAverage, firingRates, trialNum, dataname, ...
        'combinedParams', combinedParams, ...
        'simultaneous', ifSimultaneousRecording, ...
        'numRep', 2, ...  % (2) increase this number to ~10 for better accuracy ***
        'filename', 'tmp_optimalLambdas.mat', ...
        'display', false);
    
    Cnoise = dpca_getNoiseCovariance(firingRatesAverage, ...
        firingRates, trialNum, 'simultaneous', ifSimultaneousRecording);
    
    [W,V,whichMarg] = dpca(firingRatesAverage, num_comp, ... % numComp=20
        'combinedParams', combinedParams, ...
        'lambda', optimalLambda, ...
        'Cnoise', Cnoise);

    explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
    'combinedParams', combinedParams);

    if data_save_file~=""
        save(data_save_file, "W", "whichMarg","firingRatesAverage", "num_comp", "explVar");
    end
else
    if data_save_file~=""
        load(data_save_file);
    else
        disp(['Error: choose a file to load data from']);
        return;
    end
end

%**************************************************************************

% Select dPCs relative to time marginalization
timeMarg = find(whichMarg == 3); 
timeMargExplVar = explVar.componentVar(timeMarg);

[num_neurons, num_stimuli, num_decisions, num_time_points] = size(firingRatesAverage);

%% From dpca_plot.m *******************************************************
X = firingRatesAverage(:,:)';
Xcen = bsxfun(@minus, X, mean(X));
dPCs = Xcen * W;
dPCs = reshape(dPCs', num_comp, num_stimuli, num_decisions, num_time_points);
%**************************************************************************

%% Plot normalized mean spike count on each time dPC **********************
fig = figure('Position', [100,-50,600,591]);
epochNames = {'Fixation', 'Plan', 'Reach', 'Hold'};
epochTimes = [1 1000; 1631 2331; 2332 2592; 2593 3392];
timeMarg = find(whichMarg == 3); 

for j=1:length(timeMarg)
    marg = timeMarg(j);
    for i=1:length(epochTimes)
        epochName = epochNames(i);
        epochTime = epochTimes(i,:);
        nmsc_3x3 = stimuli_grid_dpca(firingRatesAverage, W, epochTime(1), epochTime(2), marg);
        if j==1
            plotTitle = epochNames(i);
        else
            plotTitle = "";
        end
        gradientmap(nmsc_3x3, 'title', plotTitle, 'subp_params', [4 5 (j-1)*5+i], 'withCb', false);
    end
end

subplot(3,5,5);
custom_cmap = [0 0 255; 115 164 255; 128 255 183; 251 218 83; 255 255 0] / 255; 
n_colors = 256;
cmap = interp1(linspace(1, n_colors, size(custom_cmap, 1)), custom_cmap, 1:n_colors, 'linear');
cb = colorbar();
colormap(cmap);
caxis([0,1]);
set(cb, 'Ticks', [0,1]);
set(cb, 'TickLabels', {'-', '+'});
ylabel(cb, sprintf('Normalized mean\nspike count'));
cb.Label.Rotation = -90;
ax = gca;
pos = ax.Position;
axis off;
cb.Position = [pos(1) pos(2) 0.012 pos(4)]; 
cb.Label.Position = [1.809776013547724,0.498269949127866,0];

%exportgraphics(fig, 'gradientMaps.pdf', 'ContentType', 'vector');
%set(fig, 'Color', 'none');

%print(fig, 'gradientMaps.png', '-dpng', '-r500', '-opengl');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bar plot with projected variances %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = figure('Position', [100,-50,645,328]);
axBar = subplot(1,1,1);
hold on
axis([0 num_comp+1 0 12.5])
ylabel('Component variance (%)')
xlabel('Number of components')
b = bar(explVar.margVar(:,1:num_comp)' , 'stacked', 'BarWidth', 0.75);

margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

for idx = 1:numel(b)
    b(idx).FaceColor = margColours(idx,:);
end 

caxis([1 length(margColours)+256]);
set(gca, 'FontSize', 14);

%print(fig, 'barPlot.png', '-dpng', '-r500', '-image');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pie Chart %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = figure;
axPie = subplot(1,1,1);

d = explVar.totalMarginalizedVar / explVar.totalVar * 100;
roundedD = floor(d);
while sum(roundedD) < 100
    [~, ind] = max(d-roundedD);
    roundedD(ind) = roundedD(ind) + 1;
end

margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;
margNames = {'Target', 'Hand', 'Condition-independent', 'T/H Interaction'};

for i=1:length(d)
    margNamesPerc{i} = [margNames{i} ' ' num2str(roundedD(i)) '%'];
end

p = pie(d, ones(size(d)));

for k = 1:length(p)
    if strcmp(get(p(k), 'Type'), 'text')
        set(p(k), 'String', ''); % Hide the text
    end
end

% Set colors
numSegments = length(d);
numColours = size(margColours, 1);
for k = 1:numSegments
    if strcmp(get(p(2*k-1), 'Type'), 'patch') % Handle patch segments
        colorIndex = mod(k-1, numColours) + 1; % Loop through colors if needed
        set(p(2*k-1), 'FaceColor', margColours(colorIndex, :));
    end
end

lgd = legend(margNamesPerc, 'Location', 'bestoutside');
set(lgd, 'FontSize', 14);

%print(fig, 'pieChart.png', '-dpng', '-r500', '-image');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scatter Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g_sc = [];
g_means = [];
g_stds = [];
for j=1:length(timeMarg)
    marg = timeMarg(j);
    all_meanSpikeCount = [];
    for i=1:length(epochNames)
        epochData = firingRatesAverage(:,:,:,epochTimes(i,1):epochTimes(i,2));
        [num_neurons, num_stimuli, num_decisions, num_time_points] = size(epochData);
        num_comp = 15;
        X = epochData(:,:)';
        Xcen = bsxfun(@minus, X, mean(X));
        dPCs = Xcen * W;
        dPCs = reshape(dPCs', num_comp, num_stimuli, num_decisions, num_time_points);
        dPCmarg = squeeze(dPCs(marg, :, :, :));
        dPCmarg_near = squeeze(dPCmarg(:, 1, :));
        meanSpikeCount = mean(dPCmarg_near, 2);
        all_meanSpikeCount = [all_meanSpikeCount; meanSpikeCount]; 
    end
    % mean and std for each dPC
    g_means = [g_means; mean(all_meanSpikeCount)];
    g_stds = [g_stds; std(all_meanSpikeCount)];
    g_sc = [g_sc; all_meanSpikeCount];
end


scatterdata = [];
for j = 1:length(timeMarg)
    marg = timeMarg(j);
    
    for i = 1:length(epochNames)
        s = ((j-1)*36+1)+(i-1)*9;
        e = ((j-1)*36+1)+i*9-1;
        normalizedMeanSpikeCount = (g_sc(s:e) - g_means(j)) / g_stds(j);
        scatterdata(j, i, :) = normalizedMeanSpikeCount;%normalizedMeanSpikeCount; % Store normalized values
    end
end

x = 1:length(epochNames);

colors = [
    0.8500, 0.3250, 0.0980;  
    0.9290, 0.6940, 0.1250;  
    0.4940, 0.1840, 0.5560;  
    0.4660, 0.6740, 0.1880;  
    0.3010, 0.7450, 0.9330;  
    0.6350, 0.0780, 0.1840;  
    0.0000, 0.4470, 0.7410;  
    0.8500, 0.3250, 0.0980;  
    0.8, 0.6, 0.4;           
];


fig = figure("Position",[0 0 1660 468]);


for plotIdx = 1:4
    subplot(1, 4, plotIdx);
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
    title(sprintf("dPC #%d [%.1f %%]", timeMarg(plotIdx),  timeMargExplVar(plotIdx)))
    ylim([-3, 3]);
    if plotIdx == 1
        ylabel("Norm. spk count");
    end

    axis square;
    set(gca, 'FontSize', 14);
    hold off;
end

hold off;
%print(fig, 'scatterPlot.png', '-dpng', '-r500', '-image');


%% scatter legend
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

% Set the axis limits
xlim([0.5 3.5]);
ylim([0.5 3.5]);

% Set the axis labels
set(gca, 'XTick', [1, 2, 3], 'XTickLabel', {'Left', 'Center', 'Right'});
set(gca, 'YTick', [1, 2, 3], 'YTickLabel', {'Near', 'Intermediate', 'Far'});

% Add grid lines
grid on;

hold off

