close all;
clc;

dPC_comb = [3 5; 3 8; 3 12; 5 8; 5 12; 8 12];

% Time trajectory for a single stimulus
figure;
for i = 1:size(dPC_comb, 1)
    dPC_x = dPC_comb(i, 1);
    dPC_y = dPC_comb(i, 2);
    
    subplot(2, 3, i);
    hold on;
    stimulus = 5; 
    dPC1 = squeeze(dPCs(dPC_x, stimulus, 2, :));
    dPC2 = squeeze(dPCs(dPC_y, stimulus, 2, :));
    
    scatter(dPC1, dPC2, 15, 1:3400, 'filled');
    colorbar;
    xlabel(['dPC #' num2str(dPC_x)]);
    ylabel(['dPC #' num2str(dPC_y)]);
    grid on;
    hold off;
end
sgtitle('Time trajectories for a single stimulus')

% Similarity in time trajectory for different stimuli
figure;
for i = 1:size(dPC_comb, 1)
    dPC_x = dPC_comb(i, 1);
    dPC_y = dPC_comb(i, 2);
    
    subplot(2, 3, i);
    hold on;
    colors = lines(9);
    
    for stimulus = 1:9
        dPC1 = squeeze(dPCs(dPC_x, stimulus, 2, 1:10:end));
        dPC2 = squeeze(dPCs(dPC_y, stimulus, 2, 1:10:end));
        scatter(dPC1, dPC2, 15, 'filled', 'MarkerFaceColor', colors(stimulus, :), 'MarkerFaceAlpha', 0.6);
    end
    
    xlabel(['dPC #' num2str(dPC_x)]);
    ylabel(['dPC #' num2str(dPC_y)]);
    legend(arrayfun(@(x) ['Stimolo ' num2str(x)], 1:9, 'UniformOutput', false), 'Location', 'bestoutside');
    grid on;
    hold off;
end
sgtitle('Time trajectories of all stimuli')



figure

hold on;
colors = lines(9);
stimulus = 5; 
dPC1 = squeeze(dPCs(1, stimulus, 2, :));
dPC2 = squeeze(dPCs(2, stimulus, 2, :));

scatter(dPC1, dPC2, 15, 1:3400, 'filled');
colorbar;
xlabel(['dPC #1']);
ylabel(['dPC #2']);
grid on;
hold off;