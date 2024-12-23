function dpca_plot_default(data, time, yspan, explVar, compNum, events, signif, ~)

% Modify this function to adjust how components are plotted.
%
% Parameters are as follows:
%   data      - data matrix, size(data,1)=1 because it's only one component
%   time      - time axis
%   yspan     - y-axis spab
%   explVar   - variance of this component
%   compNum   - component number
%   events    - time events to be marked on the time axis
%   signif    - marks time-point where component is significant
%   marg      - marginalization number


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% displaying legend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%panel positions
panel_pos=["F2 - near sx" "F3 - near centre" "F4 - near dx" "F7 - interm. sx" "F8 - interm. centre " "F9 - interm. dx" "F12 - far sx" "F13 - far centre" "F14 - far dx"];
%1290	1292	1296	1304	1306	1308	1305	1307	1310
% F2    F3      F4      F7      F8      F9      F12     F13     F14


if strcmp(data, 'legend')
    % if there is only time and no other parameter - do nothing
    if length(time) == 2
        return

        % if there is one parameter
    elseif length(time) == 3
        numOfStimuli = time(2); % time is used to pass size(data) for legend
        %colors = lines(numOfStimuli);
        %colors = colormap(brighten(hsv(9), -0.5));
        colors = [5 0 149; 28 177 203; 10 245 215; 25 218 7; 235 218 0; 245 137 13; 246 11 0; 225 72 168; 153 0 226]/256;

        hold on

        for f = 1:numOfStimuli
            plot([0.5 1], [f f], 'color', colors(f,:), 'LineWidth', 2)
            text(1.2, f, ['Stimulus ' num2str(f)])
        end
        axis([0 3 -1 1.5+numOfStimuli])
        set(gca, 'XTick', [])
        set(gca, 'YTick', [])
        set(gca,'Visible','off')
        return

    % two parameters: stimulus and decision (decision can only have two values)
    elseif length(time) == 4 && time(3) == 2
        numOfStimuli = time(2); % time is used to pass size(data) for legend
        %colors = lines(numOfStimuli);
        %colors = colormap(brighten(hsv(9), -0.5));
        colors = [5 0 149; 28 177 203; 10 245 215; 25 218 7; 235 218 0; 245 137 13; 246 11 0; 225 72 168; 153 0 226]/256;

        hold on

        for f = 1:numOfStimuli
            plot([0.5 1], [f f], 'color', colors(f,:), 'LineWidth', 2)
            text(1.2, f, ['Stimulus ' num2str(f)])
        end

        plot([0.5 1], [-2 -2], 'k', 'LineWidth', 2)
        plot([0.5 1], [-3 -3], 'k--', 'LineWidth', 2)
        text(1.2, -2, 'Hand near')
        text(1.2, -3, 'Hand far')

        axis([0 3 -4.5 1.5+numOfStimuli])
        set(gca, 'XTick', [])
        set(gca, 'YTick', [])
        set(gca,'Visible','off')
        return

        % other cases - do nothing
    else
        return
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting up the subplot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(time)
    time = 1:size(data, ndims(data));
end
%setting the axis !!
axis([time(1) time(end) yspan])

hold on

if ~isempty(explVar)
    title(['Component #' num2str(compNum) ' [' num2str(explVar,'%.1f') '%]'])
else
    title(['Component #' num2str(compNum)])
end

%plot line on the alignment event
%events =['KeyDown'; 'LED1'; 'Saccade-Off'; 'GO'; 'KeyUp'; 'TOUCH1'; 'RedOff'; 'TOUCH2'; 'keydown']




if ~isempty(events)
    plot([events; events], yspan, 'Color', [0.6 0.6 0.6])
end

if ~isempty(signif)
    signif(signif==0) = nan;
    plot(time, signif + yspan(1) + (yspan(2)-yspan(1))*0.05, 'k', 'LineWidth', 3)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting the component
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ndims(data) == 2
    % only time - plot it
    plot(time, squeeze(data(1, :)), 'k', 'LineWidth', 2)

elseif ndims(data) == 3
    % different stimuli in different colours
    numOfStimuli = size(data, 2);
    plot(time, squeeze(data(1,:,:)), 'LineWidth', 2)

elseif ndims(data) == 4 && size(data,3)==2
    %!! different stimuli in different colours and binary condition as
    % solid/dashed
    numOfStimuli = size(data, 2);
    %colors = lines(numOfStimuli);
    %colors = colormap(brighten(hsv(9), -0.5));
    colors = [5 0 149; 28 177 203; 10 245 215; 25 218 7; 235 218 0; 245 137 13; 246 11 0; 225 72 168; 153 0 226]/256;

    for f=1:numOfStimuli
        %hand near
        plot(time, squeeze(data(1, f, 1, :)), 'color', colors(f,:), 'LineWidth', 2)
        %hand far - same color for same postiiton but dashed line
        plot(time, squeeze(data(1, f, 2, :)), '--', 'color', colors(f,:), 'LineWidth', 2)
        %lines returns the lines colormap as a three-column array with the same number of
        %rows as the colormap for the current figure.
        % Each row in the array contains the red, green, and blue intensities for a specific color.
        % The intensities are in the range [0,1], and the color scheme matches the default
        % ColorOrder property of the Axes.
    end
else
    % in all other cases pool all conditions and plot them in different
    % colours
    data = squeeze(data);
    dims = size(data);
    data = permute(data, [numel(dims) 1:numel(dims)-1]);
    data = reshape(data, size(data,1), []);
    data = data';
    plot(time, data, 'LineWidth', 2)
end