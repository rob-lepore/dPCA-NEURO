function gradientmap(data, varargin)
% Function to plot a 2D matrix of numerical values as a smooth gradient
% map.
% Optional parameters:
% - subp_params: In case the gradient map has to be displayed as a subplot,
%                specifies the size of the subplot grid and the position of
%                the plot, as required by the 'subplot' function. The
%                value [1 1 1] (default) creates a single figure.
% - title: Title of the figure or subplot.
% - X_tick_labels: Tick labels along the X axis.
% - Y_tick_labels: Tick labels along the Y axis.
% - X_label: Label of the X axis.
% - Y_label: Label of the Y axis.




    % default input parameters
    options = struct('subp_params',       [1 1 1], ...
                     'title', "", ...
                     'X_tick_labels', [], ...
                     'Y_tick_labels', [], ...
                     'X_label', '', ...
                     'Y_label', '',...
                     'withCb', true);

    % read input parameters
    optionNames = fieldnames(options);
    if mod(length(varargin),2) == 1
	    error('Please provide propertyName/propertyValue pairs')
    end
    for pair = reshape(varargin,2,[])
	    if any(strcmp(pair{1}, optionNames))
            options.(pair{1}) = pair{2};
        else
            error('%s is not a recognized parameter name', pair{1})
	    end
    end


    subp_params = options.subp_params;
    if all(subp_params==1)
        figure;
    else
        subplot(subp_params(1),subp_params(2),subp_params(3));
    end


    % cross points
    [size_y, size_x] = size(data);
    values = data(:);
    y_pos = 1:size_y; 
    x_pos = 1:size_x; 

    % custom colors
    custom_cmap = [0 0 255; 115 164 255; 128 255 183; 251 218 83; 255 255 0] / 255; 
    n_colors = 256;
    cmap = interp1(linspace(1, n_colors, size(custom_cmap, 1)), custom_cmap, 1:n_colors, 'linear');

    % gradient grid
    norm_values = (values - min(values)) / (max(values) - min(values));
    color_indices = floor(norm_values * (size(cmap, 1) - 1)) + 1;
    colors_grid = reshape(cmap(color_indices, :), [size_y, size_x, 3]);

    [X, Y] = meshgrid(x_pos, y_pos);
    [XX, YY] = meshgrid(linspace(1, size_x, 100), linspace(1, size_y, 100));
    gradient_colors = zeros(100, 100, 3);

    for i = 1:3
        gradient_colors(:,:,i) = griddata(X, Y, squeeze(colors_grid(:,:,i)), XX, YY, 'natural');
    end

    imagesc([1, size_x], [1, size_y], gradient_colors);
    axis square;
    axis xy; 
    hold on;

    % crosses grid
    [X_cross, Y_cross] = meshgrid(1:size_x, 1:size_y);
    plot(X_cross(:), Y_cross(:), 'k+', 'MarkerSize', 6, 'LineWidth', 1);

    % labels
    xlabel(options.X_label);
    ylabel(options.Y_label);
    xlim([0.5, size_x + 0.5]);
    ylim([0.5, size_y + 0.5]);
    set(gca, 'XTick', 1:size_x, 'YTick', 1:size_y);
    set(gca, 'XTickLabel', options.X_tick_labels);
    set(gca, 'YTickLabel', options.Y_tick_labels);

    if options.withCb
        cb = colorbar();
        colormap(cmap);
        caxis([min(values), max(values)]);
        set(cb, 'Ticks', [min(values), max(values)]);
        set(cb, 'TickLabels', {'-', '+'});
        ylabel(cb, sprintf('Normalized mean\nspike count'));
        cb.Label.Rotation = -90; % Rotate the label text clockwise
    end

    title(options.title);

    grid on;
    hold off;

end
