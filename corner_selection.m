
% Close all figures, clear the command window, and clear all variables
close all;
clear;
clc;
files = dir('*.png');

% Loop through the list of files
for i = 1:length(files)
    % Load the image
    image = imread(files(i).name);

    % Display the image
    figure
    imshow(image)


    % Add a button to the figure
    button = uicontrol('Style', 'pushbutton', 'String', 'Use saved ginputs', ...
        'Position', [20 20 150 20], 'Callback', @use_saved_ginputs);

    % Try to load the saved ginputs from a file
    try
        load(['saved_ginputs_' files(i).name(1:end-4) '.mat'], 'x', 'y');
        % If the ginputs were successfully loaded, close the figure with the button
        close(gcf)
    catch
        % If the ginputs could not be loaded, get them from the user
        [x, y] = ginput(4);
        close(gcf)
        % Save the ginputs to a file
        save(['saved_ginputs_' files(i).name(1:end-4) '.mat'], 'x', 'y');
    end

    % Improve the corner points
    [x_improved, y_improved] = improve_corner_points(image, x, y);

    % Compute the grid points
    [x_grid, y_grid] = compute_grid_points(x_improved, y_improved);

    % Compute the brightness values
    brightness = compute_brightness(image, x_grid, y_grid);

    % Display the results
    display_results(image, x, y, x_improved, y_improved, x_grid, y_grid, brightness, files);
    
%     % Load the gamestate from the .mat file
%     load(['gamestate_' files(i).name(1:end-4) '.mat'], 'gamestate');
%     
%     % Display the fake goban
%     display_fake_goban(gamestate);

end

% function display_fake_goban(gamestate)
%     % Set the stone colors
%     black_color = [0 0 0];
%     white_color = [1 1 1];
%     empty_color = [0.5 0.5 0.5];
% 
%     % Initialize the image
%     image = zeros(19*20, 19*20, 3);
%     
%     % Loop through the grid points
%     for i = 1:19
%         for j = 1:19
%             % Determine the color of the current grid point
%             if gamestate(i, j) == 2
%                 color = black_color;
%             elseif gamestate(i, j) == 1
%                 color = white_color;
%             else
%                 color = empty_color;
%             end
% 
%             % Fill the current grid point with the appropriate color
%             color_image = repmat(color, [20 20]);
%             image((i-1)*20+1:i*20, (j-1)*20+1:j*20, :) = color_image;
%         end
%     end
% 
%     % Display the image
%     imshow(image)
% end



% Callback function for the button
function use_saved_ginputs(src, event)
% Load the saved ginputs from a file
load(['saved_ginputs_' files(i).name(1:end-4) '.mat'], 'x', 'y');

% Close the figure with the button
close(gcf)
end

function [x_improved, y_improved] = improve_corner_points(image, x, y)
% Define a window size around each corner point
window_size = round(min(size(image, 1), size(image, 2)) / 20);

% Define a grid of points in the window
[X, Y] = meshgrid(-window_size:window_size, -window_size:window_size);

% Initialize improved corner points
x_improved = zeros(size(x));
y_improved = zeros(size(y));

% Loop through the four corners
for i = 1:4
    % Extract the window around the current corner point
    y_start = max(1, round(y(i))-window_size);
    y_end = min(size(image, 1), round(y(i))+window_size);
    x_start = max(1, round(x(i))-window_size);
    x_end = min(size(image, 2), round(x(i))+window_size);
    window = image(y_start:y_end, x_start:x_end);

    % Detect corners in the window
    corner_points = corner(window, 'MinimumEigenvalue');

    % If there are any corner points detected, select the closest one to the original corner point
    if ~isempty(corner_points)
        [~, index] = min(sum((corner_points - [window_size+1, window_size+1]).^2, 2));
        x_improved(i) = corner_points(index, 1) + x_start - 1;
        y_improved(i) = corner_points(index, 2) + y_start - 1;
        % Otherwise, use the original corner point
    else
        x_improved(i) = x(i);
        y_improved(i) = y(i);
    end
end
end

function [x_grid, y_grid] = compute_grid_points(x_improved, y_improved)
% Compute the grid points
x_grid = linspace(min(x_improved), max(x_improved), 19);
y_grid = linspace(min(y_improved), max(y_improved), 19);
end



function brightness = compute_brightness(image, x_grid, y_grid)
% Compute the average brightness of each grid point
brightness = zeros(19, 19);
window_size = 2;
for i = 1:19
    for j = 1:19
        y_start = max(1, round(y_grid(i))-window_size);
        y_end = min(size(image, 1), round(y_grid(i))+window_size);
        x_start = max(1, round(x_grid(j))-window_size);
        x_end = min(size(image, 2), round(x_grid(j))+window_size);
        window = image(y_start:y_end, x_start:x_end);
        brightness(i, j) = mean(window(:));
    end
end
end

function [black_points, white_points, thresholds] = classify_stones(brightness)
% Set the lower threshold for black stones to be 150
black_threshold = 150;

% Compute the histogram of the brightness values
[counts, bins] = hist(brightness(:), 100);

% Find the right-most bin with a count of 0
rightmost_zero_bin = find(counts == 0, 1, 'last');

% Set the upper threshold for white stones to be the right-most bin with a count of 0
white_threshold = bins(rightmost_zero_bin);

% Set the thresholds for black and white stones
thresholds = [black_threshold, white_threshold];

% Classify the stones on the go board as black or white based on the thresholds
black_points = brightness < black_threshold;
white_points = brightness > white_threshold;
end




function display_histogram(brightness)
% Display a histogram of the brightness values
figure
histogram(brightness(:), 'BinWidth', 1.5)
hold on
end

function display_results(image, x, y, x_improved, y_improved, X_grid, Y_grid, brightness, files)
% Display all grid points and the black and white stone grid points on the image
[black_points, white_points] = classify_stones(brightness);


figure
imshow(image)
hold on
for i = 1:19
    for j = 1:19
        % Determine the color of the current grid point
        if black_points(i,j)
            color = 'go';
        elseif white_points(i,j)
            color = 'mo';
        else
            color = 'k+';
        end


        % Plot the current grid point
        plot(X_grid(j), Y_grid(i), color, 'MarkerSize', 10)
    end

    % Initialize the gamestate matrix with all zeros
    gamestate = zeros(19, 19);
    
    % Set the elements of gamestate that correspond to black stones to 2
    gamestate(black_points) = 2;
    
    % Set the elements of gamestate that correspond to white stones to 1
    gamestate(white_points) = 1;

    % Save the gamestate matrix to a file
    save(['gamestate_' files(i).name(1:end-4) '.mat'], 'gamestate');
end

end


