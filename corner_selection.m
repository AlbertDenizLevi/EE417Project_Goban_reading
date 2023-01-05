% This code processes a series of images of a Go board, or "goban," and
% identifies the locations of the black and white stones on the board. It
% also displays the identified stone locations on the original image. The
% code consists of a main loop that processes each image in a directory,
% and several functions that perform different tasks within the loop.
% 
% The main loop begins by calling the dir function to get a list of all the
% files in the current directory with a .png file extension. It then loops
% through the list of files, and for each file, it does the following:
% 
% Loads the image using the imread function. 
% Displays the image using the imshow function. 
% Adds a button to the figure using the uicontrol
% function, which allows the user to use previously saved "ginput" points,
% which are points clicked on the image by the user. Tries to load
% previously saved ginput points from a file using the load function. If
% the ginput points were successfully loaded, the figure with the button is
% closed. If the ginput points could not be loaded, the code prompts the
% user to click on the four corner points of the goban in the image using
% the ginput function. The clicked points are then saved to a file. Calls
% the improve_corner_points function to improve the accuracy of the corner
% points. This function takes as input the image and the corner points, and
% returns improved corner points as output. Calls the compute_grid_points
% function to compute the grid points of the goban. This function takes as
% input the improved corner points, and returns the grid points as output.
% Calls the compute_brightness function to compute the brightness values at
% the grid points. This function takes as input the image and the grid
% points, and returns the brightness values as output. Calls the
% display_results function to display the results. This function takes as
% input the image, the original and improved corner points, the grid
% points, and the brightness values, and displays the results by plotting
% the points on the image. The improve_corner_points function improves the
% accuracy of the corner points by defining a window around each corner
% point and searching for the point within the window with the highest
% gradient magnitude. It returns the improved corner points as output.
% 
% The compute_grid_points function computes the grid points of the goban by
% creating an ideal grid and transforming it to match the improved corner
% points using the fitgeotrans and transformPointsForward functions. It
% returns the grid points as output.
% 
% The compute_brightness function computes the brightness values at the
% grid points by interpolating the image values at the grid points using
% the interp2 function. It returns the brightness values as output.
% 
% The display_results function displays the results by plotting the
% original and improved corner points, the grid points, and the black and
% white stone points on the image. It also saves the gamestate, or the
% state of the goban with the black and white stone points, to a .mat file.
% 
% The use_saved_ginputs function is a callback function for the button in
% the main loop. It loads the saved ginput points from a file and closes
% the figure with the button.
% 
% The display_fake_goban function displays a fake goban with filled circles
% as the stones to show the saved gamestate. It takes the gamestate as
% input and displays




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
    
end




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


