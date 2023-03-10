function [maskedImage] = AreaSelection_Circle_Mod_SiIn(originalImage,Center)

% Demo to have the user click and draw a circle over an image, then blacken outside the circle and crop out the circular portion into a new image.
clc;% Clear the command window.
fprintf('Beginning to run %s.m ...\n', mfilename);
close all;  % Close all figures (except those of imtool.)
imtool close all;  % Close all imtool figures.
workspace;  % Make sure the workspace panel is showing
fontSize = 15;

% Get image.
[rows, columns, ~] = size(originalImage);

% Get a mask of the circle
mask = uint16(poly2mask(xunit, yunit, rows, columns));

% Mask the image.
maskedImage = originalImage.*mask;


% Crop the image to the bounding box.
% props = regionprops(mask, 'BoundingBox');
% maskedImage = imcrop(maskedImage, props.BoundingBox);

fprintf('Done running %s.m ...\n', mfilename);

end