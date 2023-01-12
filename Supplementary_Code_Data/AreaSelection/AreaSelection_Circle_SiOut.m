% clearvars -except R_Crop_Bothnm Center x_crop_Bothnm y_crop_Bothnm
%%
function [maskedImage, Center, radius, xunit, yunit] = AreaSelection_Circle_SiOut(originalImage)

% originalImage = R_Crop_Bothnm(:,:,1);

% imagesc(originalImage)
[rows cols ~] = size(originalImage)
imshow(originalImage)
axis('on', 'image');
% Maximize the window to make it easier to draw.
g = gcf;
g.WindowState = 'maximized';

uiwait(helpdlg('Select ROI enclosing inner circle.'));

h.Radius = 0;
while h.Radius == 0
	h = drawcircle('Color','k','FaceAlpha',0);
	if h.Radius == 0
		uiwait(helpdlg('You double-clicked.  You need to single click, then drag, then single click again.'));
	end
end
% Get coordinates of the circle.,'DrawingArea',[x,y,w,h]
angles = linspace(0, 2*pi, 10000);
x = cos(angles) * h.Radius + h.Center(1);
y = sin(angles) * h.Radius + h.Center(2);
%%
C_x = ( max(x) + min(x) ) / 2;
C_y = ( max(y) + min(y) ) / 2;
Center = [C_x C_y]

%% Demo to have the user click and draw a circle over an image, then blacken outside the circle and crop out the circular portion into a new image.
% clc;% Clear the command window.
% fprintf('Beginning to run %s.m ...\n', mfilename);
% close all;  % Close all figures (except those of imtool.)
% imtool close all;  % Close all imtool figures.
% workspace;  % Make sure the workspace panel is showing
% fontSize = 15;
% 
% %% Get image.
% [rows, columns, ~] = size(originalImage);
% figure(1)
% imshow(originalImage);
% axis('on', 'image');
% %% Maximize the window to make it easier to draw.
% g = gcf;
% g.WindowState = 'maximized';
% 
% uiwait(helpdlg('Select outer circle of ROI for only Silicon.'));
% 
% h.Radius = 0;
% while h.Radius == 0
% 	h = drawcircle('Color','k','FaceAlpha',0);
% 	if h.Radius == 0
% 		uiwait(helpdlg('You double-clicked.  You need to single click, then drag, then single click again.'));
% 	end
% end
% 
% %% Get coordinates of the circle.,'DrawingArea',[x,y,w,h]
% angles = linspace(0, 2*pi, 10000);
% x = cos(angles) * h.Radius + h.Center(1);
% y = sin(angles) * h.Radius + h.Center(2);
% 
% %% Get a mask of the circle
% 
% % if x > y
% %    x = y
% % elseif y > x
% %    y = x
% % end
% 
% mask = poly2mask(x, y, rows, columns);
% 
% %% Mask the image with the circle.
% maskedImage = bsxfun(@times, originalImage, cast(mask, class(originalImage)));
% 
% 
% %% Crop the image to the bounding box.
% props = regionprops(mask, 'BoundingBox');
% maskedImage = imcrop(maskedImage, props.BoundingBox);
% % 
% end
% % 
% % figure
% % plot(x,y,'r')
% % xlim([min(x) max(x)])
% % ylim([min(y) max(y)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clearvars -except R_Crop_Bothnm Center cropped_im_R_0nm_out
% 
% originalImage = R_Crop_Bothnm(:,:,1);
% otherIm = cropped_im_R_0nm_out; 

%% Demo to have the user click and draw a circle over an image, then blacken outside the circle and crop out the circular portion into a new image.
clc;% Clear the command window.
fprintf('Beginning to run %s.m ...\n', mfilename);
close all;  % Close all figures (except those of imtool.)
imtool close all;  % Close all imtool figures.
workspace;  % Make sure the workspace panel is showing
fontSize = 15;

%% Get image.
[rows, columns, ~] = size(originalImage);

fig = figure(1)
% imshow(originalImage);
% axis('on', 'image');
%% Maximize the window to make it easier to draw.

axis ij;
%
axis manual;
%

%
[m n] = size(originalImage);  %%d1 is the image which I want to display
%
axis([0 m 0 n])

imshow(originalImage);
axis on
uiwait(helpdlg('Click on anywhere to draw a circle for outer silicon ROI.'));

% Enlarge figure to full screen.

set(gcf, 'units','normalized','outerposition',[0 0 1 1]);

coordinates_input = ginput(1)

rowCenter= round(coordinates_input(2));
columnCenter = round(coordinates_input(1));
fprintf('You clicked on coordinates (row, column) = (%f, %f)\nWhich is the pixel in row %d, column %d\n', ...
coordinates_input(1), coordinates_input(2), rowCenter, columnCenter);
radius = sqrt((abs(coordinates_input(1) - Center(1)))^2 + (abs(coordinates_input(2) - Center(2))^2))
croppi = viscircles(Center,radius)
pause(1)
uiwait(helpdlg('Select OK if you confirm circle.'));
pause(1)
close(fig)

%% function h = circle(x,y,r)
% hold on
th = linspace(0, 2*pi, 10000);
xunit = radius * cos(th) + Center(1);
yunit = radius * sin(th) + Center(2);
% h = plot(xunit, yunit);
% hold off

%% 


% g = gcf;
% g.WindowState = 'maximized';
% 
% set(gcf,'Position',[0 0 columns rows]); %fullscreen
% uiwait(helpdlg('Click on point .'));
% waitforbuttonpress 
% point1 = get(gcf,'CurrentPoint') 
%%
% uiwait(helpdlg('Select inner circle of ROI for only Silicon .'));
% 
% h.Radius = 0;
% while h.Radius == 0
% 	h = drawcircle('Color','k','FaceAlpha',0);
% 	if h.Radius == 0
% 		uiwait(helpdlg('You double-clicked.  You need to single click, then drag, then single click again.'));
% 	end
% 
% end
%% Get coordinates of the circle.,'DrawingArea',[x,y,w,h]
% angles = linspace(0, 2*pi, 10000);
% x = cos(angles) * h.Radius + h.Center(1);
% y = sin(angles) * h.Radius + h.Center(2);

% if x > y
%    x = y;
% elseif y > x
%    y = x;
% end

% mask = uint16(poly2mask(xunit, yunit, m, n))
% 
% maskedImage = bsxfun(@times, originalImage, cast(mask, class(originalImage)));


%% Crop the image to the bounding box.
% props = regionprops(mask, 'BoundingBox');
% maskedImage = imcrop(maskedImage, props.BoundingBox);

%% Get a mask of the circle
mask = uint16(poly2mask(xunit, yunit, rows, columns));

%% Mask the image with the circle.
%% 	maskedImage = bsxfun(@times, originalImage, cast(mask, class(originalImage)));
maskedImage = originalImage.*mask;

% Crop the image to the bounding box.
% props = regionprops(mask, 'BoundingBox');
% maskedImage = imcrop(maskedImage2, props.BoundingBox);

end

