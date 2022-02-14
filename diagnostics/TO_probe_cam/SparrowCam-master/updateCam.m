function updateCam(~,event,hImage)

global hImageAxes hXSlice hYSlice hLineSliceX hLineSliceY hCrosshairX hCrosshairY;
global edit_BRan_x edit_BRan_y edit_FWHM_x edit_FWHM_y edit_maxpos_x edit_maxpos_y edit_maxval edit_fps;
global settings imageRes pixsize_x pixsize_y margin prev_toc data;
global cmap background backgroundData setBackgroundData;

% This callback function updates the displayed frame and the histogram.

margin = str2double(settings.margin);
% smoothing of 1D slices is CPU costy?
toggle_smoothing = str2double(settings.toggle_smoothing);

toggle_pixel_scaling = 0;

% subtracting the background
if background
    if setBackgroundData
        backgroundData = event.Data;
        setBackgroundData = 0;
    end
    data = event.Data - backgroundData;
else
    data = event.Data;
end

% Display the current image frame. 
set(hImage, 'CData', ind2rgb(data, cmap));

mod_data = medfilt2(data);

% Slice at the center of imageRes(1).
[val,ind] = max(mod_data(:));

[Y,X] = ind2sub(size(data),ind); % Position of maximum(inversed).

maxint = 255; % Upper limit of intensity.
slice_x = data(Y,:); % Slice on x direction. 
slice_y = data(:,X); % Slice on y direction.
thr_2 = 1/2; % Threshold of peak(1/2).
thr_e2 = 1/exp(2); % Threshold of peak(1/(e^2)).


if toggle_pixel_scaling
    ccd_res_x = 1920; % Actual resolution for x direction(DMK 23UX174).
    ccd_res_y = 1200; % Actual resolution for y direction(DMK 23UX174).
    c_x = pixsize_x*(ccd_res_x/imageRes(2)); % Length for 1 pixel on x direction.
    c_y = pixsize_y*(ccd_res_y/imageRes(1)); % Length for 1 pixel on y direction.
else
    c_x = pixsize_x;
    c_y = pixsize_y;
end

% specifiing coordinates for image corners (must coincide with axes limits)
set(hImage, 'XData', [0 imageRes(2)*c_x]);
set(hImage, 'YData', [0 imageRes(1)*c_y]);

set(hImageAxes, 'Visible', 'on');

set(hImageAxes, 'XLimMode', 'manual');
set(hImageAxes, 'YLimMode', 'manual');

% adjusting the limits of axes for the image
hImageAxes.XLim = [0 imageRes(2)*c_x];
hImageAxes.YLim = [0 imageRes(1)*c_y];
hImageAxes.TickDir = 'in';

% Keeps ticks in the picture
set(hImageAxes, 'XTickMode', 'auto');
set(hImageAxes, 'YTickMode', 'auto');

% deactivating Tick labels
hImageAxes.YTickLabel = '';
hImageAxes.XTickLabel = '';

set(hImageAxes, 'XGrid', settings.xgrid);
set(hImageAxes, 'YGrid', settings.ygrid);

colormap(cmap);
try
 colorbar(hImageAxes);
catch
    'Sorry colorbar is not available in this system.'
    
end
% Draw crosshair at the maximum
hCrosshairX.XData = [0 imageRes(2)*c_x];
hCrosshairX.YData = [Y*c_y Y*c_y];

hCrosshairY.XData = [X*c_x X*c_x];
hCrosshairY.YData = [0 imageRes(1)*c_y];

% Calculating widths of peak
half_x = slice_x > (val*thr_2); % Judgement if it is included in peak or not on x direction(1/2).
tot_half_x = sum(half_x); % Number of pixels on x direction.
area_x = tot_half_x*c_x; % Total area on x direction.

half_y = slice_y > (val*thr_2); % Judgement if it is included in peak or not on y direction(1/2).
tot_half_y = sum(half_y); % Number of pixels on x direction.
area_y = tot_half_y*c_y; % Total area on x direction.

two_e_x = slice_x > (val*thr_e2); % Judgement if it is included in peak or not on x direction(1/(e^2)).
tot_e_x = sum(two_e_x); % Number of pixels on x direction.
range_x = tot_e_x*c_x; % Total area on x direction.

two_e_y = slice_y > (val*thr_e2); % Judgement if it is included in peak or not on y direction(1/(e^2)).
tot_e_y = sum(two_e_y); % Number of pixels on x direction.
range_y = tot_e_y*c_y; % Total area on x direction.

% mask for smoothing 1D slices
conv_num = str2double(settings.mask_size); 
conv_c = (1/conv_num)*ones(1,conv_num);

function xplot()

    grid_x = linspace(0, imageRes(2)*c_x, imageRes(2)); %Grid for x direction.
    
    set(hImageAxes.Parent, 'CurrentAxes', hXSlice);
    hXSlice.Position = [hImageAxes.Position(1) hXSlice.Position(2) hImageAxes.Position(3) hXSlice.Position(4)];
    hLineSliceX.XData = grid_x;
    hLineSliceX.YData = slice_x;
%     title('Position x and its intensity');
    xlabel('x [\mum]');
%     ylabel('Intensity');
    
    axis([0,imageRes(2)*c_x,0,maxint]);
    
    if toggle_smoothing
        hold(hXSlice,'on');

        smooth_x = conv(double(slice_x),double(conv_c),'same');
        plot(grid_x,smooth_x);
        hold(hXSlice,'off');
    end

end

% Slice at the center of imageRes(2).

function yplot()
    
    grid_y = linspace(0, imageRes(1)*c_y, imageRes(1)); %Grid for y direction.
    set(hImageAxes.Parent, 'CurrentAxes', hYSlice);
    hYSlice.Position = [ hYSlice.Position(1) hImageAxes.Position(2) hYSlice.Position(3) hImageAxes.Position(4)];
    hLineSliceY.XData = slice_y;
    hLineSliceY.YData = grid_y;
    hYSlice.YDir = 'reverse';

    ylabel('y [\mum]');

    axis([0,maxint,0,imageRes(1)*c_y]);
    if toggle_smoothing
        hold(hYSlice,'on');

        smooth_y = conv(double(slice_y),double(conv_c),'same');
        plot(smooth_y,grid_y);
        hold(hYSlice,'off');
    end

end

% updating values in monitors
set(edit_maxval, 'String', sprintf('%3d',val));
set(edit_maxpos_x, 'String', sprintf('%4.0f um',X*c_x));
set(edit_maxpos_y, 'String', sprintf('%4.0f um',Y*c_y));
set(edit_FWHM_x, 'String', sprintf('%4.0f um',area_x));
set(edit_FWHM_y, 'String', sprintf('%4.0f um',area_y));
set(edit_BRan_x, 'String', sprintf('%4.0f um',range_x));
set(edit_BRan_y, 'String', sprintf('%4.0f um',range_y));

xplot();
yplot();

% FPS counting:

time = toc;
fps = 1./(time - prev_toc);
prev_toc = time;

set(edit_fps, 'String', sprintf('%2.1f',fps));


% Refresh the display.
drawnow


end


% function [inds, vals] = center_of_mass(X, Y, data)
% 
%     t1 = data.*X;
%     t2 = data.*Y;
% 
%     mass = sum(data(:));
%     x_pos = sum(t1(:))/mass;
%     y_pos = sum(t2(:))/mass;
%     vals = [x_pos y_pos];
%     
%     [~,x_ind] = min(abs(X(1,:) - x_pos));
%     [~,y_ind] = min(abs(Y(:,1) - y_pos));
%     inds = [x_ind, y_ind];
%     
% end
% 
% [inds, vals] = center_of_mass(X, Y, data);