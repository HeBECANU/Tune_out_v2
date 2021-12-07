function [h]=plot3d_errorbars(x, y, z, ex, ey, ez,color, varargin)
%from https://stackoverflow.com/questions/23654296/multi-dimensional-2d-better-3d-scatter-plot-with-different-errorbars-in-matlab
% by dan https://stackoverflow.com/users/1011724/dan

% create the standard 3d scatterplot
%hold off;
%h=plot3(x, y, z, varargin{:});

% looks better with large points
%set(h, 'MarkerSize', 25);



% now draw the vertical errorbar for each point
for i=1:length(x)
        xV = [x(i); x(i)];
        yV = [y(i); y(i)];
        zV = [z(i); z(i)];
        if ~isempty(ex)
            xMin = x(i) + ex(i);
            xMax = x(i) - ex(i);
            xB = [xMin, xMax];
            h=plot3(xB, yV, zV, '-k');
            set(h, 'LineWidth', 2);
            set(h, 'color', color);
        end
        if ~isempty(ey)
            yMin = y(i) + ey(i);
            yMax = y(i) - ey(i);
            yB = [yMin, yMax];
            h=plot3(xV, yB, zV, '-k');
            set(h, 'LineWidth', 2);
             set(h, 'color', color);
        end
        if ~isempty(ez)
            zMin = z(i) + ez(i);
            zMax = z(i) - ez(i);
            zB = [zMin, zMax];
            % draw error bars
            if range(zB)~=0
                h=plot3(xV, yV, zB, '-k');
                set(h, 'LineWidth', 2);
                set(h, 'color', color);
            end
        end
end


end