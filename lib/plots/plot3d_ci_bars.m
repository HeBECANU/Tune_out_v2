function [h]=plot3d_ci_bars(x, y, z, ci_neg_x,ci_pos_x, ci_neg_y,ci_pos_y, ci_neg_z,ci_pos_z,color, varargin)
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
        if ~isempty(ci_neg_x)
            if isempty(ci_pos_x)
                ci_pos_x=-ci_neg_x;
            end
            xMin = x(i) + ci_neg_x(i);
            xMax = x(i) + ci_pos_x(i);
            xB = [xMin, xMax];
            h=plot3(xB, yV, zV, '-k');
            set(h, 'LineWidth', 2);
            set(h, 'color', color);
        end
        if ~isempty(ci_neg_y)
            if isempty(ci_pos_y)
                ci_pos_y=-ci_neg_y;
            end
            yMin = y(i) + ci_neg_y(i);
            yMax = y(i) + ci_pos_y(i);
            yB = [yMin, yMax];
            h=plot3(xV, yB, zV, '-k');
            set(h, 'LineWidth', 2);
             set(h, 'color', color);
        end
        if ~isempty(ci_neg_z)
            if isempty(ci_pos_z)
                ci_pos_z=-ci_neg_z;
            end
            zMin = z(i) + ci_neg_z(i);
            zMax = z(i) + ci_pos_z(i);
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