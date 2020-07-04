function [c, c_light, c_dark] = palette(n_colors,bg,L)
% Create a pretty palette of distinguishable colors for graphics
%   n_colors: number of colors to generate
%   bg: background. e.g. {'w','k'}, 1x3 RGB; (default 'w'hite)
%   L: lightness multiplier for accessories (default 1.7)
%   
%   c: n_colors X 3 array
%   c_light: n_colors X 3 array
%   c_dark: n_colors X 3 array
%
%
% DKS
% 2018-06-04

% parse inputs
if ~exist('bg','var')
    bg='w';
end
if ~exist('L','var')
    L=1.7;
end

% create colors
c = distinguishable_colors(n_colors,bg);    % main colors
[c_light,c_dark] = colshades(c,L);          % light/dark shades
    
end