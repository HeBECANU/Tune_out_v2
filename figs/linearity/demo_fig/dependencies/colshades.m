function [c_light, c_dark] = colshades(col,k_lumin)
% Creates lighter and darker shade of given RGB-triplet by scaling
% luminescence
%
%   col:        N-by-3 RBG array
%   k_lumin:    scaling factor for luminescence (default is ---)
%
%   c_light:    lighter shades
%   c_dark:     darker shades
%   
%
% DKS
% 2018-08-17

if ~exist('k_lumin','var')
    k_lumin=1.7;        % default L scaling factor
end

c_hsl=colorspace('RGB->HSL',col);   % transform color to HSL space

% create light/dark shades and represent in RGB
c_light=colorspace('HSL->RGB',[1,1,k_lumin].*c_hsl);
c_dark = colorspace('HSL->RGB',[1,1,1/k_lumin].*c_hsl);

end