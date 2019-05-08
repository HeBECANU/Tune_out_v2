function [st_pts,level_xing]=stpt_and_level_xing(tdat,xdat,smoothing_time,level)
%TODO
% build test script

%change input into col vectors
tdat=tdat(:);
xdat=xdat(:);

if size(tdat,2)~=1 || size(xdat,2)~=1 || ~isequal(size(xdat),size(tdat))
    error('input the wrong size')
end
    
if ~isnan(smoothing_time) && smoothing_time~=0 && ~isempty(smoothing_time)
    xdat_filt=gaussfilt(tdat,xdat,smoothing_time);
else
    xdat_filt=xdat;
end


xdat_mean=mean(xdat_filt);

xderiv=gradient(xdat_filt,tdat); %calculate the single sided derivative
%find the zero crossings of the derivative (st pts)
[xing_idx,xing_t]=crossing(xderiv,tdat,0,[]);
xing_idx=xing_idx(:);
xing_t=xing_t(:);
st_pts.idx=xing_idx;
st_pts.time=xing_t;

%determine if these crossings are above or below the mean value with the assumption that if the zero derivative pt
%is above the mean then the voltage is turning arround to be a negative slope
st_pts.above_mean=xdat(xing_idx)>xdat_mean;

% find the curvature of these crossings
x2deriv=gradient(xderiv,tdat);
%determine if these crossings they are increasing or decreasing in derivative
st_pts.positive_curvature=x2deriv(xing_idx)>0;

%optionaly compute the level crossing
if nargout>1
    if isnan(level) || isempty(level)
        xing_level=xdat_mean;
    end
    [xing_idx,xing_t]=crossing(xdat,tdat,xing_level,[]);
    xing_idx=xing_idx(:);
    xing_t=xing_t(:);
    level_xing.idx=xing_idx;
    level_xing.time=xing_t;
    level_xing.positive_deriv=xderiv(xing_idx)>0;
end



end