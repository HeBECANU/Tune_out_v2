
function out=sum_sine_waves(x,components,offset)
% components is a Nx3 array of amp,freq,phase
% offset shifts eveything
% example usage
%sine_params=[[1,2,3];[4,5,6]];
%sum_sine_waves(linspace(0,1,1e3),sine_params,0)
%component_amp=rowfun(@(amp,freq,phase) amp.*sin(x.*2.*pi*freq+phase),components);
x=x(:);
component_amp=col_row_fun_mat(@(params) params(1).*sin(x.*2.*pi*params(2)+params(3)),components,2);
out=sum(component_amp,1)+offset;
out=out(:);
end