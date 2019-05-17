
function out=sum_sine_waves_chirped(x,components,offset,freq_series)
% components is a Nx3 array of amp,harmonic,phase
% offset shifts eveything
% example usage
%sine_params=[[1,2,3];[4,5,6]];
%sum_sine_waves(linspace(0,1,1e3),sine_params,0)
%component_amp=rowfun(@(amp,freq,phase) amp.*sin(x.*2.*pi*freq+phase),components);

if size(components,2)~=3
    error('components wrong size')
end
x=x(:);
sine_fun=@(params) params(1).*sin(x.*taylor_series(x,freq_series,0).*2.*pi*params(2)+params(3));
component_amp=row_col_fun_mat(sine_fun,components,1);
out=sum(component_amp,1)+offset;
out=out(:);
end

