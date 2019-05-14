function out=sum_sine_waves_amp_freq_chirped(x,components,offset,freq_series)
% components is a (Nx(2+J)] array of harmonic,phase,amp_tay_series[J terms]
% offset shifts eveything
% example usage
%sine_params=[[1,2,3];[4,5,6]];
%sum_sine_waves(linspace(0,1,1e3),sine_params,0)
%component_amp=rowfun(@(amp,freq,phase) amp.*sin(x.*2.*pi*freq+phase),components);

if size(components,2)<3
    error('components wrong size need at least harmonic,phase,amp_d0')
end
x=x(:);
sine_fun=@(params) taylor_series(x,params(3:end),0).*sin(x.*taylor_series(x,freq_series,0).*2.*pi*params(1)+params(2));
component_amp=col_row_fun_mat(sine_fun,components,2);
out=sum(component_amp,1)+offset;
out=out(:);
end

