time=linspace(0,5,1e5);
sine_params=[[1,1,0];[0.1,3,0]];
ydat=sum_sine_waves_chriped(time,sine_params,0,[0,50]);
plot(time,ydat)