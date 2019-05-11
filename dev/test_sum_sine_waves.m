%test_sum_sine_waves(chirped)

tdat=linspace(0,0.1,1e4)';
sine_params=[[1,50,0];[0.5,150,0]];
unchriped=sum_sine_waves(tdat,sine_params,0);
sine_params=[[1,1,0];[0.5,3,0]];
chriped=sum_sine_waves_chirped(tdat,sine_params,0,[50,1]);
plot(tdat,unchriped)
hold on
plot(tdat,chriped)
hold off



%%

tdat=linspace(0,1,1e4)';
sine_params=[[1,50,0];[0.5,150,0]];
unchriped=sum_sine_waves(tdat,sine_params,0);
sine_params=[[1,0,1,0];[3,0,0.5,0]];
chriped=sum_sine_waves_amp_freq_chirped(tdat,sine_params,0,[50,0]);
plot(tdat,unchriped)
hold on
plot(tdat,chriped)
hold off

