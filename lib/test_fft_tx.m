% test script for fft_tx
times=linspace(0,1e3,1e6);
amp_strong=10;
freq_strong=100;
amp_weak=1;
freq_weak=133;
noise_amp=0;
testfun=@(t) amp_strong*sin(2*pi*t*freq_strong+pi)+amp_weak*sin(2*pi*t*freq_weak+pi)+noise_amp*rand(size(t));
val=testfun(times);
out=fft_tx(times,val,'padding',10,'window','gauss','win_param',{5});
figure(5)
plot(out(1,:),abs(out(2,:)))
[max_amp,nearest_idx]=max(abs(out(2,:)));
max_freq=out(1,nearest_idx);
logic_str = {'FAIL', 'pass'};
amp_tol=1e-2; %fractional
freq_tol=1e-2; %absolute
fprintf('INFO: peak freq error      : %g\n',max_freq-freq_strong)
fprintf('INFO: peak amp  frac error : %g\n',max_amp/amp_strong-1)
fprintf('TEST: peak freq within tolerance: %s\n',logic_str{1+(abs(max_freq-freq_strong)<freq_tol)})
fprintf('TEST: peak amp  within tolerance: %s\n',logic_str{1+(abs(max_amp/amp_strong-1)<amp_tol)})



%%
times=linspace(0,1e3,1e6);
testfun=@(t) 100*sin(2*pi*t*100+pi)+1*sin(2*pi*t*133+pi);
val=testfun(times);
out=fft_tx(times,val,'padding',10,'window','gauss','win_param',{5});
figure(5)
plot(out(1,:),abs(out(2,:)))