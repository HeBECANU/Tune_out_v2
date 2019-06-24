close all;
freq = 50;
amp = 1;
noise_amp = 0;

%Create a ~500Hz signal that runs for a seconds
sig_len =1;
t = linspace(0,sig_len,1e3); %sample rate really low until I get workflow
signal = amp*sin(freq*t);
% noise = noise_amp*randn(size(t)); %Gaussian noise

figure();
plot(t,signal)

samp_freq = 150;
samp_dur = 0.2;
samp_t0 = 0.1;
samp_times = [samp_t0:1/samp_freq:samp_t0+samp_dur];
samples = signal(samp_times);
plot(samples, 'ko')
