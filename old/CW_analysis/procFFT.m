function procFFT(amplitude, t)
	persistent fig_fft_handle	
	%t=linspace(0, 5, 1000);	amplitude=10*sin(2*pi*10*t)+1*rand(size(t)); %test	
	%Function that computes the power spectrum of data
	%Inputs
	% amplitude - array of amplitude data
	% t - array of time data
	% Nfft- length of data that the FFT is averaged over (higher number
	%           gives a higher resolution
	% plt - plots the power spectrum
	% fstart - start freq. for plot
	% fstop - stop freq. for plot

	%Outputs
	%P-Amplitude From FFT (sqrt(Hz))
	%f-Freq from FFT (Hz)
	%len-Length of FFT
	Nfft=2^nextpow2(length(t));
	plt=1; 
	if ishandle(fig_fft_handle)		
		figure(fig_fft_handle);
		x_lims=xlim;
		y_lims=ylim;
	else
		fig_fft_handle=figure;
	end
	fstart=0; fstop=1/(t(2)-t(1));
	[P,f,len]=power_spec(amplitude,t,Nfft,plt,fstart,fstop);
	save('C:\Program Files\MATLAB\R2011b\Pow.mat', 'P', 'f')
	[pks,locs] = findpeaks(P, 'MINPEAKHEIGHT',max(P)/3, 'THRESHOLD', max(P)/100);
	hold on;
	if exist('x_lims', 'var')
		xlim(x_lims);
		ylim(y_lims);		
	end	
	plot(f(locs),pks,'k^')
	hold off;
	title('FFT(TOF) hello')
	