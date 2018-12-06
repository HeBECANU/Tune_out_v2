function[out]=fft_tx(t,x,varargin)
% fft_tx - general purpose fft
% built to hande most use cases: uneven data, windowing, padding to increase the freq sampling. 
% all built to keep the returned amplitudes correct, this will mean that the power is wrong.
% https://dsp.stackexchange.com/questions/7788/setup-frequency-array-properly
% uses the spread in the sampling times to dynamicaly change the resampling
% ratio

% Syntax:   freq_amp=fft_tx(times,val)
%           freq_amp=fft_tx(times,val,'padding',10,'window','gauss','win_param',{5});
%
% Inputs:
%   t          - vector, sample times (will deal with unmatched vector dimension)
%   x          - vector, the signal   (will deal with unmatched vector dimension)
%   Optional Name Value Pairs:
%       'padding'        - float, value>=1, factor to padd the data by to increase the frequeny sampling (but not resolvability)
%       'window'         - string, windowing function to apply to the data before processing.
%                       ['none','hamming','gauss','blackman','hanning']
%       win_param       - cell array, extra inputs after the length argument to the windowing function
%                           - gauss ,win_param={3}, reciprocal of the standard deviation
%                           - all others (besides none) 'symmetric'(recomended) or 'periodic'

% Outputs:
%    freq_amp - 2*L matrix
%               - freq_amp(1,:) are the frequency bins
%               - freq_amp(2,:) are the (complex) amplitudes
%
% Example: 
%         times=linspace(0,1e3,1e6);
%         testfun=@(t) 100*sin(2*pi*t*100+pi)+1*sin(2*pi*t*133+pi);
%         val=testfun(times);
%         out=fft_tx(times,val,'padding',10,'window','gauss','win_param',{5});
%         figure(5)
%         plot(out(1,:),abs(out(2,:)))

% Other m-files required: none
% Also See: test_fft_tx
% Subfunctions: none
% MAT-files required: none
%
% Known BUGS/ Possible Improvements
%    - chosable resampling stratagey
%
% Author: Bryce Henson
% email: Bryce.Henson@live.com
% Last revision:2018-11-21

%------------- BEGIN CODE --------------

%adaptively deal with the data if its in row or col fromat
if size(size(t),2)==2
    if size(t,1)~=1 && size(t,2)==1
        t=t';
    elseif size(t,1)<=1 && size(t,2)<=1
        error('thats not a vector in t')
    end
else
    error('you have tried to input the wrong shape in t')
end
if size(size(x),2)==2
    if size(x,1)~=1 && size(x,2)==1
        x=x';
    elseif size(x,1)==1 && size(x,2)==1
        error('thats not a vector in x')
    end
else
    error('you have tried to input the wrong shape in x')
end  


p = inputParser;
addParameter(p,'window','none',@(x) sum(contains({'none','hamming','gauss','blackman','hanning'},x))==1);
addParameter(p,'win_param',{},@(x) true);
addParameter(p,'padding','',@(x) isfloat(x) && x>=1);
parse(p,varargin{:});
pad=p.Results.padding;
window_fun=p.Results.window;
window_param=p.Results.win_param;

%first sort the data
[t,i]=sort(t);
x=x(i);
%find the number of sampling times present
sample_times=uniquetol(t(1:end-1)-t(2:end),1e-10); %avoid machine error
if size(sample_times,1)>1 %if not uniform
    fprintf('resampling\n');
    dyn_res_factor=10+500*abs(std(sample_times)/mean(sample_times));
    t_resample=linspace(min(t),max(t),size(t,1)*dyn_res_factor);
    x=interp1(t,x,t_resample,'spline')';
    t=t_resample;
end
dt = (t(2)-t(1));             % Sampling period
%disp(num2str(dt))
fs=1/dt;
len_before_pad = numel(x);             % Length of signal

%apply windowing function
switch window_fun
    case 'hamming'
        win=hamming(len_before_pad,window_param{:})';
    case 'gauss'
        win=gausswin(len_before_pad,window_param{:})'; 
    case 'blackman'
        win=blackman(len_before_pad,window_param{:})'; 
    case 'hanning'
        win=hann(len_before_pad,window_param{:})'; 
    case 'none'
        %do nothing
end
if isequal(window_fun,'none')
    win=win.*(len_before_pad/sum(win));
    x=x.*win;
end



if pad<1
    error('pad cant be less than one')
elseif pad~=1
    x=[x,zeros(1,round(len_before_pad*(pad-1)))];
end

len = numel(x);             % Length of signal
fs=fs;%*len_before_pad/len;
%disp(num2str(L))
%t = (0:L-1)*T;        % Time vector
y = fft(x);
amp = y/len_before_pad;
if  mod(len,2) %odd case
    niq=ceil(len/2);
    amp = amp(1:niq);
    amp(2:end) = 2*amp(2:end);
    f = fs/len*((0:niq-1));
else
    niq=len/2+1;
    amp = amp(1:niq);
    amp(2:end) = 2*amp(2:end); %multiply all non DC values by 2 because of half cut
    f = fs*(0:(len/2))/len;  
end


out=[f;amp];


end
