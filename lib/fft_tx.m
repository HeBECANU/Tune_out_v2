function[out]=fft_tx(t,x,pad)
%a fft function that can handle unevenly spaced data
%uses the spread in the sampling times to dynamicaly change the resampling
%ratio
%https://dsp.stackexchange.com/questions/7788/setup-frequency-array-properly
%also can zero padd the data to decrease freq bin size
%current normalization set so that amplitude is unchanged by padding, this will mean that the power is wrong


%inputs
%   t       time vector [1,n]
%   x       signal vector [1,n]
%   padd    factor to pad by

%version 2
%fixed that t,x had to be row and col vectors
%supports zero padding data for increased freq resolution

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
if pad<1
    error('pad cant be less than zero')
elseif pad~=1
    x=[x,zeros(1,round(len_before_pad*(pad-1)))];
end

len = numel(x);             % Length of signal
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
