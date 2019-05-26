%demonstrate and test two different smoothing approaches based on gaussian convolution



smooth_time=0.001;
xdat=col_vec(linspace(0,1,1e5));
ydat_clean=0*cos(2*pi*10.*xdat+2*pi*1)...
            +10*(xdat>0.5);
[~,idx]=closest_value(xdat,0.2);
ydat_clean(idx)=ydat_clean(idx)+100;
ydat_noisy=ydat_clean+...
    rand(numel(xdat),1)*1e-3+...
    randn(numel(xdat),1)*1e-3;

sr=1/range(xdat([1,2]));

tic
kernel_sd=20;
klen=ceil(kernel_sd*smooth_time*sr);
alpha =(klen-1)/(smooth_time*sr*2);
kernel = gausswin(klen,alpha);
kernel=kernel/sum(kernel);
ydat_smooth1= col_vec(nanconv(ydat_clean,kernel,'edge','1d'));
toc
tic
ydat_smooth2=gaussfilt(xdat,ydat_noisy,smooth_time);
toc
stfig('smooth compare');
clf
plot(xdat,ydat_smooth1,'b')
hold on
plot(xdat,ydat_smooth2,'r')

%ydat_smooth1 max (0.2,0.7998) fwhm 2*0.00058 gives sigma at 0.0027

%ydat_smooth2 max (0.2,0.39935) fwhm 2*0.00118 gives sigma at 0.0028



%%  check that the impulse response is the right width

smooth_time=0.01;
xdat=col_vec(linspace(0,1,1e5));
ydat_clean=0*xdat;
[~,idx]=closest_value(xdat,0.5);
ydat_clean(idx)=ydat_clean(idx)+100;

sr=1/range(xdat([1,2]));

tic
kernel_sd=20;
klen=ceil(kernel_sd*smooth_time*sr);
alpha =(klen-1)/(smooth_time*sr*2);
kernel = gausswin(klen,alpha);
kernel=kernel/sum(kernel);
ydat_smooth1= col_vec(nanconv(ydat_clean,kernel,'edge','1d'));
sum_ydat_smooth1=sum(ydat_smooth1);

ydat_smooth2=col_vec(gaussfilt(xdat,ydat_clean,smooth_time));
sum_ydat_smooth2=sum(ydat_smooth2);

(sum_ydat_smooth2-sum_ydat_smooth1)./mean([sum_ydat_smooth1,sum_ydat_smooth2])

modelfun = @(b,x) b(3)*gauss1d(b(1),b(2),x); %simple linear model
opts = statset('TolFun',1e-10,'TolX',1e-10,...
            'MaxIter',1e4,... %1e4
            'UseParallel',1);
beta0 = [0.5,smooth_time*1.5,1]; %intial guesses
fit_mdl = fitnlm(xdat,ydat_smooth1,modelfun,beta0,'Options',opts)
meth1_fit_data=predict(fit_mdl,xdat);


fit_mdl = fitnlm(xdat,ydat_smooth2,modelfun,beta0,'Options',opts)
meth2_fit_data=predict(fit_mdl,xdat);


stfig('smooth compare');
clf
subplot(2,1,1)
plot(xdat,ydat_smooth1,'b')
hold on
plot(xdat,meth1_fit_data,'g')
plot(xdat,ydat_smooth2,'r')
plot(xdat,meth2_fit_data,'m')
subplot(2,1,2)
plot(xdat,ydat_smooth1-ydat_smooth2,'k')


