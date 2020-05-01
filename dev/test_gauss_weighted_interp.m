%

sigma=0.05;

x=col_vec(rand(1,10));
x=sort(x);
y_noisefree=taylor_series(x,[0,4,50],0.5);
y=y_noisefree+randn(size(y_noisefree,1),size(y_noisefree,2));

stfig('test gauss interp 1')
clf
plot(x,y,'x')
hold on
%x_samp=linspace(min(x),max(x),1e3);
x_samp=linspace(-0.1,1.1,1e3);
y_smooth=gauss_weighted_interp(x,y,x_samp,sigma);
plot(x_samp,y_smooth,'-k')
y_smooth=exp_weighted_interp(x,y,x_samp,sigma);
plot(x_samp,y_smooth,'-r')
plot(x,gaussfilt(x,y,sigma),'-')
plot(x,y_noisefree,'-')
plot(x,gaussfilt(x,y_noisefree,sigma),'-')
legend('data','gauss interp','exp interp','gauss filt','noiseless','filt noiseless')

stfig('test gauss interp 2')
clf
y=x*0;
y(round(numel(y)/2))=1;
plot(x,y,'x')
hold on
%x_samp=linspace(min(x),max(x),1e3);
x_samp=linspace(-1,2,1e3);

y_smooth=gauss_weighted_interp(x,y,x_samp,sigma);
plot(x_samp,y_smooth,'-k')
y_smooth=exp_weighted_interp(x,y,x_samp,sigma);
plot(x_samp,y_smooth,'-r')


stfig('test gauss interp 2')
clf
x=linspace(0,1,30);
y=x*0;
y=cat(2,y,1);
x=cat(2,x,1.1);
plot(x,y,'x')
hold on
%x_samp=linspace(min(x),max(x),1e3);
x_samp=linspace(-1,2,1e3);
y_smooth=gauss_weighted_interp(x,y,x_samp,sigma);
plot(x_samp,y_smooth,'-k')
y_smooth=exp_weighted_interp(x,y,x_samp,sigma);
plot(x_samp,y_smooth,'-r')


%%

x=[0,1,2,3,4,5];
y=[1,0,0,0,0,0];
xq=[-10,-1,0,3,10];
gauss_weighted_interp(x,y,xq,1)

