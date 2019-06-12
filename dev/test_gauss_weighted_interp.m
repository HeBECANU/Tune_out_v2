%

sigma=0.05;

x=col_vec(rand(1,10));
y=taylor_series(x,[0,4,50],0.5);
y=y+randn(size(y,1),size(y,2));

stfig('test gauss interp')
clf
subplot(3,1,1)
plot(x,y,'x')
hold on
%x_samp=linspace(min(x),max(x),1e3);
x_samp=linspace(-1,2,1e3);
y_smooth=gauss_weighted_interp(x,y,x_samp,sigma);
plot(x_samp,y_smooth,'-k')
y_smooth=exp_weighted_interp(x,y,x_samp,sigma);
plot(x_samp,y_smooth,'-r')

subplot(3,1,2)
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


subplot(3,1,3)
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

