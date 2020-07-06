
%% Here I wish to develop the fit function for a damped anharmonic oscillator
% we define our potential as a taylor series about the origin where u(0)=0,u'(0)=0
% U =taylor_series(x,[0,0,u2,u3,u4,u5...],0)
% for a conservative potential this is pretty easy to solve
% F=- dU(x)/dx = m d^2x/d t^2
% but now we also wish to add damping to the system
% F=- dU(x)/dx -damping_c* dx/dt = m d^2x/dt^2

% so solve this system we do a finite step method in both velocity and in acceleration
% input is v,x
% dx = v
% dv = (-dU(x)/dx -damping_c*v )/m
% then the ode solver will update x,v

% now lets make a nice closed form of dU/dx
%U =taylor_series(x,[0,0,u2,u3,u4,u5...],0)
% to take the taylor series the taylor series just shifts down by 1
%dU/dx=taylor_series(x,[0,u2,u3,u4,u5...],0)
%dU/dx=deriv_taylor_series(pos_vel(1),trap_derivs,0,1)

% clear all





%% Set up the system
pos_vel_start = [0 10e-3];
trap_freq_hz = 430;
trap_freq=trap_freq_hz*2*pi;
t_lims=[0,1.3];
mass=1;
damping_time=1;
% end user var

%trap_derivs=[0,0,mass*trap_freq^2,mass*trap_freq^2,mass*trap_freq^2];
% we can calculate the spring constant as
% omega_0=sqrt(k/m) and k=d^2 u / d x^2
trap_derivs=[0,0,mass*trap_freq^2,0,1e12*mass*trap_freq^2]; %starting from the zeroth derivative
damping_ratio=1/(trap_freq*damping_time);
damping_coef=damping_ratio*2*sqrt(mass*trap_derivs(3));

%sqrt(du2(x)/m)
%sqrt(du2(0)+du3*x+du4*x^2/2)

trap_de(pos_vel_start,mass,damping_coef,trap_derivs)

t_samp = linspace(0,t_lims(2),range(t_lims)*trap_freq*100);
%system_param = [1 trap_freq^2 1 0 0];%[m, spring, damping, cubic, quartic] 

ode_opts = odeset('RelTol',1e-9,'AbsTol',1e-10,'Stats','on','InitialStep',trap_freq/1000,'MaxStep',trap_freq/200,'Vectorized',1);

%ode45
%ode113
[t,y] = ode113(@(t,X) trap_de(X,mass,damping_coef,trap_derivs),t_samp,pos_vel_start,ode_opts);
%[t,y_apx] = ode45(@(t,X) trap_DE(t,X,P0),T_win,X0);
txvdata = [t,y(:,1),y(:,2)];
%%
stfig('generated oscillaton');
clf
subplot(2,1,1)
rmsx=rms(txvdata(:,2));
plot(txvdata(:,1),txvdata(:,2)/rmsx,'k') %xval
hold on
rmsv=rms(txvdata(:,3));
plot(txvdata(:,1),txvdata(:,3)/rmsv,'r') %vel
hold off
legend(sprintf('x/%f',rmsx),sprintf('v/%f',rmsv))

subplot(2,1,2)
out_fft_x=fft_tx(txvdata(:,1),txvdata(:,2),'padding',10,'f_lim',[0,trap_freq_hz*5],'window','chebyshev','win_param',{130});
amp_fft_x=abs(out_fft_x(2,:));
plot(out_fft_x(1,:),amp_fft_x./rmsx,'k')


hold on
out_fft_v=fft_tx(txvdata(:,1),txvdata(:,3),'padding',10,'f_lim',[0,trap_freq_hz*5],'window','chebyshev','win_param',{130});
amp_fft_v=abs(out_fft_v(2,:));
plot(out_fft_v(1,:),amp_fft_v./rmsv,'r')
hold off
legend('x','v')

%%
ax = gca;
ax.YScale = 'log';


%%
stfig('oscillaton spectrogram');
fs=1/range(txvdata([1,2],1));
%spectrogram(txvdata(:,3),10000,[],linspace(0,trap_freq_hz*5,1e3),fs,'yaxis')
%spectrogram(txvdata(:,3),gausswin(10000,2),[],linspace(0,trap_freq_hz*5,1e3),fs,'yaxis')

pspectrum(txvdata(:,3),fs,'spectrogram','FrequencyResolution',10, ...
    'OverlapPercent',90,'MinTHreshold',-100,'Reassign',true,'FrequencyLimits',[trap_freq_hz/5, trap_freq_hz*10])

% subsampled_data=txvdata(1:100:end,:);
% [amp_fsst,freq_fsst,time_fsst]=fsst(subsampled_data(:,3),seconds(range(subsampled_data([1,2],1))),kaiser(1000,20));
%%
ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';
yl=ax.YLim;
ylim([trap_freq_hz/2, trap_freq_hz*10]);
axis tight

%%
[sp,fp,tp] =pspectrum(txvdata(:,3),fs,'spectrogram','Leakage',0.2,'OverlapPercent',95, ...
   'FrequencyResolution',10,'FrequencyLimits',[0, trap_freq_hz*10]);
mesh(tp,fp,log10(sp))
view(-15,60)
xlabel('Time (s)')
ylabel('Frequency (Hz)')
colormap(viridis)
%pspectrum(txvdata(:,3),fs,'spectrogram','Leakage',1,'OverlapPercent',0, ...
%    'MinThreshold',-10,'TimeResolution', 100e-3,'FrequencyLimits',[0, trap_freq_hz*10])


%% now lets try some fitting
% this approach works perfectly when there is no anharmonicity and gives agreement to num error

 fit_tlims=[0.01,0.1];
fit_mask=fast_sorted_mask(txvdata(:,1),fit_tlims(1),fit_tlims(2));
sub_polling_factor=round(1/(range(txvdata(1:2,1))*500*20));

t_subset=txvdata(fit_mask(1):sub_polling_factor:fit_mask(2),1);
x_subset=txvdata(fit_mask(1):sub_polling_factor:fit_mask(2),3);

dom_opt=[];
dom_opt.num_components=1;
dom_out=dominant_freq_components(t_subset,x_subset,dom_opt);
fit_freq=dom_out.freq(1);
fit_phase=dom_out.phase(1);
fit_amp=dom_out.amp(1);


cof_names={'amp','freq','phase','offset','damp','grad'};

 gf_opt=[];
gf_opt.domain=[[0.1,10]*rmsv;...   %amp
               [200,600];...       %freq
               [-1,1]*2*pi;...        %phase
               [-50,50]*1e-3;... %offset
               [-1e-2,10];...    %damp
               [-10,10]*1e-3;... % gradient
               ];        
gf_opt.start=[fit_amp, fit_freq,fit_phase,0,1,0];
gf_opt.rmse_thresh=rmsv*1e-4;%2e-6;
gf_opt.plot=true;
gf_opt.verbose=5;
gf_opt.level=2;
gf_out=global_fit(t_subset,x_subset,...
            @damped_sine_wave,gf_opt);

%%
 fit_tlims=[0.01,0.1];
fit_mask=fast_sorted_mask(txvdata(:,1),fit_tlims(1),fit_tlims(2));      
        
beta0  =gf_out.params;

fit_opt = statset('MaxIter',1e3,...,
            'TolFun',1e-5,...%1e4
            'TolX',1e-5,...
            'UseParallel',1);
fit_obj_decay_sine = fitnlm(txvdata(fit_mask(1):fit_mask(2),1),txvdata(fit_mask(1):fit_mask(2),3),...
    @damped_sine_wave,beta0,'CoefficientNames',cof_names,'options',fit_opt)

fit_vel=predict(fit_obj_decay_sine,txvdata(fit_mask(1):fit_mask(2),1),'Alpha',1-erf(1/sqrt(2)));
stfig('fit to osc')
clf
subplot(2,1,1)
plot(txvdata(fit_mask(1):fit_mask(2),1),txvdata(fit_mask(1):fit_mask(2),3),'r')
hold on
plot(txvdata(fit_mask(1):fit_mask(2),1),fit_vel,'k')
%plot the inital guess
plot(txvdata(fit_mask(1):fit_mask(2),1),damped_sine_wave(beta0,txvdata(fit_mask(1):fit_mask(2),1) ),'b')
legend('data','fit soln','initial guess')
hold off
subplot(2,1,2)
plot(txvdata(fit_mask(1):fit_mask(2),1),txvdata(fit_mask(1):fit_mask(2),3)-fit_vel,'k')
ylabel('residuals (vel)')


% ok so this seems to fail as time gets bigger

%% TRY AGAIN WITH MORE COMPLEX MODEL
fit_tlims=[0,1];
fit_mask=fast_sorted_mask(txvdata(:,1),fit_tlims(1),fit_tlims(2));
sub_polling_factor=round(1/(range(txvdata(1:2,1))*500*20));

t_subset=txvdata(fit_mask(1):sub_polling_factor:fit_mask(2),1);
x_subset=txvdata(fit_mask(1):sub_polling_factor:fit_mask(2),3);

dom_opt=[];
dom_opt.num_components=1;
dom_out=dominant_freq_components(t_subset,x_subset,dom_opt);
fit_freq=dom_out.freq(1);
fit_phase=dom_out.phase(1);
fit_amp=dom_out.amp(1);
            
fit_model=@(b,x) amp_freq_cpl_sine(b,x)       ;    
cof_names={'amp','freq','phase','offset','damp','grad','afc1','afc2'};
% damp_sine = @(b,x) exp(-x(:,1).*max(0,b(7))).*b(1).*...
%             sin(   taylor_series(x(:,1),[b(2),b(9),b(10)])  .*x(:,1)*pi*2   +   b(3)*pi*2)...
%             +b(4)+b(5)*x(:,2)+b(8)*x(:,1)+b(6)*x(:,3)  ;

 gf_opt=[];
gf_opt.domain=[[0.1,10]*rmsv;...   %amp
               [200,600];...       %freq
               [-1,1]*2*pi;...        %phase
               [-50,50]*1e-3;... %offset
               [-1e-2,10];...    %damp
               [-1,1]*1e3;... % gradient
               [-1,1]*1e-3;...
               [-1,1]*1e-6;...
               ];        
           %trap_freq_hz
gf_opt.start=[fit_amp, fit_freq,fit_phase,0,1,0,0.0003,-0.00000006];
gf_opt.rmse_thresh=rmsv*1e-2;%2e-6;
gf_opt.plot=true;
gf_opt.verbose=5;
gf_opt.level=1;
gf_out=global_fit(txvdata(fit_mask(1):fit_mask(2),1),txvdata(fit_mask(1):fit_mask(2),3),...
            fit_model,gf_opt);

%%
fit_model(gf_out.params,col_vec(linspace(0,2,1e3)))
        
        
%% 
 fit_tlims=[0,1];
fit_mask=fast_sorted_mask(txvdata(:,1),fit_tlims(1),fit_tlims(2));      
        
beta0  =gf_out.params;

fit_opt = statset('MaxIter',1e5,... %1e4
            'TolFun',1e-7,...
            'TolX',1e-10,...
            'UseParallel',1);
fit_obj = fitnlm(txvdata(fit_mask(1):fit_mask(2),1),txvdata(fit_mask(1):fit_mask(2),3),...
    fit_model,beta0,'CoefficientNames',cof_names,'options',fit_opt)

fit_vel=predict(fit_obj,txvdata(fit_mask(1):fit_mask(2),1),'Alpha',1-erf(1/sqrt(2)));
stfig('fit to osc')
subplot(2,1,1)
plot(txvdata(fit_mask(1):fit_mask(2),1),txvdata(fit_mask(1):fit_mask(2),3),'r')
hold on
plot(txvdata(fit_mask(1):fit_mask(2),1),fit_vel,'k')
hold off
subplot(2,1,2)
plot(txvdata(fit_mask(1):fit_mask(2),1),txvdata(fit_mask(1):fit_mask(2),3)-fit_vel,'k')
ylabel('residuals (vel)')


%ok so we can see that that fit works better but not perfectly, lets just bite the bulltet and try to fit to the
%numerical osc simulation directly to the data

%% Test the fit fun
% for harmonic functions it works great

tsamp=linspace(0,0.05,1e3);
u34=[1e6,1e3];
params=[1,420,0,0,2,u34];
vel_out=ode_fit_model(params,tsamp);
clf
plot(tsamp,vel_out,'b')
hold on
params=[1,420,pi/2,0,2,u34];
vel_out=ode_fit_model(params,tsamp);
plot(tsamp,vel_out,'k')
params=[1,420,pi/4,0,2,u34];
vel_out=ode_fit_model(params,tsamp);
plot(tsamp,vel_out,'r')
hold off

%% Phase dependence
% but for anharmonic potentials there is an amplitude phase codepenedednce
tsamp=linspace(0.01,0.05,1e3);
u34=[1e6,1e8];
params=[1,420,0,0,2,u34];
vel_out=ode_fit_model(params,tsamp);
clf
plot(tsamp,vel_out,'b')
hold on
params=[1,420,pi,0,2,u34];
vel_out=ode_fit_model(params,tsamp);
plot(tsamp,vel_out,'k')
params=[1,420,pi/4,0,2,u34];
vel_out=ode_fit_model(params,tsamp);
plot(tsamp,vel_out,'r')
hold off

%% timing
tic
a=ode_fit_model(beta0,txvdata(:,1) );
toc

%% fit with this numerical model (making von neumann cry)

 fit_tlims=[0.01,0.1];
fit_mask=fast_sorted_mask(txvdata(:,1),fit_tlims(1),fit_tlims(2));      
x_dat_subset=txvdata(fit_mask(1):fit_mask(2),3);
t_dat_subset=txvdata(fit_mask(1):fit_mask(2),1);
beta0  =cat(1,fit_obj_decay_sine.Coefficients.Estimate(1:5),0,trap_derivs(5)*0.5);

fit_opt = statset('MaxIter',0,... %1e4
            'TolFun',1e-2,...
            'TolX',1e-1,...
            'UseParallel',1,...
            'DerivStep',1e-2*col_vec([200,500,3,1,1,1e3,1e3])');
        
cof_names={'amp','freq','phase','offset','damp','trap_d3','trap_d4'};

fit_obj = fitnlm(t_dat_subset,x_dat_subset,...
    @ode_fit_model,beta0,'CoefficientNames',cof_names,'options',fit_opt)

fit_vel=predict(fit_obj,txvdata(fit_mask(1):fit_mask(2),1),'Alpha',1-erf(1/sqrt(2)));
stfig('fit to osc')
clf
subplot(2,1,1)
plot(txvdata(fit_mask(1):fit_mask(2),1),txvdata(fit_mask(1):fit_mask(2),3),'r')
hold on
plot(txvdata(fit_mask(1):fit_mask(2),1),fit_vel,'k')
%plot the inital guess
plot(txvdata(fit_mask(1):fit_mask(2),1),ode_fit_model(beta0,txvdata(fit_mask(1):fit_mask(2),1) ),'b')
legend('data','fit soln','initial guess')
hold off
subplot(2,1,2)
plot(txvdata(fit_mask(1):fit_mask(2),1),txvdata(fit_mask(1):fit_mask(2),3)-fit_vel,'k')
ylabel('residuals (vel)')

%% 
predictor=t_dat_subset;
response=x_dat_subset;
cost_fun=@(x) sqrt(sum(abs(ode_fit_model(col_vec(x),predictor)-response).^2)/(numel(response)-numel(x)));

params_matrix=repmat(beta0',50,1);
params_matrix(:,7)=linspace(trap_derivs(5)/100,trap_derivs(5),size(params_matrix,1));
cost_out=col_row_fun_mat(cost_fun,params_matrix,2);
plot(params_matrix(:,7)/trap_derivs(5),cost_out)

%%
fit_tlims=[0.01,0.01+4/fit_freq];
fit_mask=fast_sorted_mask(txvdata(:,1),fit_tlims(1),fit_tlims(2));
%sub_polling_factor=round(1/(range(txvdata(1:2,1))*500*20));
sub_polling_factor=1;

t_subset=txvdata(fit_mask(1):sub_polling_factor:fit_mask(2),1);
x_subset=txvdata(fit_mask(1):sub_polling_factor:fit_mask(2),3);

dom_opt=[];
dom_opt.num_components=1;
dom_out=dominant_freq_components(t_subset,x_subset,dom_opt);
fit_freq=dom_out.freq(1);
fit_phase=dom_out.phase(1);
fit_amp=dom_out.amp(1);
            
%cof_names={'amp','freq','phase','offset','damp','trap_d3','trap_d4'};

 gf_opt=[];
gf_opt.domain=[[0.1,10]*rmsv;...   %amp
               [200,700];...       %freq
               [-1,1]*2*pi;...        %phase
               [-50,50]*1e-3;... %offset
               [-1e-2,10];...    %damp
               [-1,1]*1e3;...
               [-0.05,10]*trap_derivs(5);...
               ];        
           %trap_freq_hz
gf_opt.start=cat(1,fit_amp,fit_freq-50,fit_phase,0,2,0,trap_derivs(5)*0.5);
gf_opt.rmse_thresh=rmsv*1e-2;%2e-6;
gf_opt.plot=true;
gf_opt.verbose=5;
gf_opt.level=2;
gf_out=global_fit(t_subset,x_subset,...
            @ode_fit_model,gf_opt);
        


%% fit the ode numeical model to all the data

 fit_tlims=[0.01,1];
fit_mask=fast_sorted_mask(txvdata(:,1),fit_tlims(1),fit_tlims(2));      
        
beta0  =fit_obj.Coefficients.Estimate;

fit_opt = statset('MaxIter',0,... %1e4
            'TolFun',1e-7,...
            'TolX',1e-10,...
            'UseParallel',1,...
            'DerivStep',1e-2*col_vec([200,500,3,1,1,1e3,1e3])');
        
cof_names={'amp','freq','phase','offset','damp','trap_d3','trap_d4'};

fit_obj = fitnlm(txvdata(fit_mask(1):fit_mask(2),1),txvdata(fit_mask(1):fit_mask(2),3),...
    @ode_fit_model,beta0,'CoefficientNames',cof_names,'options',fit_opt)

fit_vel=predict(fit_obj,txvdata(fit_mask(1):fit_mask(2),1),'Alpha',1-erf(1/sqrt(2)));
stfig('fit to osc')
clf
subplot(2,1,1)
plot(txvdata(fit_mask(1):fit_mask(2),1),txvdata(fit_mask(1):fit_mask(2),3),'r')
hold on
plot(txvdata(fit_mask(1):fit_mask(2),1),fit_vel,'k')
%plot the inital guess
plot(txvdata(fit_mask(1):fit_mask(2),1),ode_fit_model(beta0,txvdata(fit_mask(1):fit_mask(2),1) ),'b')
legend('data','fit soln','initial guess')
hold off
subplot(2,1,2)
plot(txvdata(fit_mask(1):fit_mask(2),1),txvdata(fit_mask(1):fit_mask(2),3)-fit_vel,'k')
ylabel('residuals (vel)')

%%

function vel_out=ode_fit_model(p,t_in)
%cof_names={'amp','freq','phase','offset','damp','trap_d3','trap_d4'};

fprintf('%s\n',sprintf('%.3f,',p))
trap_freq_hz = p(2);
trap_freq=trap_freq_hz*2*pi;
%% phase shift
% starting the oscillation at a given phase
% for a harmonic trap this works great, but for an anharomic one phase is not that well defined (would have to
% numericaly int. to find it)
% ke_start/m=v^2+ w^2 x^2
%pos_vel_start = [sin(p(3)),trap_freq*cos(p(3))]*p(1); %amplitude is in units of velocity

%% the other way is to apply a phase delay
pos_vel_start=[0,1]*p(1); 
t_in=t_in+(p(3)+pi/2)/trap_freq;


p=col_vec(p);
t_in=col_vec(t_in);
if t_in(1)<0
    error('t must be zero or positive')
end
is_first_elm_zero=t_in(1)==0;
if ~is_first_elm_zero
    t_in=cat(1,0,t_in);
end

mass=1;
damping_time=p(5);
% end user var

%trap_derivs=[0,0,mass*trap_freq^2,mass*trap_freq^2,mass*trap_freq^2];
% we can calculate the spring constant as
% omega_0=sqrt(k/m) and k=d^2 u / d x^2
trap_derivs=[0,0,mass*trap_freq^2,p(6),p(7)]; %starting from the zeroth derivative
damping_ratio=1/(trap_freq*damping_time);
damping_coef=damping_ratio*2*sqrt(mass*trap_derivs(3));
% set up the ode solver
ode_opts = odeset('AbsTol',0.01*1e-3,'Stats','off','InitialStep',trap_freq/1000,'MaxStep',trap_freq/20,'Vectorized',1);
%'RelTol',1e-9,
%ode45   
%ode113
[t,y] = ode113(@(t,X) trap_de(X,mass,damping_coef,trap_derivs),t_in,pos_vel_start,ode_opts);
vel_out=y(:,2)+p(4);
vel_out=col_vec(vel_out);
if ~is_first_elm_zero
    vel_out=vel_out(2:end);
end

end

        
function y = damped_sine_fundemental(p,t)
    y = exp(-abs(p(3))*t).*sin(2*pi*p(1)*t-p(2))+p(4);
end

function y = damped_sine_wave(p,t)
    %cof_names={'amp','freq','phase','offset','damp','grad'};
    undamped_angular_freq=2*pi*p(2);
    damping_time=max([0,p(5)]); %prevent a growing oscillation
    damping_ratio=1/(undamped_angular_freq*damping_time);
    damped_angular_freq=undamped_angular_freq*sqrt(1-damping_ratio^2);
    if damping_time>0
        damp_term=exp(-t./damping_time);
    else
        damp_term=1;
    end
    amp=p(1).*damp_term;
    y = amp.*sin(damped_angular_freq*t-p(3))+p(4)+p(6).*t;
end

function dx = trap_de(pos_vel,mass,damping,trap_derivs)
    dx(1) = pos_vel(2); 
    dx(2) = (-deriv_taylor_series(pos_vel(1),trap_derivs,0,1) - damping*pos_vel(2))./mass;
    dx = dx';
end %function


% function dx = trap_de(vel_acc,system_param)
%     dx(1) = vel_acc(2)- system_param(3)*vel_acc(1); %comment out after the first term to turn off damping
%     dx(2) = (-system_param(2)*vel_acc(1) - (system_param(4)*vel_acc(1)^2)/2 - (system_param(5)*vel_acc(1)^3)/6)/system_param(1);
%     dx = dx';
% end %function


% function y = mod_sine(p,t)
%     y = exp(-p(3)*t).*sin(p(1).*t-p(2))+p(4);
% end

% function dx = trap_de(vel_acc,system_param)
%     dx(1) = vel_acc(2)- system_param(3)*vel_acc(1); %comment out after the first term to turn off damping
%     dx(2) = (-system_param(2)*vel_acc(1) - (system_param(4)*vel_acc(1)^2)/2 - (system_param(5)*vel_acc(1)^3)/6)/system_param(1);
%     dx = dx';
% end %function

% function f = objfun(P,data,IC)
%     [t,y_apx] = ode45(@(t,X) trap_de(t,X,P),data(:,1),IC);
%     f = sum((data(:,2)-y_apx(:,1)).^2); %L2 error
% end

function P = down_sample(func, P,T)
    P = func(P,T);
end


function out=amp_freq_cpl_sine(b,x)           
%cof_names={'amp','freq','phase','offset','damp','grad','afc1','afc2'};

amplitude=exp(-x.*max(0,b(5))).*b(1);
frequency=taylor_series(amplitude,[b(2),b(7),b(8)]);
out=amplitude.*sin(frequency.*x*pi*2 +b(3))+b(4)+b(6).*x;

end
