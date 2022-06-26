
addpath('../lib/Core_BEC_Analysis/lib/') %add the path to set_up_project_path
set_up_project_path('..')

hebec_constants
close('all')

%%
trap_afreq=426*2*pi;
beam_waist=20e-6;
gauss_sigma=beam_waist/2;

%gauss_amp=-3e-30;


beam_power=0.1;%0.120;
peak_intensity= 2*beam_power/(pi*(beam_waist^2));
detuning_hz=6e12;
to_polz_si_deriv=1.8017e-53;
dipole_pot=- (1/(2*const.epsilon0*const.c))*to_polz_si_deriv(1)*peak_intensity*detuning_hz;
gauss_amp=dipole_pot

%%
m=const.mhe;
pot_fun=@(x,trap_afreq,gamp,gsig) (1/2) *m *trap_afreq^2 * x.^2 + gaussian_function_1d(x,gsig,0,gamp);
xrange=gauss_sigma*3;
x_samp=linspace(-xrange,xrange,1e3);
pot_samp_harm=pot_fun(x_samp,trap_afreq,0,0);
pot_samp_comb=pot_fun(x_samp,trap_afreq,gauss_amp,gauss_sigma)-pot_fun(0,trap_afreq,0,gauss_sigma);



xfactor=1e6;
yfactor=1e30;
yfactor=10^(round(log10(yfactor)));
stfig('Trap Potential')
plot(x_samp*xfactor,pot_samp_comb*yfactor)
hold on
plot(x_samp*xfactor,pot_samp_harm*yfactor)
hold off
xlabel('Position, $y$ ($\mu$m)')
ylabel(sprintf('Potential, $U(y)$ ($10^{%d}$ J)',round(log10(1/yfactor))))


x_max=10e-3/trap_afreq;
%gauss_amp=1e-20;
excitaiton_energy=pot_fun(x_max,trap_afreq,gauss_amp,gauss_sigma);
period_intagrand=@(x) 1./sqrt(excitaiton_energy-pot_fun(x,trap_afreq,gauss_amp,gauss_sigma) );
num_period=sqrt(2*m)*integral(period_intagrand,-x_max,x_max);
num_freq=1/num_period ;





%%
gauss_amp=-1e-32;
gauss_sigma=10e-6/2;


anh_dep=[];
anh_dep.osc_x_amp=linspace(gauss_sigma/100,gauss_sigma*6,1e3);
anh_dep.osc_afreq=anh_dep.osc_x_amp*nan;

gauss_afreq=gauss_trap_freq(gauss_sigma,gauss_amp,m);
ho_comb_afreq=sqrt(gauss_afreq^2+trap_afreq^2);

iimax=numel(anh_dep.osc_x_amp);
for ii=1:iimax
    x_max=anh_dep.osc_x_amp(ii);
    %gauss_amp=1e-20;
    excitaiton_energy=pot_fun(x_max,trap_afreq,gauss_amp,gauss_sigma);
    period_intagrand=@(x) 1./sqrt(excitaiton_energy-pot_fun(x,trap_afreq,gauss_amp,gauss_sigma) );
    num_period=sqrt(2*m)*integral(period_intagrand,-x_max,x_max);
    num_freq=1/num_period ;
    anh_dep.osc_afreq(ii)=num_freq*2*pi;
end


colors_main=[[[214,72,154];
[102,181,69];
[151,90,214];
[208,153,44];
[212,73,58]]]./255;
hsv=colorspace('RGB->HSV',colors_main(:,:));
hsv(:,2)=hsv(:,2);
colors_shaded=colorspace('HSV->RGB',hsv);
font_name='cmr10';
linewidth=1.5;
font_size=11;



stfig('ahn gauss')
clf
set(gca, 'FontName', font_name)
set(gca, 'FontSize', font_size)
%
subplot(1,2,1)
xfactor=1e6;
yfactor=1/(2*pi);
plot(anh_dep.osc_x_amp*xfactor,anh_dep.osc_afreq*yfactor,'LineWidth',linewidth,'Color',colors_main(1,:))
hold on
yline(ho_comb_afreq*yfactor,'--','LineWidth',linewidth*1.2,'Color',colors_main(3,:))
yline(trap_afreq*yfactor,':','LineWidth',linewidth*1.2,'Color',colors_main(4,:))
xline(gauss_sigma*2*xfactor,'-.','LineWidth',linewidth*1.2,'Color',colors_main(2,:))
xline(osc_amp_x_operation*xfactor,'--','LineWidth',linewidth*1.5,'Color',colors_main(5,:))
hold off
xlabel('Oscillation Amplitude ($\mu$m)')
ylabel('Combined Trap Frequency $\Omega_{Net}$ (Hz)')
legend('Simuation','Harm. Approx. $\Omega$','Magnetic Trap $\Omega$','Beam Waist','Exp. Osc. Amp.')
legend('Location','east')
yrange=abs(trap_afreq-ho_comb_afreq);
ylim([trap_afreq-yrange*0.05,ho_comb_afreq+yrange*0.05]*yfactor)

sp_pos=get(gca,'Position');
an=annotation('textbox', [sp_pos(1)+0.27, sp_pos(2)+sp_pos(4)*0.90, 0, 0], 'string', '(A)','FontSize',round(font_size*1.5),'FontName',font_name);
set(gca, {'XColor', 'YColor'}, {'k', 'k'});
set(gca,'linewidth', 1.0)
set(gca,'TickLength',[0.015,0])
set(gca, 'FontName', font_name)
set(gca, 'FontSize', font_size)


infered_gauss_freq=sqrt(anh_dep.osc_afreq.^2 - trap_afreq.^2);
subplot(1,2,2)
plot(anh_dep.osc_x_amp*xfactor,infered_gauss_freq*yfactor,'LineWidth',linewidth,'Color',colors_main(1,:))
hold on
yline(gauss_afreq*yfactor,'--','LineWidth',linewidth*1.2,'Color',colors_main(3,:))
xline(gauss_sigma*2*xfactor,'-.','LineWidth',linewidth*1.2,'Color',colors_main(2,:))
xline(osc_amp_x_operation*xfactor,'--','LineWidth',linewidth*1.5,'Color',colors_main(5,:))
hold off
xlabel('Oscillation Amplitude ($\mu$m)')
ylabel('Probe Trap Frequency, $\Omega_{Probe}$ (Hz)')
legend('Simuation','Harm. Approx. $\Omega$','Beam Waist','Exp. Osc. Amp.')
legend('Location','best')
ylim([0,gauss_afreq*1.05]*yfactor)

sp_pos=get(gca,'Position');
an=annotation('textbox', [sp_pos(1)+0.27, sp_pos(2)+sp_pos(4)*0.90, 0, 0], 'string', '(B)','FontSize',round(font_size*1.5),'FontName',font_name);


set(gca, {'XColor', 'YColor'}, {'k', 'k'});
set(gca,'linewidth', 1.0)
set(gca,'TickLength',[0.015,0])
set(gca, 'FontName', font_name)
set(gca, 'FontSize', font_size)


set(gcf,'Units','Pixels')
%set(gcf,'Position',[1068,355,676,453])
fig_width_px=1000;
fig_aspect_ratio=0.3; %0.67;
set(gcf,'Position',[700,355,fig_width_px,fig_width_px*fig_aspect_ratio])       
hold off

%%
fig_name='combined_trap_freq_sim';
%fig_name='pal_mean_v_acc_dyn_static';
fig_dir='./figs/thesis_figs';
export_fig(fullfile(fig_dir,strcat(fig_name,'.svg')))
pause(0.1)
export_fig(fullfile(fig_dir,strcat(fig_name,'.pdf')))


%
fprintf('Probe freq  is %.3f in zero amp limit \n',gauss_afreq/(2*pi))
fprintf('Combined freq  is %.3f in zero amp limit \n',ho_comb_afreq/(2*pi))
fprintf('Probe freq with osc amp of %f um was a factor of %.3f than zero amp limit \n',...
    gauss_sigma*2*1e6,interp1(anh_dep.osc_x_amp,infered_gauss_freq,gauss_sigma*2)/gauss_afreq)


%%
osc_amp_v_operation=12e-3;
osc_amp_x_operation=osc_amp_v_operation/trap_afreq;

stfig('oscc sens')
clf
xfactor=1e6;
yfactor=1e4/(2*pi);
font_size=12;

sens_metric=infered_gauss_freq.*anh_dep.osc_x_amp;
[max_sens_metric,idx_max_sens]=max(sens_metric);
max_sens_osc=anh_dep.osc_x_amp(idx_max_sens);
plot(anh_dep.osc_x_amp*xfactor,sens_metric*yfactor,'LineWidth',linewidth,'Color',colors_main(1,:))
hold on
plot(max_sens_osc*xfactor,max_sens_metric*yfactor,'o','MarkerSize',10,'LineWidth',2,'Color',colors_main(3,:))
harm_sens_approx=anh_dep.osc_x_amp*gauss_afreq;
plot(anh_dep.osc_x_amp*xfactor,harm_sens_approx*yfactor,'-.','LineWidth',linewidth,'Color',colors_main(4,:))
%yline(trap_afreq*yfactor,'--')
xline(gauss_sigma*2*xfactor,':','LineWidth',linewidth*1.5,'Color',colors_main(2,:))
xline(osc_amp_x_operation*xfactor,'--','LineWidth',linewidth*1.5,'Color',colors_main(5,:))
hold off
xlabel('Oscillation Amplitude ($\mu$m)')
ylabel('Sensitivty Metric , $\Omega_P \cdot A_y$ (arb. u.)')
legend('Simuation','Optimum','Harmonic Approx.','Beam Waist','Exp. Osc. Amp.')
legend('Location','best')
legend('FontSize',font_size*1)
ylim([0,max(sens_metric)*1.8]*yfactor)
set(gca, 'FontName', font_name)
set(gca, 'FontSize', font_size)
set(gca, {'XColor', 'YColor'}, {'k', 'k'});
set(gca,'linewidth', 1.0)
set(gca,'TickLength',[0.015,0])

set(gcf,'Units','Pixels')
%set(gcf,'Position',[1068,355,676,453])
fig_width_px=500;
fig_aspect_ratio=0.67; %0.67;
set(gcf,'Position',[100,355,fig_width_px,fig_width_px*fig_aspect_ratio])       
hold off


%%
fig_name='probe_freq_sens_sim';
%fig_name='pal_mean_v_acc_dyn_static';
fig_dir='./figs/thesis_figs';
export_fig(fullfile(fig_dir,strcat(fig_name,'.svg')))
export_fig(fullfile(fig_dir,strcat(fig_name,'.pdf')))


fprintf('Optimum osc amp %.3f of waist \n', max_sens_osc/(gauss_sigma*2))
fprintf('Harm approx at optimum better by factor of %.3f \n', gauss_afreq*max_sens_osc/(max_sens_metric))
fprintf('Operating with an osc amp of %.1f mm/s, %.1f um was a factor of %.3f the optimum \n',...
    osc_amp_v_operation*1e3,osc_amp_x_operation*1e6,interp1(anh_dep.osc_x_amp,sens_metric,osc_amp_x_operation)/max_sens_metric)

fprintf('Operating with an osc amp of %.1f mm/s, %.1f um was a factor of %.3f the harm approx \n',...
    osc_amp_v_operation*1e3,osc_amp_x_operation*1e6,...
    interp1(anh_dep.osc_x_amp,sens_metric,osc_amp_x_operation)/...
    interp1(anh_dep.osc_x_amp,harm_sens_approx,osc_amp_x_operation))

%%


trap_afreq=412*2*pi;
gauss_sigma=10e-6/2;
gauss_amp=-1e-30;

gauss_afreq=gauss_trap_freq(gauss_sigma,gauss_amp,m);

m=const.mhe;
pot_fun=@(x,trap_afreq,gamp,gsig) (1/2) *m *trap_afreq^2 * x.^2 + gaussian_function_1d(x,gsig,0,gamp);
xrange=gauss_sigma*3;
x_samp=linspace(-xrange,xrange,1e3);
pot_samp_harm=pot_fun(x_samp,trap_afreq,0,0);
pot_samp_comb=pot_fun(x_samp,trap_afreq,gauss_amp,gauss_sigma)-pot_fun(0,trap_afreq,0,gauss_sigma);
pot_samp_gauss=pot_samp_comb-pot_samp_harm;


stfig('Trap Potential');
set(gca, 'FontName', font_name)
set(gca, 'FontSize', font_size)
font_size=12;

xfactor=1e6;
yfactor=1e30;
yfactor=10^(round(log10(yfactor)));
plot(x_samp*xfactor,pot_samp_harm*yfactor,'--','LineWidth',linewidth,'Color',colors_main(4,:))
hold on
plot(x_samp*xfactor,pot_samp_gauss*yfactor,':','LineWidth',linewidth,'Color',colors_main(3,:))
plot(x_samp*xfactor,pot_samp_comb*yfactor,'-','LineWidth',linewidth,'Color',colors_main(1,:))
xline([-1,1]*gauss_sigma*2*xfactor,'-.','LineWidth',linewidth*1.2,'Color',colors_main(2,:))
xline([-1,1]*osc_amp_x_operation*xfactor,'--','LineWidth',linewidth*1.5,'Color',colors_main(5,:))
hold off
xlabel('Position, $y$ ($\mu$m)')
ylabel(sprintf('Potential, $U(y)$ ($10^{%d}$ J)',round(log10(1/yfactor))))
%yl=ylim;
ylim([1,-3].*min(pot_samp_comb)*1.1*yfactor)
xlim([x_samp(1),x_samp(end)].*xfactor)
legend('Magnetic Trap','Probe','Combined','Beam Waist','','Exp. Osc. Amp.')
legend('Location','north')
legend('FontSize',font_size*1)

set(gca, {'XColor', 'YColor'}, {'k', 'k'});
set(gca,'linewidth', 1.0)
set(gca,'TickLength',[0.015,0])
set(gca, 'FontName', font_name)
set(gca, 'FontSize', font_size)

set(gcf,'Units','Pixels')
%set(gcf,'Position',[1068,355,676,453])
fig_width_px=600;
fig_aspect_ratio=0.67; %0.67;
set(gcf,'Position',[100,355,fig_width_px,fig_width_px*fig_aspect_ratio])       
hold off

%%
fig_name='potential_diagram_comb_trap';
%fig_name='pal_mean_v_acc_dyn_static';
fig_dir='./figs/thesis_figs';
export_fig(fullfile(fig_dir,strcat(fig_name,'.svg')))
export_fig(fullfile(fig_dir,strcat(fig_name,'.pdf')))


fprintf('probe beam trap freq %f Hz \n', gauss_afreq/(2*pi))


%% lets look at if the linearity is any different with oscillation amplitude
% damm looks like its perfectly linear even with anharmonic osc

osc_amp_x_samp=[1,10,20].*1e-3./trap_afreq;
beam_detune_samp=linspace(-1,1,1e2)*1e9;

beam_waist=10e-6;
beam_sigma=beam_waist/2;
beam_power=10e-3;
peak_intensity= 2*beam_power/(pi*(beam_waist^2));
detuning_hz=30e6;

anh_lin=[];
anh_lin.detuning=beam_detune_samp;
anh_lin.osc_amp=osc_amp_x_samp;
anh_lin.ho_approx_afreq=nan(numel(osc_amp_x_samp),numel(beam_detune_samp));
anh_lin.num_afreq=anh_lin.ho_approx_afreq*nan;

% from tune_out_sensitivity_calc
to_polz_si_deriv_1=1.80168e-53;



iimax=numel(anh_lin.osc_amp);
jjmax=numel(anh_lin.detuning);
for ii=1:iimax
    x_max_tmp=anh_lin.osc_amp(ii);
    for jj=1:jjmax
        beam_detune_tmp=anh_lin.detuning(jj);
        dipole_pot=- (1/(2*const.epsilon0*const.c))*to_polz_si_deriv_1*peak_intensity*beam_detune_tmp;
        % find the predicted freq
        gauss_afreq=gauss_trap_freq(beam_sigma,dipole_pot,m);
        ho_comb_afreq=sqrt(gauss_afreq^2+trap_afreq^2);
        anh_lin.ho_approx_afreq(ii,jj)=ho_comb_afreq;
        %gauss_amp=1e-20;
        excitaiton_energy=pot_fun(x_max_tmp,trap_afreq,dipole_pot,beam_sigma);
        period_intagrand=@(x) 1./sqrt(excitaiton_energy-pot_fun(x,trap_afreq,dipole_pot,beam_sigma) );
        num_period=sqrt(2*m)*integral(period_intagrand,-x_max_tmp,x_max_tmp);
        num_freq=1/num_period ;
        anh_lin.num_afreq(ii,jj)=num_freq*2*pi;

    end
end

%
stfig('anharmonic linearity')
clf
hold on
yscale=1/(2*pi);
xscale=1e-9;
for ii=1:iimax
    signal_ho_tmp=trap_afreq.^2 - anh_lin.ho_approx_afreq(ii,:).^2;
    plot(anh_lin.detuning*xscale,signal_ho_tmp*yscale)
    signal_anh_tmp=trap_afreq.^2 - anh_lin.num_afreq(ii,:).^2;
    plot(anh_lin.detuning*xscale,signal_anh_tmp*yscale)
end
hold off
xlabel('detuning (GHz)')
ylabel('signal $\Omega^2-\Omega^2$ (Hz$^2$)')
box on


%%

beam_waist=11e-6;
beam_sigma=beam_waist/2;
beam_power=0.100;
peak_intensity= 2*beam_power/(pi*(beam_waist^2));
beam_detune_tmp=4e9;

dipole_pot=- (1/(2*const.epsilon0*const.c))*to_polz_si_deriv_1*peak_intensity*beam_detune_tmp;
dipole_pot/const.kb
gauss_trap_freq(beam_sigma,dipole_pot,m)^2*(1/(2*pi))
power_dep*peak_intensity*beam_detune_tmp*(1/(2*pi))
signal_sq_per_detune=gauss_trap_freq(beam_sigma,dipole_pot,m)^2/beam_detune_tmp

36e-6*2*(1/signal_sq_per_detune)

%%
x_max=10e-3/trap_afreq;
%gauss_amp=1e-20;
excitaiton_energy=pot_fun(x_max,trap_afreq,gauss_amp,gauss_sigma);
period_intagrand=@(x) 1./sqrt(excitaiton_energy-pot_fun(x,trap_afreq,gauss_amp,gauss_sigma) );
num_period=sqrt(2*m)*integral(period_intagrand,-x_max,x_max);
num_freq=1/num_period ;








%%

function afreq=gauss_trap_freq(sigma,amp,mass)
    curvature=-amp* (1./( sigma.^2) );
    afreq=sqrt((1/mass)*curvature);
end
