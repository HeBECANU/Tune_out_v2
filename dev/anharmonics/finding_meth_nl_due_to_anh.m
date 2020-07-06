%finding_meth_nl_due_to_anh


% set upt he options for the numerical integeration
opts_template=[];
opts_template.trap_freq=424;
opts_template.mass=const.mhe;
opts_template.gauss_radius=9e-6*2;
opts_template.gauss_displacement=2e-6*2;
opts_template.find_new_min=true;
opts_template.plot_find_new_min=true;
%opts.tlims=[0,(1/opts.trap_freq)*100];
opts_template.tlims=[0,1.1];
% damping_ratio=1/(trap_freq*damping_time);
% damping_coef=damping_ratio*2*sqrt(mass*trap_derivs(3));
opts_template.damping_time=0.5; % 0.6;
opts_template.inital_state=[0,13e-3];
beam_power=80e-3;

osc_e_in_j=(1/2)*const.mhe*opts_template.inital_state(2)^2;
osc_e_in_k=osc_e_in_j/const.kb;
fprintf('oscillation amplitude %.2f mm/s \n', 1e3* opts_template.inital_state(2))
fprintf('            energy    %.2f nK \n', 1e9*osc_e_in_k)
% wine bottle threshold
% when the net trap frequency is zero
dpolz_domega=1.79e-53;
power=beam_power;
detuning_hz=-20e9;
probe_trap_freq=sqrt((1/(2*const.epsilon0*const.c*const.mhe))...
                    *detuning_hz*dpolz_domega...
                    *( (2*power)/(pi*(opts_template.gauss_radius^2)) ) ...
                    *( (4)/(opts_template.gauss_radius^2) ) );
abs(probe_trap_freq)
probe_trap_freq/opts_template.trap_freq
power*detuning_hz


%%
do_plots_of_osc=false;
%beam_detunings=col_vec(linspace(-4e9,4e9,5));
%beam_detunings=col_vec(linspace(-2e9,6e9,100));
beam_detunings=col_vec(linspace(-50e9,50e9,20)); % detuning from the TO in Hz
beam_detunings=cat(1,0,beam_detunings); %add on the zero detuning case as the baseline
iimax=numel(beam_detunings);
fit_results=cell(iimax,1);
parfor_progress_imp(iimax);
for ii=1:iimax
    %fprintf('simulating motion %04u of %04u \n',ii,iimax)
    this_fit_st=[];
    
    opts=opts_template;
    
    % find the potential dpeth of the gaussian from the power and detuning

    %detuning_hz=-3e9; % best uncert from measurement -0.03e9;
    detuning_hz=beam_detunings(ii);
    this_fit_st.probe.detuning_hz=detuning_hz;
    this_fit_st.probe.power_w=beam_power;
    u_dip_0=-(1/(2*const.epsilon0*const.c))...
            *detuning_hz*dpolz_domega...
            *( (2*power)/(pi*(opts.gauss_radius^2)) ) ;
    opts.gauss_amp=u_dip_0;
    
    motion_solution=motion_in_harmonic_w_gaussian(opts);
    this_fit_st.motion_solution=motion_solution;

    % sample like the atom laser does
    n_pulses=136;
    pulsedt=8.000e-3;
    t_al_pulse=0.01 +(0:136)*pulsedt ;
    response = interp1(motion_solution.x, motion_solution.y(2,:), t_al_pulse);
    predictor=t_al_pulse;
    
    this_fit_st.al_samp.time=predictor;
    this_fit_st.al_samp.vel=response;

    %

    % err_fit_dim=xyzerr_tmp(:,anal_opts_osc_fit.dimesion);
    % weights=1./(err_fit_dim.^2);
    % weights(isnan(weights))=1e-20; %if nan then set to smallest value you can
    % weights=weights/sum(weights);

    dom_opt=[];
    dom_opt.components_min_amp=1e-3;
    dom_out=dominant_freq_components(predictor,response,dom_opt);
    fit_freq=dom_out.freq(1);
    fit_phase=dom_out.phase(1);
    fit_amp=dom_out.amp(1);
    %fft_phase=angle(out(2,nearest_idx))+0.535;

    % do a global fit to get a good place to start the fit routine
    % before we call fitnlm we will try to get close to the desired fit parameters using a robust global
    % optimizer
    modelfun = @(b,x) exp(-x.*max(0,b(6))).*b(1).*sin(b(2)*x*pi*2+b(3)*pi*2)+b(4)+b(5)*x;
    % simple model
    %modelfun_simple = @(b,x) exp(-x(:,1).*max(0,b(6))).*b(1).*sin(b(2)*x(:,1)*pi*2+b(3)*pi*2)+b(4)+b(5)*x(:,1);
    gf_opt=[];
    gf_opt.domain=[[1,50]*1e-3;...   %amp
                   [10,80];...       %freq
                   [-2,2]*pi;...     %phase
                   [-50,50]*1e-3;... %offset
                   [-10,10]*1e-3;... %gradient
                   [-1e-2,10];...    %damp
                   ];        
    gf_opt.start=[fit_amp, fit_freq, fit_phase, 0,0,2];
    gf_opt.rmse_thresh=2e-3;
    gf_opt.plot=false;
    gf_out=global_fit(predictor,response,modelfun,gf_opt);
    
    
    % do the fit with the start conditions from global_fit
    cof_names={'amp','freq','phase','offset','grad','damp'};
    opt = statset('TolFun',1e-10,'TolX',1e-10,...
        'MaxIter',1e4,... %1e4
        'UseParallel',1);
    beta0=gf_out.params;

    % 'Weights',weights,
    fitobject=fitnlm(predictor,response,modelfun,beta0,...
        'options',opt,...
        'CoefficientNames',cof_names);
    this_fit_st.osc_fit.model=fitobject;
    this_fit_st.osc_fit.fitparam=fitobject.Coefficients;
    this_fit_st.osc_fit.rmse=fitobject.RMSE;
    % we will later use this fit to find the undamped frequency of the oscillator

    fit_freq_unfolded=3*(1/pulsedt)+fitobject.Coefficients.Estimate(2);
    damping_ratio=fitobject.Coefficients.Estimate(6)./(2*pi*fit_freq_unfolded);
    undamp_fit_freq=fit_freq_unfolded./sqrt(1-(damping_ratio).^2);

    this_fit_st.osc_fit.fit_freq_unfolded=fit_freq_unfolded;
    this_fit_st.osc_fit.damping_ratio=damping_ratio;
    this_fit_st.osc_fit.undamp_fit_freq=undamp_fit_freq;
    % write the output
    fit_results{ii}=this_fit_st;
    %%
    if do_plots_of_osc
        stfig('motion in trap');
        subplot(2,1,1)
        plot(motion_solution.x,motion_solution.y(1,:))
        xlabel('time (s)')
        ylabel('position')
        subplot(2,1,2)
        plot(motion_solution.x,motion_solution.y(2,:))
        ylabel('velocity')


        font_name='cmr10';
        font_size_global=10;
        font_size_label=10;
        colors_main=[[255,0,0];[33,188,44];[0,0,0]]./255;
        lch=colorspace('RGB->LCH',colors_main(:,:));
        lch(:,1)=lch(:,1)+20;
        colors_detail=colorspace('LCH->RGB',lch);
        %would prefer to use srgb_2_Jab here
        color_shaded=colorspace('RGB->LCH',colors_main(3,:));
        color_shaded(1)=50;
        color_shaded=colorspace('LCH->RGB',color_shaded);


        detuning_samp_vals=linspace(min(predictor),...
            max(predictor),1e5)';
        [prediction,ci]=predict(fitobject,detuning_samp_vals,'Alpha',1-erf(1/sqrt(2)),'Prediction','observation');
        stfig('fit osc');
        clf;
        time_start=min(predictor);

        shaded_ci_lines=false;

        if shaded_ci_lines
            patch([predictor', fliplr(predictor')]-time_start, [ci(:,1)', fliplr(ci(:,2)')]*1e3, color_shaded,'EdgeColor','none')  %[1,1,1]*0.80
            hold on
        else
            plot(detuning_samp_vals-time_start,ci(:,1)*1e3,'-','LineWidth',1.5,'Color',color_shaded)
            hold on
            plot(detuning_samp_vals-time_start,ci(:,2)*1e3,'-','LineWidth',1.5,'Color',color_shaded)
        end  
        plot(detuning_samp_vals-time_start,prediction*1e3,'-','LineWidth',1.0,'Color',colors_main(3,:))
        ax = gca;
        set(ax, {'XColor', 'YColor'}, {'k', 'k'});
        plot(predictor-time_start,response*1e3,'o','MarkerSize',5,'Color',colors_main(1,:),'MarkerFaceColor',colors_detail(1,:),'LineWidth',1.5) 
        set(gcf,'Color',[1 1 1]);
        xlabel('Time (s)','FontSize',font_size_label)
        ylabel('$V_{x}$ (mm/s)','FontSize',font_size_label)
        title(sprintf('amp=%.3f$\\pm$%.3f mm/s,$\\omega$=$2\\pi$ %.3f$\\pm$%.3f Hz,$\\lambda$=%.3f$\\pm$%.3f s',...
            fitobject.Coefficients.Estimate(1)*1e3,fitobject.Coefficients.SE(1)*1e3,...
             undamp_fit_freq,fitobject.Coefficients.SE(2),...
             fitobject.Coefficients.Estimate(6),fitobject.Coefficients.SE(6)),'Interpreter','latex')
        hold off
        ax = gca;
        %xlim([predictorplot(1,1)-0.01,predictorplot(end,1)+0.01]-time_start)
        set(ax, {'XColor', 'YColor'}, {'k', 'k'});
        set(gca,'linewidth',1.0)
        set(gca,'FontSize',font_size_global,'FontName',font_name)
        %saveas(gca,sprintf('%sfit_dld_shot_num%04u.png',anal_opts_osc_fit.global.out_dir,dld_shot_num))
        pause(0.01)
    end
    
    parfor_progress_imp;
end
parfor_progress_imp(0);


%%


fit_results_arr=cell_array_of_struct_to_struct_of_array(fit_results);

if (fit_results_arr.probe.detuning_hz(1)*fit_results_arr.probe.power_w(1))~=0
    error('first value should be zero')
end
probe_freq_squared=fit_results_arr.osc_fit.undamp_fit_freq(2:end).^2- ...
                        fit_results_arr.osc_fit.undamp_fit_freq(1).^2 ;
detuning=fit_results_arr.probe.detuning_hz(2:end);

set(0,'defaulttextInterpreter','latex')
stfig('method linearity');
xscale=1e-9;
plot(detuning*xscale,probe_freq_squared)
xlabel('$\omega-\omega_{\mathrm{TO}}$ ($2\pi$ GHz)')
ylabel('$\Omega_{\mathrm{Net}}^2-\Omega_{\mathrm{Trap}}^2$ (Hz$^2$)')


% fit the dependence a linear function

predictor=detuning*xscale;
response=probe_freq_squared;

fit_in=cat(2,predictor,response,response*nan);
meth_lin_fit=fit_poly_with_int(fit_in,1,0,0);
fprintf('fit intercept %f MHz \n',meth_lin_fit.x_intercept.val*1e3)
%

font_name='cmr10';
font_size_global=10;
font_size_label=10;
colors_main=[[255,0,0];[33,188,44];[0,0,0]]./255;
lch=colorspace('RGB->LCH',colors_main(:,:));
lch(:,1)=lch(:,1)+20;
colors_detail=colorspace('LCH->RGB',lch);
%would prefer to use srgb_2_Jab here
color_shaded=colorspace('RGB->LCH',colors_main(3,:));
color_shaded(1)=50;
color_shaded=colorspace('LCH->RGB',color_shaded);

detuning_samp_vals=linspace(min(predictor),...
                    max(predictor),1e5)';
[prediction,ci]=predict(meth_lin_fit.fit_mdl,detuning_samp_vals,'Alpha',1-erf(1/sqrt(2)),'Prediction','observation');

stfig('polz from osc freq fit')
clf;
shaded_ci_lines=false;
if shaded_ci_lines
    patch([predictor', fliplr(predictor')], [ci(:,1)', fliplr(ci(:,2)')], color_shaded,'EdgeColor','none')  %[1,1,1]*0.80
    hold on
else
    plot(detuning_samp_vals,ci(:,1),'-','LineWidth',1.5,'Color',color_shaded)
    hold on
    plot(detuning_samp_vals,ci(:,2),'-','LineWidth',1.5,'Color',color_shaded)
end  
plot(detuning_samp_vals,prediction,'-','LineWidth',1.0,'Color',colors_main(3,:))
plot(predictor,response,'x')


%% fit the dependence with a higher order 
