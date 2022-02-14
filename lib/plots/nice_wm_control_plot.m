%% nice wavemeter control plot

% clf
% plot(data.labview.setpoint-data.labview.setpoint(1),'xb')
% hold on
% plot(data.wm_log.proc.probe.freq.act.mean*1e6-data.labview.setpoint(1),'rx')

%%



%%
font_name='cmr10';
linewidth=1.5;
font_size=12;

colors_main=[[98,136,63];[95,109,187];[180,72,117]]./255;
hsv=colorspace('RGB->HSV',colors_main(:,:));
hsv(:,2)=hsv(:,2);
colors_shaded=colorspace('HSV->RGB',hsv);

stfig('camera fit B')
clf
set(gca, 'FontName', font_name)
set(gca, 'FontSize', font_size)


% [fit_params_out,resnorm,residual,exitflag]=lsqcurvefit(,fit_params_guess,...
%     pos_val,subframe,options);

    
%subplot(2,4,6)
%imagesc(xvals_sub,yvals_sub,D2GaussFunctionRot(fit_params_out,pos_data)) %plot fit


subplot(1,2,1)
%x line profile

%
x_offset=-4.3;
x_lim=[520,590];
x_lim_ref=[x_lim(1),x_lim(1)]+[20,50];
x_lim_inset1=[400,1900];
x_lim_inset2=[x_lim(1),x_lim(1)]+[3.8,7];

freq_meas=data.wm_log.raw.feedback.actual;
freq_set=data.wm_log.raw.feedback.setpt;
freq_err=freq_meas-freq_set;
time=data.wm_log.raw.feedback.posix_time;
time=time-time(1);
plot_mask=time>x_lim(1) & time<x_lim(2);
time_masked=time(plot_mask);
freq_meas_masked=freq_meas(plot_mask);
freq_set_masked=freq_set(plot_mask);
freq_err_masked=freq_err(plot_mask);

mask_inset=time>x_lim_ref(1) & time<x_lim_ref(2);
freq_ref=mean(freq_meas(mask_inset));


plot(time_masked-time_masked(1)+x_offset,freq_set_masked-freq_ref,'--','LineWidth',linewidth,'Color',colors_main(3,:))
hold on
plot(time_masked-time_masked(1)+x_offset,freq_meas_masked-freq_ref,'-','LineWidth',linewidth,'Color',colors_main(2,:))
hold off
xlabel('Time (s)')
ylabel('Frequency (MHz)')
box on
lnh=legend('set point','measurment');
legend('Location','northwest')
lnh.Position=lnh.Position-[0.07,0.0,0,0];
xlim(x_lim-time_masked(1)+x_offset)

set(gca,'linewidth', 1.1)
set(gca,'TickLength',[0.015,0])
set(gca, {'XColor', 'YColor'}, {'k', 'k'});

sp_pos=get(gca,'Position');
an=annotation('textbox', [sp_pos(1)+0.285, sp_pos(2)+sp_pos(4)*0.90, 0, 0], 'string', '(A)','FontSize',round(font_size*1.2),'FontName',font_name);


axes('Position',[.2 .25 .2 .24])
yfactor=1e-3;
mask_inset=time>x_lim_inset1(1) & time<x_lim_inset1(2);
time_masked_inset=time(mask_inset);
freq_meas_masked_inset=freq_meas(mask_inset);
freq_set_masked_inset=freq_set(mask_inset);
freq_ref_inset=mean(freq_meas_masked_inset);
plot(time_masked_inset-time_masked_inset(1),(freq_meas_masked_inset-freq_ref_inset)*yfactor,'LineWidth',linewidth*0.8,'Color',colors_main(2,:))
xlabel('Time (s)')
ylabel('Frequency (GHz)')
xlim(x_lim_inset1-time_masked_inset(1))
ylim([-2.1,2.5])

axes('Position',[.24 .65 .15 .24])

mask_inset2=time>x_lim_inset2(1) & time<x_lim_inset2(2);
time_masked_inset=time(mask_inset2);
freq_meas_masked_inset=freq_meas(mask_inset2);
freq_set_masked_inset=freq_set(mask_inset2);
plot(time_masked_inset-time_masked(1)+x_offset,freq_set_masked_inset-freq_set_masked_inset(1),'--','LineWidth',linewidth,'Color',colors_main(3,:))
hold on
plot(time_masked_inset-time_masked(1)+x_offset,freq_meas_masked_inset-freq_set_masked_inset(1),'-','LineWidth',linewidth,'Color',colors_main(2,:))
hold off
xlabel('Time (s)')
ylabel('Frequency (MHz)')
xlim([time_masked_inset(1),time_masked_inset(end)]-time_masked(1)+x_offset)
ylim([-20,370])


subplot(1,2,2)
std_region=[20,55];
std_region_std=[x_lim(1),x_lim(1)]+std_region;
std_region_mask=time>std_region_std(1) & time<std_region_std(2);
err_mean=mean(freq_err(std_region_mask));
err_std=std(freq_err(std_region_mask));

plot(time_masked-time_masked(1)+x_offset,freq_err_masked,'-','LineWidth',linewidth,'Color',colors_main(2,:))
xline(std_region(1)+x_offset,':','LineWidth',linewidth*1.3,'Color',colors_main(1,:))
xline(std_region(2)+x_offset,':','LineWidth',linewidth*1.3,'Color',colors_main(1,:))
yline(err_mean+err_std,'--','LineWidth',linewidth*1.3,'Color',colors_main(3,:))
yline(err_mean-err_std,'--','LineWidth',linewidth*1.3,'Color',colors_main(3,:))
xlabel('Time (s)')
ylabel('Frequency Error (MHz)')
box on
y_lim=[-1,1]*2;
ylim(y_lim)
legend('error','meas region','','std')
legend('Location','south')
xlim(x_lim-time_masked(1)+x_offset)

sp_pos=get(gca,'Position');
an=annotation('textbox', [sp_pos(1)+0.285, sp_pos(2)+sp_pos(4)*0.90, 0, 0], 'string', '(B)','FontSize',round(font_size*1.2),'FontName',font_name);


set(gca,'linewidth', 1.0)
set(gca,'TickLength',[0.015,0])
set(gca, {'XColor', 'YColor'}, {'k', 'k'});

%



set(gcf,'Units','Pixels')
%set(gcf,'Position',[1068,355,676,453])
fig_width_px=1000;
fig_aspect_ratio=0.3; %0.67;
set(gcf,'Position',[-1800,355,fig_width_px,fig_width_px*fig_aspect_ratio])       
hold off
%%

fig_name='laser_wm_feedback_meas';
%fig_name='pal_mean_v_acc_dyn_static';
fig_dir='./figs/thesis_figs';
export_fig(fullfile(fig_dir,strcat(fig_name,'.svg')))
export_fig(fullfile(fig_dir,strcat(fig_name,'.pdf')))

