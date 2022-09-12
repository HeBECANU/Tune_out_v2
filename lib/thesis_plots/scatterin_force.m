
addpath('../../lib/Core_BEC_Analysis/lib/') %add the path to set_up_project_path
set_up_project_path('../..')

hebec_constants
close('all')

%%
% scat rate per decay rate
% detune relative to the decay rate
% intensity relative to the sat intensity
scatt_rate=@(rel_int,rel_detune) (1/2).*(rel_int./(1+((2*rel_detune).^2)+ rel_int));

%norm_vel=vel*gamma/wave_vec
%norm_detune=delta/gamma
% scattering force is force/( \Gamma hbar * wave vector )
scatt_force=@(norm_vel,norm_detune,rel_int) -scatt_rate(rel_int,norm_detune-norm_vel) + scatt_rate(rel_int,norm_detune+norm_vel);

%%

colors_main=[[[214,72,154];
[102,181,69];
[70,90,214];
[208,153,44];
[212,73,58]]]./255;
hsv=colorspace('RGB->HSV',colors_main(:,:));
hsv(:,2)=hsv(:,2);
colors_shaded=colorspace('HSV->RGB',hsv);
font_name='cmr10';
linewidth=1.5;
font_size=11;


%%




vel_sample=linspace(-1,1,1e4)*10;






stfig('ahn gauss')
clf
set(gca, 'FontName', font_name)
set(gca, 'FontSize', font_size)
%
subplot(1,2,1)
scatt_foce_sample=scatt_force(vel_sample,1,1);
plot(vel_sample,scatt_foce_sample,'--','LineWidth',linewidth,'Color',colors_main(1,:))
hold on
scatt_foce_sample=scatt_force(vel_sample,1,10);
plot(vel_sample,scatt_foce_sample,'-','LineWidth',linewidth,'Color',colors_main(2,:))
scatt_foce_sample=scatt_force(vel_sample,1,30);
plot(vel_sample,scatt_foce_sample,':','LineWidth',linewidth,'Color',colors_main(3,:))
scatt_foce_sample=scatt_force(vel_sample,1,100);
plot(vel_sample,scatt_foce_sample,'-.','LineWidth',linewidth,'Color',colors_main(4,:))
% yline(ho_comb_afreq*yfactor,'--','LineWidth',linewidth*1.2,'Color',colors_main(3,:))
% yline(trap_afreq*yfactor,':','LineWidth',linewidth*1.2,'Color',colors_main(4,:))
% xline(gauss_sigma*2*xfactor,'-.','LineWidth',linewidth*1.2,'Color',colors_main(2,:))
% xline(osc_amp_x_operation*xfactor,'--','LineWidth',linewidth*1.5,'Color',colors_main(5,:))
hold off
xlabel('Norm. Veclocity $v\cdot 2\pi f/(c \Gamma)$')
ylabel('Norm. Acceleration $F \cdot  m c /(f h \Gamma )$')
legend('$I=I_s$', ...
    '$I=10\: I_s$', ...
    '$I=30 \: I_s$', ...
    '$I=100 \: I_s$')
legend('Location','best')
sp_pos=get(gca,'Position');
an=annotation('textbox', [sp_pos(1)+0.02, sp_pos(2)+sp_pos(4)*0.95, 0, 0], 'string', '(A)','FontSize',round(font_size*1.5),'FontName',font_name);
set(gca, {'XColor', 'YColor'}, {'k', 'k'});
set(gca,'linewidth', 1.0)
set(gca,'TickLength',[0.015,0])
set(gca, 'FontName', font_name)
set(gca, 'FontSize', font_size)

subplot(1,2,2)
scatt_foce_sample=scatt_force(vel_sample,0.5,10);
plot(vel_sample,scatt_foce_sample,':','LineWidth',linewidth,'Color',colors_main(1,:))
hold on
scatt_foce_sample=scatt_force(vel_sample,1,10);
plot(vel_sample,scatt_foce_sample,'--','LineWidth',linewidth,'Color',colors_main(3,:))
scatt_foce_sample=scatt_force(vel_sample,2,10);
plot(vel_sample,scatt_foce_sample,'-','LineWidth',linewidth,'Color',colors_main(2,:))
scatt_foce_sample=scatt_force(vel_sample,5,10);
plot(vel_sample,scatt_foce_sample,'-.','LineWidth',linewidth,'Color',colors_main(4,:))
% yline(ho_comb_afreq*yfactor,'--','LineWidth',linewidth*1.2,'Color',colors_main(3,:))
% yline(trap_afreq*yfactor,':','LineWidth',linewidth*1.2,'Color',colors_main(4,:))
% xline(gauss_sigma*2*xfactor,'-.','LineWidth',linewidth*1.2,'Color',colors_main(2,:))
% xline(osc_amp_x_operation*xfactor,'--','LineWidth',linewidth*1.5,'Color',colors_main(5,:))
hold off
xlabel('Norm. Veclocity $v\cdot 2\pi f/(c \Gamma)$')
ylabel('Norm. Acceleration $F \cdot  m c /(f h \Gamma )$')
legend('$\Delta=\Gamma/2$', ...
    '$\Delta=-\Gamma$', ...
   '$\Delta=-2\Gamma$', ...
    '$\Delta=-5\Gamma$')
legend('Location','best')
sp_pos=get(gca,'Position');
an=annotation('textbox', [sp_pos(1)+0.02, sp_pos(2)+sp_pos(4)*0.95, 0, 0], 'string', '(B)','FontSize',round(font_size*1.5),'FontName',font_name);
set(gca, {'XColor', 'YColor'}, {'k', 'k'});
set(gca,'linewidth', 1.0)
set(gca,'TickLength',[0.015,0])
set(gca, 'FontName', font_name)
set(gca, 'FontSize', font_size)

set(gcf,'Units','Pixels')
%set(gcf,'Position',[1068,355,676,453])
fig_width_px=1200;
fig_aspect_ratio=0.35; %0.67;
set(gcf,'Position',[700,355,fig_width_px,fig_width_px*fig_aspect_ratio])       
hold off

%%
fig_name='laser_cool_force';
%fig_name='pal_mean_v_acc_dyn_static';
fig_dir='./figs/thesis_figs';
export_fig(fullfile(fig_dir,strcat(fig_name,'.svg')))
pause(0.1)
export_fig(fullfile(fig_dir,strcat(fig_name,'.pdf')))



%% scale conversion
gamma=2*pi*1.62e6;
f=f2wl(1083e-9);
yscale=const.c*const.mhe/(f*const.h*gamma);
1/yscale
xcale=2*pi*f/(const.c*gamma)
1/xcale

fprintf('1 on the x scale corresponds to %.3f ,m/s\n',1/xcale)
fprintf('1 on the y scale corresponds to %.3g m/s^2 \n',1/yscale)