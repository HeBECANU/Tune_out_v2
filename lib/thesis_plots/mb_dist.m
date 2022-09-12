addpath('../../lib/Core_BEC_Analysis/lib/') %add the path to set_up_project_path
set_up_project_path('../..')

hebec_constants
close('all')

%%

mb_e=@(energy,temp) 2*sqrt(energy./pi).*((1/(const.kb*temp))).^(3/2) .* exp(-energy./(const.kb*temp) );

mb_e(1,1)



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


inital_temperature=0.1e-3;
cut_temperature=inital_temperature*2.0;
% lets find the population below the cut
% and the mean energy of the atoms below the cut

pop_before_cut=integral(@(energy) mb_e(energy,inital_temperature),0,inital_temperature*1000*const.kb,'AbsTol',1e-12);
if abs(frac_diff(pop_before_cut,1))>1e-6
    error('mb dist not norm ')
end
mean_enegy_before_cut=integral(@(energy) energy.*mb_e(energy,inital_temperature),0,inital_temperature*1000*const.kb,'AbsTol',1e-12);
if abs(frac_diff(inital_temperature,mean_enegy_before_cut/((3/2)*const.kb)))>1e-6
    error('not the inital temp')
end
fprintf('\n\n')
fprintf('inital temp %g uk \n',inital_temperature*1e6)
fprintf('cut temp %g uk \n',cut_temperature*1e6)
pop_after_cut=integral(@(energy) mb_e(energy,inital_temperature),0,cut_temperature*(3/2)*const.kb);
fprintf('remaining population %g \n',pop_after_cut)
mean_enegy_after_cut=(1/pop_after_cut)*integral(@(energy) energy.*mb_e(energy,inital_temperature),0,cut_temperature*(3/2)*const.kb,'AbsTol',1e-12);
temp_after_cut=mean_enegy_after_cut/((3/2)*const.kb);
fprintf('temperature %g uK \n',temp_after_cut*1e6)
fprintf('frac change %g \n',temp_after_cut/inital_temperature)

%%

en_sample=linspace(0,1,1e4)*const.kb*0.6e-3;





xfactor=1/((3/2)*const.kb) * 1e3;
yfactor=1e-26;
stfig('forced evap')
clf
set(gca, 'FontName', font_name)
set(gca, 'FontSize', font_size)
%

en_sample_cut=linspace(cut_temperature*((3/2)*const.kb),en_sample(end),1e2);
y_cut_region=mb_e(en_sample_cut,inital_temperature);
ph_removed=patch([en_sample_cut fliplr(en_sample_cut)]*xfactor, [0*en_sample_cut fliplr(y_cut_region)]*yfactor, [1,1,1]*0.8,...
    'EdgeColor','none')
hold on

density_sampl=mb_e(en_sample,inital_temperature);
ph_ini=plot(en_sample*xfactor,density_sampl*yfactor,'-','LineWidth',linewidth,'Color',colors_main(3,:))
density_sampl=mb_e(en_sample,temp_after_cut);
ph_after=plot(en_sample*xfactor,density_sampl*yfactor*pop_after_cut,'--','LineWidth',linewidth,'Color',colors_main(2,:))

% scatt_foce_sample=scatt_force(en_sample,1,10);
% plot(en_sample,scatt_foce_sample,'-','LineWidth',linewidth,'Color',colors_main(2,:))
% scatt_foce_sample=scatt_force(en_sample,1,30);
% plot(en_sample,scatt_foce_sample,':','LineWidth',linewidth,'Color',colors_main(3,:))
% scatt_foce_sample=scatt_force(en_sample,1,100);
% plot(en_sample,scatt_foce_sample,'-.','LineWidth',linewidth,'Color',colors_main(4,:))
% yline(ho_comb_afreq*yfactor,'--','LineWidth',linewidth*1.2,'Color',colors_main(3,:))
% yline(trap_afreq*yfactor,':','LineWidth',linewidth*1.2,'Color',colors_main(4,:))
% xline(gauss_sigma*2*xfactor,'-.','LineWidth',linewidth*1.2,'Color',colors_main(2,:))
ph_cut=xline(cut_temperature*1e3,':','LineWidth',linewidth*1.5,'Color',colors_main(5,:))
hold off
xlabel('Energy$/k_B$ (mK)')
ylabel('Density (arb. u.)')
legend([ph_ini,ph_after,ph_removed,ph_cut], ...
    {'Inital', ...
    'After Cut', ...
    'Removed Popultion', ...
    'Cut Energy'})
legend('Location','best','FontSize',round(font_size*1.1))
sp_pos=get(gca,'Position');
an=annotation('textbox', [sp_pos(1)+sp_pos(3)*0.2, sp_pos(2)+sp_pos(4)*0.7, 0.2, 0],'interpreter','latex', ...
    'color',colors_main(3,:),'string', ...
    sprintf('$T=%.0f$ $\\mu$K',inital_temperature*1e6),'FontSize',round(font_size*1.5),'FontName',font_name,'LineStyle','none');
an=annotation('textbox', [sp_pos(1)+sp_pos(3)*0.11, sp_pos(2)+sp_pos(4)*0.3, 0.2, 0],'interpreter','latex', ...
    'color',colors_main(2,:),'string', ...
    sprintf('$T=%.0f$ $\\mu$K',temp_after_cut*1e6),'FontSize',round(font_size*1.5),'FontName',font_name,'LineStyle','none');
an=annotation('textbox', [sp_pos(1)+sp_pos(3)*0.6, sp_pos(2)+sp_pos(4)*0.2, 0.2, 0],'interpreter','latex', ...
    'color',[1,1,1]*0.6,'string', ...
    sprintf('$N=%.3f$',1-pop_after_cut),'FontSize',round(font_size*1.5),'FontName',font_name,'LineStyle','none');
set(gca, {'XColor', 'YColor'}, {'k', 'k'});
set(gca,'linewidth', 1.0)
set(gca,'TickLength',[0.015,0])
set(gca, 'FontName', font_name)
set(gca, 'FontSize', font_size)
xlim([0,max(en_sample)]*xfactor)
ylim([0,4.5])
box on


set(gcf,'Units','Pixels')
%set(gcf,'Position',[1068,355,676,453])
fig_width_px=800;
fig_aspect_ratio=0.5; %0.67;
set(gcf,'Position',[700,355,fig_width_px,fig_width_px*fig_aspect_ratio])       
hold off

%%
fig_name='forced_evap';
%fig_name='pal_mean_v_acc_dyn_static';
fig_dir='./figs/thesis_figs';
exportgraphics(gcf,fullfile(fig_dir,strcat(fig_name,'.eps')))
%export_fig(fullfile(fig_dir,strcat(fig_name,'.svg')))
pause(0.1)
export_fig(fullfile(fig_dir,strcat(fig_name,'.pdf')))
