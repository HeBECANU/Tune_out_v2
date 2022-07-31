% plot the vector and tensor components of the polarizability
% want to justify assuming tensor is constant


fsamp=linspace(180,800,1e6)*1e12;
[polz_au,polz_si,higher_terms]=aprox_he_polz(fsamp);



font_name='cmr10';
linewidth=1.5;
font_size=13;
colors_main=[[[102,181,69];
[212,73,58];
[214,72,154];
[208,153,44];
[151,90,214];]]./255;
hsv=colorspace('RGB->HSV',colors_main(:,:));
hsv(:,2)=hsv(:,2)-0.3;
hsv(:,3)=hsv(:,3)+0.25;
colors_shaded=colorspace('HSV->RGB',hsv);

% check that the new calcluation agrees
% stfig('polz au');
% clf
% plot(fsamp*1e-12,polz_au)
% yrange=1e6;
% ylim([-1,1]*yrange)
% yline(0)
% xlim([fsamp(1),fsamp(end)]*1e-12)
% xlabel('Frequency, $f$ (THz)')
% ylabel('Polarizability, $\alpha$ (atom. u.)')

%
fh_wide=stfig('polz_terms si');
clf


yscale=1e-38;
yscale=10^round(log10(yscale));
scale_vector=1e3;
scale_tensor=1e6;
scale_vector=10^round(log10(scale_vector));
scale_tensor=10^round(log10(scale_tensor));

plot(fsamp*1e-12,higher_terms.si_polz.scalar/yscale,'LineWidth',linewidth,'Color',colors_main(1,:))
hold on
plot(fsamp*1e-12,higher_terms.si_polz.vector*scale_vector/yscale,'--','LineWidth',linewidth,'Color',colors_main(2,:))
plot(fsamp*1e-12,higher_terms.si_polz.tensor*scale_tensor/yscale,':','LineWidth',linewidth,'Color',colors_main(5,:))
%plot(fsamp*1e-12,polz_si/yscale,'-','LineWidth',linewidth,'Color',colors_main(4,:))


hold off
xlim([fsamp(1),fsamp(end)]*1e-12)
yrange=1;
ylim([-1,1]*yrange)
line(xlim(), [0,0], 'LineWidth', 1, 'Color', 'k','LineStyle','--');
xlabel('Frequency, $f$ (THz)')
ylabel('Polarizability, $\alpha$ (SI)')
legend('Scalar, $\alpha^{s}$', ...
    sprintf('Vector, $\\alpha^{v} \\cdot 10^{%.0g}$',log10(scale_vector)), ...
     sprintf('Tensor, $\\alpha^{T} \\cdot 10^{%.0g}$',log10(scale_tensor)), ...
    '')
legend('Box','off')
legend('FontSize',font_size*0.9)
ln=legend('Location','south');
ln.Position=ln.Position+[0.1,0,0,0];
ylabel(sprintf('$\\alpha^{(s,v,T)}$, Polarizability ($10^{%.0f} \\mathrm{C}\\cdot \\mathrm{m}^{2} \\cdot\\mathrm{V}^{-1}$)',log10(yscale)))
set(gca, 'FontName', font_name)
set(gca, 'FontSize', font_size)

%%
%%
set(gcf,'Units','Pixels')
%set(gcf,'Position',[1068,355,676,453])
fig_width_px=600;
fig_aspect_ratio=0.67; %0.67;
set(gcf,'Position',[500,355,fig_width_px,fig_width_px*fig_aspect_ratio])
%%
pause(0.1)
fig_name='polz_scal_vec_tens';
fig_dir='./figs/level_diagram_and_polz';
export_fig(fullfile(fig_dir,strcat(fig_name,'.svg')),fh_wide)
%export_fig(fullfile(fig_dir,strcat(fig_name,'.pdf')))
exportgraphics(fh_wide,fullfile(fig_dir,strcat(fig_name,'.pdf')))


%%

% tune out freq
fminopt = optimset('TolX',1,'TolFun',1e-12);  
tune_out_f_approx=fminsearch(@(x) abs(aprox_he_polz(x)),700e12,fminopt);  

fsamp=tune_out_f_approx+linspace(-10,5,1e6)*1e12;
[polz_au,polz_si,higher_terms]=aprox_he_polz(fsamp);


% check that the new calcluation agrees
% stfig('polz au');
% clf
% plot(fsamp*1e-12,polz_au)
% yrange=1e6;
% ylim([-1,1]*yrange)
% yline(0)
% xlim([fsamp(1),fsamp(end)]*1e-12)
% xlabel('Frequency, $f$ (THz)')
% ylabel('Polarizability, $\alpha$ (atom. u.)')

%
fh_narrow=stfig('polz_terms si narrow');
clf


yscale=1e-40;
scale_vector=1e3;
scale_tensor=1e4;
scale_vector=10^round(log10(scale_vector));
scale_tensor=10^round(log10(scale_tensor));
plot((fsamp-tune_out_f_approx)*1e-12,higher_terms.si_polz.scalar/yscale,'LineWidth',linewidth,'Color',colors_main(1,:))
hold on
plot((fsamp-tune_out_f_approx)*1e-12,higher_terms.si_polz.vector*scale_vector/yscale,'--','LineWidth',linewidth,'Color',colors_main(2,:))
plot((fsamp-tune_out_f_approx)*1e-12,higher_terms.si_polz.tensor*scale_tensor/yscale,':','LineWidth',linewidth,'Color',colors_main(5,:))
%plot(fsamp*1e-12,polz_si/yscale,'-','LineWidth',linewidth,'Color',colors_main(4,:))


hold off
xlim(([fsamp(1),fsamp(end)]-tune_out_f_approx)*1e-12)

ylim([-2,1])
line(xlim(), [0,0], 'LineWidth', 1, 'Color', 'k','LineStyle','--');
xlabel('Frequency, $f$ (THz)')
ylabel('Polarizability, $\alpha$ (SI)')
legend('Scalar, $\alpha^{s}$', ...
    sprintf('Vector, $\\alpha^{v} \\cdot 10^{%.0g}$',log10(scale_vector)), ...
     sprintf('Tensor, $\\alpha^{T} \\cdot 10^{%.0g}$',log10(scale_tensor)), ...
    '')
legend('Box','off')
legend('FontSize',font_size*0.9)
legend('Location','south')
ylabel(sprintf('$\\alpha^{(s,v,T)}$, Polarizability ($10^{%.0f} \\mathrm{C}\\cdot \\mathrm{m}^{2} \\cdot\\mathrm{V}^{-1}$)',log10(yscale)))
set(gca, 'FontName', font_name)
set(gca, 'FontSize', font_size)

%%
query_pts=[-5,5]*1e9;
polz_vector_scan_lims=interp1((fsamp-tune_out_f_approx),higher_terms.si_polz.vector,query_pts);
polz_tensor_scan_lims=interp1((fsamp-tune_out_f_approx),higher_terms.si_polz.tensor,query_pts);

frac_diff(polz_vector_scan_lims(1),polz_vector_scan_lims(2))
frac_diff(polz_tensor_scan_lims(1),polz_tensor_scan_lims(2))


%%
set(gcf,'Units','Pixels')
%set(gcf,'Position',[1068,355,676,453])
fig_width_px=600;
fig_aspect_ratio=0.67; %0.67;
set(gcf,'Position',[500,355,fig_width_px,fig_width_px*fig_aspect_ratio])
%%
pause(0.1)
fig_name='polz_scal_vec_tens_narrow';
fig_dir='./figs/level_diagram_and_polz';
export_fig(fullfile(fig_dir,strcat(fig_name,'.svg')),fh_narrow)
%export_fig(fullfile(fig_dir,strcat(fig_name,'.pdf')))
exportgraphics(fh_narrow,fullfile(fig_dir,strcat(fig_name,'.pdf')))



%% demo sign flips inside manifold

fsamp=linspace(276.72,276.77,1e6)*1e12;
%fsamp=linspace(276.76,276.765,1e3)*1e12;
[polz_au,polz_si,higher_terms]=aprox_he_polz(fsamp);



% check that the new calcluation agrees
% stfig('polz au');
% clf
% plot(fsamp*1e-12,polz_au)
% yrange=1e6;
% ylim([-1,1]*yrange)
% yline(0)
% xlim([fsamp(1),fsamp(end)]*1e-12)
% xlabel('Frequency, $f$ (THz)')
% ylabel('Polarizability, $\alpha$ (atom. u.)')


fh_23p=stfig('polz_terms si in manifold');
clf

%
xoffset=276e12;
xscale=1e9;
yscale=1e-34;
yscale=10^round(log10(yscale));
scale_vector=1e0;
scale_tensor=1e1;
scale_vector=10^round(log10(scale_vector));
scale_tensor=10^round(log10(scale_tensor));

plot((fsamp-xoffset)/xscale,higher_terms.si_polz.scalar/yscale,'LineWidth',linewidth,'Color',colors_main(1,:))
hold on
plot((fsamp-xoffset)/xscale,higher_terms.si_polz.vector*scale_vector/yscale,'--','LineWidth',linewidth,'Color',colors_main(2,:))
plot((fsamp-xoffset)/xscale,higher_terms.si_polz.tensor*scale_tensor/yscale,':','LineWidth',linewidth,'Color',colors_main(5,:))
%plot(fsamp*1e-12,polz_si/yscale,'-','LineWidth',linewidth,'Color',colors_main(4,:))


hold off
%xlim([fsamp(1),fsamp(end)]*1e-12)
yrange=5;
ylim([-1,1]*yrange)
line(xlim(), [0,0], 'LineWidth', 1, 'Color', 'k','LineStyle','--');
xlabel('Frequency, $f-276000$  (GHz)')
ylabel('Polarizability, $\alpha$ (SI)')
legend('Scalar, $\alpha^{s}$', ...
     sprintf('Vector, $\\alpha^{v} \\cdot 10^{%.0g}$',log10(scale_vector)), ...
     sprintf('Tensor, $\\alpha^{T} \\cdot 10^{%.0g}$',log10(scale_tensor)), ...
    '')
legend('Box','off')
legend('FontSize',font_size*0.9)
legend('Location','south')
ylabel(sprintf('$\\alpha^{(s,v,T)}$, Polarizability ($10^{%.0f} \\mathrm{C}\\cdot \\mathrm{m}^{2} \\cdot\\mathrm{V}^{-1}$)',log10(yscale)))
set(gca, 'FontName', font_name)
set(gca, 'FontSize', font_size)

%%
set(gcf,'Units','Pixels')
%set(gcf,'Position',[1068,355,676,453])
fig_width_px=600;
fig_aspect_ratio=0.67; %0.67;
set(gcf,'Position',[500,355,fig_width_px,fig_width_px*fig_aspect_ratio])

%%


pause(0.1)
fig_name='polz_23p_manifold';
fig_dir='./figs/level_diagram_and_polz';
export_fig(fullfile(fig_dir,strcat(fig_name,'.svg')),fh_23p)
%export_fig(fullfile(fig_dir,strcat(fig_name,'.pdf')))
exportgraphics(fh_23p,fullfile(fig_dir,strcat(fig_name,'.pdf')))



%%

%% does the polz diverge

fsamp=linspace(276.72,276.740,1e6)*1e12;
%fsamp=linspace(276.76,276.765,1e3)*1e12;
[polz_au,polz_si,higher_terms]=aprox_he_polz(fsamp);



% check that the new calcluation agrees
% stfig('polz au');
% clf
% plot(fsamp*1e-12,polz_au)
% yrange=1e6;
% ylim([-1,1]*yrange)
% yline(0)
% xlim([fsamp(1),fsamp(end)]*1e-12)
% xlabel('Frequency, $f$ (THz)')
% ylabel('Polarizability, $\alpha$ (atom. u.)')


fh_23p=stfig('polz_terms si in manifold');
clf

%
xoffset=276e12;
xscale=1e9;
yscale=1e-30;
yscale=10^round(log10(yscale));
scale_vector=1e0;
scale_tensor=1e1;
scale_vector=10^round(log10(scale_vector));
scale_tensor=10^round(log10(scale_tensor));

plot((fsamp-xoffset)/xscale,higher_terms.si_polz.scalar/yscale,'LineWidth',linewidth,'Color',colors_main(1,:))
hold on
plot((fsamp-xoffset)/xscale,higher_terms.si_polz.vector*scale_vector/yscale,'--','LineWidth',linewidth,'Color',colors_main(2,:))
plot((fsamp-xoffset)/xscale,higher_terms.si_polz.tensor*scale_tensor/yscale,':','LineWidth',linewidth,'Color',colors_main(5,:))
%plot(fsamp*1e-12,polz_si/yscale,'-','LineWidth',linewidth,'Color',colors_main(4,:))


hold off
%xlim([fsamp(1),fsamp(end)]*1e-12)
yrange=5;
ylim([-1,1]*yrange)
line(xlim(), [0,0], 'LineWidth', 1, 'Color', 'k','LineStyle','--');
xlabel('Frequency, $f-276000$  (GHz)')
ylabel('Polarizability, $\alpha$ (SI)')
legend('Scalar, $\alpha^{s}$', ...
     sprintf('Vector, $\\alpha^{v} \\cdot 10^{%.0g}$',log10(scale_vector)), ...
     sprintf('Tensor, $\\alpha^{T} \\cdot 10^{%.0g}$',log10(scale_tensor)), ...
    '')
legend('Box','off')
legend('FontSize',font_size*0.9)
legend('Location','south')
ylabel(sprintf('$\\alpha^{(s,v,T)}$, Polarizability ($10^{%.0f} \\mathrm{C}\\cdot \\mathrm{m}^{2} \\cdot\\mathrm{V}^{-1}$)',log10(yscale)))
set(gca, 'FontName', font_name)
set(gca, 'FontSize', font_size)
