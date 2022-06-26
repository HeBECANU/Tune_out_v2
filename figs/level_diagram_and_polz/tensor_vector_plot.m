% plot the vector and tensor components of the polarizability
% want to justify assuming tensor is constant


fsamp=linspace(180,800,1e6)*1e12;
[polz_au,polz_si,higher_terms]=aprox_he_polz(fsamp);



font_name='cmr10';
linewidth=1.5;
font_size=12;
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

yscale=1e-38;
scale_vector=1e3;
scale_tensor=1e4;
scale_vector=10^round(log10(scale_vector));
scale_tensor=10^round(log10(scale_tensor));
fh=stfig('polz_terms si');
clf
plot(fsamp*1e-12,higher_terms.si_polz.scalar/yscale,'LineWidth',linewidth,'Color',colors_main(1,:))
hold on
plot(fsamp*1e-12,higher_terms.si_polz.vector*scale_vector/yscale,'--','LineWidth',linewidth,'Color',colors_main(2,:))
plot(fsamp*1e-12,higher_terms.si_polz.tensor*scale_tensor/yscale,':','LineWidth',linewidth,'Color',colors_main(3,:))
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
legend('Location','best')
ylabel('$\alpha^{(s,v,T)}$, Polarizability ($10^{-38} \mathrm{C}\cdot \mathrm{m}^{2} \cdot\mathrm{V}^{-1}$)')



%%
pause(0.1)
fig_name='polz_scal_vec_tens';
fig_dir='./figs/level_diagram_and_polz';
export_fig(fullfile(fig_dir,strcat(fig_name,'.svg')))
%export_fig(fullfile(fig_dir,strcat(fig_name,'.pdf')))
exportgraphics(fh,fullfile(fig_dir,strcat(fig_name,'.pdf')))


%%

% tune out freq
fminopt = optimset('TolX',1,'TolFun',1e-12);  
tune_out_f_approx=fminsearch(@(x) abs(aprox_he_polz(x)),700e12,fminopt);  

fsamp=tune_out_f_approx+linspace(-10,1.5,1e6)*1e12;
[polz_au,polz_si,higher_terms]=aprox_he_polz(fsamp);


font_name='cmr10';
linewidth=1.5;
font_size=12;
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

yscale=1e-40;
scale_vector=1e3;
scale_tensor=1e4;
scale_vector=10^round(log10(scale_vector));
scale_tensor=10^round(log10(scale_tensor));
fh=stfig('polz_terms si narrow');
clf
plot((fsamp-tune_out_f_approx)*1e-12,higher_terms.si_polz.scalar/yscale,'LineWidth',linewidth,'Color',colors_main(1,:))
hold on
plot((fsamp-tune_out_f_approx)*1e-12,higher_terms.si_polz.vector*scale_vector/yscale,'--','LineWidth',linewidth,'Color',colors_main(2,:))
plot((fsamp-tune_out_f_approx)*1e-12,higher_terms.si_polz.tensor*scale_tensor/yscale,':','LineWidth',linewidth,'Color',colors_main(3,:))
%plot(fsamp*1e-12,polz_si/yscale,'-','LineWidth',linewidth,'Color',colors_main(4,:))


hold off
xlim(([fsamp(1),fsamp(end)]-tune_out_f_approx)*1e-12)

ylim([-3,1])
line(xlim(), [0,0], 'LineWidth', 1, 'Color', 'k','LineStyle','--');
xlabel('Frequency, $f$ (THz)')
ylabel('Polarizability, $\alpha$ (SI)')
legend('Scalar, $\alpha^{s}$', ...
    sprintf('Vector, $\\alpha^{v} \\cdot 10^{%.0g}$',log10(scale_vector)), ...
     sprintf('Tensor, $\\alpha^{T} \\cdot 10^{%.0g}$',log10(scale_tensor)), ...
    '')
legend('Box','off')
legend('FontSize',font_size*0.9)
legend('Location','best')
ylabel('$\alpha^{(s,v,T)}$, Polarizability ($10^{-40} \mathrm{C}\cdot \mathrm{m}^{2} \cdot\mathrm{V}^{-1}$)')

%%
query_pts=[-5,5]*1e9;
polz_vector_scan_lims=interp1((fsamp-tune_out_f_approx),higher_terms.si_polz.vector,query_pts);
polz_tensor_scan_lims=interp1((fsamp-tune_out_f_approx),higher_terms.si_polz.tensor,query_pts);

frac_diff(polz_vector_scan_lims(1),polz_vector_scan_lims(2))
frac_diff(polz_tensor_scan_lims(1),polz_tensor_scan_lims(2))


%%
pause(0.1)
fig_name='polz_scal_vec_tens_narrow';
fig_dir='./figs/level_diagram_and_polz';
export_fig(fullfile(fig_dir,strcat(fig_name,'.svg')))
%export_fig(fullfile(fig_dir,strcat(fig_name,'.pdf')))
exportgraphics(fh,fullfile(fig_dir,strcat(fig_name,'.pdf')))
