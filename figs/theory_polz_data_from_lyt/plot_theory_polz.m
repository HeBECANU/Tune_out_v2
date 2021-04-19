path='./figs/theory_polz_data_from_lyt/fig1b_data.csv';
data_table=readtable(path)

bohr_radius=5.29177210903e-11;

%%
stfig('theory polz');
clf

wl_axis=axes('Position',[.1 .1 .8 1e-12]);
set(wl_axis,'Units','normalized');
set(wl_axis,'Color','none');
xlabel(wl_axis,'Wavelength (nm)')

% axis for km/h with stem-plot
freq_axis=axes('Position',[.1 .25 .8 .65]);
set(freq_axis,'Units','normalized');

set(wl_axis,'FontSize',10);

% conversion factor from https://arxiv.org/pdf/1004.3567.pdf
conversion_factor=4*pi*const.epsilon0*(bohr_radius^3);

polz_si=data_table.polarizability_a_u__*conversion_factor;
freq=f2wl(data_table.nm*1e-9);

% add a high freq point so the plot goes off the page
freq= cat(1,1000e12,freq);
polz_si=cat(1,0,polz_si);

plot(freq_axis,freq*1e-12,polz_si*1e38,'k','LineWidth', 2)

%xlim([300,1200])
xlabel('$\omega$, Frequency ($2\pi$ THz)')
ylabel('$\alpha$, Polarizability ($10^{-37} \mathrm{C}\cdot \mathrm{m}^{2} \cdot\mathrm{V}^{-1}$)')
ylim([-1,1]*1.2)
xlim([170,980])
line(xlim(), [0,0], 'LineWidth', 1, 'Color', 'k','LineStyle','--');
xticks(200:200:900)
%set(gca,'Yscale','log')
set(gca,'LineWidth',1.2,'TickLength',[1 1]*0.02);

% % plot the same numbers on x scale
% f_lims=xlim;
% f_lims=f_lims*1e12;
% wl_lims=f2wl(f_lims);
% wl_vals_in_f=freq_axis.XTick*1e12;
% wl_vals=f2wl(wl_vals_in_f)*1e9;
% wl_vals=round(wl_vals);
% [wl_vals_in_f,sort_order]=sort(wl_vals_in_f);
% set(wl_axis,'xlim',f_lims)
% set(wl_axis,'XTick',wl_vals_in_f)
% wl_axis.XTickLabel=arrayfun(@num2str,wl_vals(sort_order),'UniformOutput',false)


% now we plot round numbers on the wl scale
f_lims=xlim;
f_lims=f_lims*1e12;
wl_lims=f2wl(f_lims);
%wl_vals=400:200:1500;
wl_vals=[350,400,500,700,900,1200];
wl_vals_in_f=f2wl(wl_vals*1e-9);
[wl_vals_in_f,sort_order]=sort(wl_vals_in_f);
set(wl_axis,'xlim',f_lims)
set(wl_axis,'XTick',wl_vals_in_f)
wl_axis.XTickLabel=arrayfun(@num2str,wl_vals(sort_order),'UniformOutput',false)

set(wl_axis,'LineWidth',1.2,'TickLength',[1 1]*0.02);


% inset plot
f_to=725736810e6;
inset_axis=axes('Position',[.39 .7 .25 0.25])
set(inset_axis,'FontSize',10);
box(inset_axis,'on')
set(inset_axis,'LineWidth',1.2,'TickLength',[1 1]*0.04);
plot(inset_axis,freq*1e-12-f_to*1e-12,polz_si*1e41,'k','LineWidth', 2)
line(xlim(), [0,0], 'LineWidth', 1, 'Color', 'k','LineStyle','--');
xlim(inset_axis,[-1,1]*1.3)
xlabel(inset_axis,'$\omega-\omega_{\mathrm{TO}}$ ($2\pi$ THz)')
ylabel(inset_axis,'$\alpha$ ($10^{-41} \mathrm{C}\cdot \mathrm{m}^{2} \cdot\mathrm{V}^{-1}$)')
set(gca,'LineWidth',1.2,'TickLength',[1 1]*0.02);

% plot(freq_axis,freq*1e-12,polz_si*1e37,'k','LineWidth', 2)
% 
% %xlim([300,1200])
% xlabel('Frequency (THz)')
% ylabel('Polarizability ($10^{-37} \mathrm{C}\cdot \mathrm{m}^{2} \cdot\mathrm{V}^{-1}$)')
% ylim([-1,1]*1)
% xlim([170,980])
% line(xlim(), [0,0], 'LineWidth', 1, 'Color', 'k','LineStyle','--');
% xticks(200:200:900)
% %set(gca,'Yscale','log')
% set(gca,'LineWidth',1.2,'TickLength',[1 1]*0.02);




%%
