path='./figs/level_diagram_and_polz/fig1b_data.csv';
data_table=readtable(path);

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

polz_au=data_table.polarizability_a_u__*conversion_factor;
freq=f2wl(data_table.nm*1e-9);

% add a high freq point so the plot goes off the page
freq= cat(1,1000e12,freq);
polz_au=cat(1,0,polz_au);

yscale=1e-38;
plot(freq_axis,freq*1e-12,polz_au/yscale,'k','LineWidth', 2)

%xlim([300,1200])
xlabel('$\omega$, Frequency ($2\pi$ THz)')
ylabel(sprintf('$\\alpha^{(s,v,T)}$, Polarizability ($10^{%.0f} \\mathrm{C}\\cdot \\mathrm{m}^{2} \\cdot\\mathrm{V}^{-1}$)',log10(yscale)))
%ylabel('$\alpha$, Polarizability ($10^{-38} \mathrm{C}\cdot \mathrm{m}^{2} \cdot\mathrm{V}^{-1}$)')
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
plot(inset_axis,freq*1e-12-f_to*1e-12,polz_au*1e41,'k','LineWidth', 2)
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


%% make an interpolating model from this data and try fitting to it

freq=f2wl(data_table.nm*1e-9);
au_polz=data_table.polarizability_a_u__;
mask=~isnan(freq) & ~isnan(au_polz) & ~isinf(freq) & ~isinf(au_polz);

freq=freq(mask);
au_polz=au_polz(mask);

polz_model_au=@(x) interp1(freq,au_polz,x,'pchip');


stfig('atomic polz');
clf
fsamp=linspace(1,1000,1e5)*1e12;
polz_au=polz_model_au(fsamp);
plot(fsamp*1e-12,polz_au)
ylim([-200,600])
yline(0)

fminopt = optimset('TolX',1e-20,'TolFun',1e-6);   %,'Display','iter'
tune_out_freq=fzero(@(x) polz_model_au(x),f2wl(413e-9),fminopt);  
au_polz_at_to=aprox_he_polz(tune_out_freq);
f2wl(tune_out_freq)
hold on
plot(tune_out_freq*1e-12,au_polz_at_to,'xr')
hold off
fprintf('found tune out at %.6f THz, %.6f nm, polz %g au \n',tune_out_freq*1e-12,f2wl(tune_out_freq)*1e9,aprox_he_polz(tune_out_freq))

%%
to_polz_au_deriv=arrayfun(@(y) derivest(@(x) aprox_he_polz(x),tune_out_freq,'DerivativeOrder',y),1:4);
fprintf('au polz deriv /hz')
to_polz_au_deriv(1)
to_polz_au_deriv(2)


%%

fsamp=col_vec(linspace(-1,1,1e5))*4e9+tune_out_freq;
polz_au=polz_model_au(fsamp);
stfig('atomic polz from lyt');
clf
xscale=1e-9;
yscale=1e0;
plot((fsamp-tune_out_freq)*xscale,polz_au*yscale)
xlabel('$\omega_{\mathrm{TO}}-\omega$ ($2\pi$ GHz)')
%ylabel(sprintf('$\\alpha$ ($10^{%g} \\mathrm{C}\\cdot \\mathrm{m}^{2} \\cdot\\mathrm{V}^{-1}$)',-log10(yscale)))
ylabel('$\alpha$')
yline(0)


predictor=(fsamp-tune_out_freq)*1e-9;
response=polz_au*yscale;

fit_in=cat(2,predictor,response,response*nan);
meth_lin_fit=fit_poly_with_int(fit_in,1,0,0);
fprintf('lin fit intercept %f MHz \n',(meth_lin_fit.x_intercept.val/xscale)*1e-6)

fit_in=cat(2,predictor,response,response*nan);
meth_lin_fit=fit_poly_with_int(fit_in,2,0,0);
fprintf('quad fit intercept %f MHz \n',(meth_lin_fit.x_intercept.val/xscale)*1e-6)
