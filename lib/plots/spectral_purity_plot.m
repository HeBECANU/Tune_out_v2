function dum = spectral_purity_plot()
%Script that scrapes the analysed data from dirs (currently messy but works)
%setup directories you wish to loop over
loop_config.dir = {
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190218_filt_dep_0\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190218_filt_dep_1\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190218_filt_dep_1_run2\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190218_filt_dep_2\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190218_filt_dep_3\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190219_filt_dep_2_run2\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190219_filt_dep_3_run2\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190220_filt_dep_1_run3\'
    };

data = to_data(loop_config);
selected_dirs = 1:numel(loop_config.dir);
filt_num = zeros(numel(data.drift.to_time),1);
shot_idx = 1;

for loop_idx=selected_dirs
    current_dir = loop_config.dir{loop_idx};
    strt = strfind(current_dir,'dep_');
    filt_num(shot_idx:(shot_idx+data.main.scan_num(loop_idx)-1),1) = ones(data.main.scan_num(loop_idx),1).*str2double(current_dir(strt+4));
    shot_idx = shot_idx+data.main.scan_num(loop_idx);
end
%%
sfigure(3000);
plot_offset = sum(data.drift.to_val{1}.*data.drift.to_val{2})/sum(data.drift.to_val{2});
errorbar(filt_num,(data.drift.to_val{1}-plot_offset).*1e-6,data.drift.to_val{2}.*1e-6,'kx')
set(gcf,'color','w')
xlabel('Filter Number')
ylabel(sprintf('Tune-out value - %.3f (MHz)',plot_offset.*1e-6))
xlim([-0.1, 3.2])

%%
%make a nice figure

to_freqs_val = data.drift.to_val{1}./1e6; %convert to blue
to_freqs_err = data.drift.to_val{2}./1e6; %convert to blue

%fit sin waves to the two sets of data
mdl_fun = @(b,x) b(1)+b(2).*x(:,1);
beta0 = [nanmean(to_freqs_val),1e5];
wlin=1./(to_freqs_err.^2);
opts = statset('MaxIter',1e4);
fit_mdl_lin = fitnlm(filt_num,to_freqs_val,mdl_fun,beta0,...
    'Options',opts,'Weights',wlin,'CoefficientNames' ,{'offset','grad'});

sfigure(8320);
disp_config.colors_main = [[75,151,201];[193,114,66];[87,157,95]];
disp_config.plot_title = '';
disp_config.x_label = 'Number of Filters';
disp_config.x_ticks= 0:3
disp_config.font_name = 'cmr10';
disp_config.font_size_global=14;
disp_config.mdl_fun = mdl_fun;
disp_config.beta0 = beta0;
disp_config.opts=statset('nlinfit');
disp_config.fig_number=3400;
disp_config.bin_tol=0.01;
%set offset to zero filter value
% disp_config.plot_offset.val=fit_mdl_lin.Coefficients.Estimate(1);
% disp_config.plot_offset.unc=fit_mdl_lin.Coefficients.SE(1);
%set offset to value found from lin pol dependence
disp_config.plot_offset.val=predict(fit_mdl_lin,3);

plot_sexy(disp_config,filt_num,to_freqs_val,wlin,fit_mdl_lin)


% errorbar(x_grouped,(y_grouped(:,1)-plot_offset).*1e-6,yneg,ypos,'o','CapSize',0,'MarkerSize',5,'Color',colors_main(1,:),...
%     'MarkerFaceColor',colors_detail(1,:),'LineWidth',1.5);
% xlabel('Filter Number')
% ylabel(sprintf('Tune-out value - %.3f (MHz)',plot_offset.*1e-6))
% set(gca,'xlim',[-0.1,3.1])
% %set(gca,'ylim',first_plot_lims(2,:))
% set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])
% set(gcf,'color','w')
% set(gca,'FontSize',font_size_global,'FontName',font_name)
% box on
end

