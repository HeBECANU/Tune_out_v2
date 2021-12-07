% function dum = spectral_purity_plot()
%Script that scrapes the analysed data from dirs (currently messy but works)
%setup directories you wish to loop over
use_dirs = {
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\unsorted\filt_dep\20190218_filt_dep_0\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\unsorted\filt_dep\20190218_filt_dep_1\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\unsorted\filt_dep\20190218_filt_dep_1_run2\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\unsorted\filt_dep\20190218_filt_dep_2\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\unsorted\filt_dep\20190219_filt_dep_2_run2\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\unsorted\filt_dep\20190219_filt_dep_3_run2\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\unsorted\filt_dep\20190220_filt_dep_1_run3\'
    };


selected_dirs = 1:numel(use_dirs);
shot_idx = 1;
data = cell(numel(use_dirs),1);
filt_num = zeros(numel(use_dirs),1);
for loop_idx=selected_dirs
    loop_config.dir = {use_dirs{loop_idx}};
    current_dir = loop_config.dir{1};
    data{loop_idx} = load_pocessed_to_data(loop_config);
    strt = strfind(loop_config.dir{1},'dep_');
    filt_num(loop_idx) = str2double(current_dir(strt+4));
end
to_vals = cell2mat(cellfun(@(x) x.main.lin.to.val,data,'uni',0));
to_unc = cell2mat(cellfun(@(x) x.main.lin.to.unc',data,'uni',0));

%%

w_fit = 1./to_unc(:,2).^2;
w_fit = w_fit/sum(w_fit);

[mdl,gof,fit_out] = fit(filt_num,to_vals,'poly1');
fit_offset = mdl.p2;
vals =[mdl.p1,mdl.p2];
Alpha = 1-erf(1/sqrt(2));
CI=confint(mdl,Alpha);
PI=predint(mdl,Alpha);

to_means = zeros(4,1);
to_ses = zeros(4,1);
for idx = 1:4
   mask = filt_num == idx - 1;
    to_means(idx) = mean(to_vals(mask));
    to_ses(idx) = vecnorm(to_unc(mask,2))/sum(mask);
end

x_plt = linspace(-0.2,3.2,300);
mdl_y = mdl(x_plt);
mdl_PI = predint(mdl,x_plt,Alpha,'observation','off')/1e6;
mdl_CI = predint(mdl,x_plt,Alpha,'functional','off')/1e6;

fprintf('Offset %.f (%.f,%.f) MHz\n',vals(2)/1e6,1e-6*(CI(:,2)'-vals(2)))
fprintf('Slope %.f (%.f,%.f) MHz\n',vals(1)/1e6,1e-6*(CI(:,1)'-vals(1)))
% fprintf('Slope &.f (%.f,%.f) MHz/filter'

pred_vals=mdl(filt_num);
pred_PI = predint(mdl,filt_num,erf(1/sqrt(2)),'observation','off')-pred_vals;
pred_std = mean(abs(pred_PI),2);
chi_square = sum((fit_out.residuals./pred_std).^2);
chi_square_per_dof = sum((fit_out.residuals./pred_std).^2)/gof.dfe

colors_main= [[233,87,0];[33,188,44];[0,165,166]]./255;
lch=colorspace('RGB->LCH',colors_main(:,:));
lch(:,1)=lch(:,1)+20;
colors_detail=colorspace('LCH->RGB',lch);

stfig('Filt dep');
clf
hold on

fill([x_plt,fliplr(x_plt)],[mdl_PI(:,1);fliplr(mdl_PI(:,2))]'-fit_offset/1e6,...
    0.7*[1,1,1],'FaceAlpha',0.3,'EdgeColor','none')

fill([x_plt,fliplr(x_plt)],[mdl_CI(:,1);fliplr(mdl_CI(:,2))]'-fit_offset/1e6,...
    0.3*[1,1,1],'FaceAlpha',0.3,'EdgeColor','none')
    
plot(x_plt,1e-6*(mdl_y-fit_offset),'k','LineWidth',2)
% errorbar(filt_num,(to_vals-fit_offset)/1e6,to_unc(:,2)/1e6,'kx')
errorbar(0:3,(to_means-fit_offset)/1e6,to_ses/1e6,'bo',...
    'Markersize',8,'MarkerFaceColor',colors_detail(1,:),...
    'MarkerEdgeColor',colors_main(1,:),'Color',colors_main(1,:),...
    'LineWidth',3)
set(gcf,'color','w')
xlabel('Number of filters')
ylabel(sprintf('Tune-out value - %.f (MHz)',fit_offset/1e6))
xlim([-0.1, 3.2])
box on
set(gca,'FontSize',18,'LineWidth',2)
%%
%make a nice figure

to_freqs_val = to_vals/1e6; %convert to blue
to_freqs_err = to_unc(:,2)/1e6; %convert to blue

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

errorbar(filt_num,to_freqs_val,to_freqs_err)


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
% end

