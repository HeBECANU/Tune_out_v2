function out=fit_temp_segmented(data,anal_opts)
%anal_opts.xylim=anal_opts.tdc_import.txylim(2:3,:);

% lets find the temperature of the BEC usign the AL pulse data
% unfortunately the thermal cloud will not be in phase with the BEC part of the AL pulse
% this is bc the thermal and BEC component have different trap frequencies
% and bc the gaussian fit width using a simple masking procedure is of a similar size ~10mm/s to the oscillation amplitude
% ~10mm/s


out=[];
% this will measure the termperature by taking a cylinder of counts in the weak axis from every atom laser pulse
cen_type='mean';  %could be mean,window or thermal 
%cen_type='seg_low';
%cen_type='seg_low_close';
%cen_type='seg_high';
blur_radius=1e-3;


%data_type='all pulse';
data_type='seg low';

global const

if strcmp(cen_type,'mean')
    vel_center=data.mcp_tdc.al_pulses.vel_zxy.mean; %file,pulse_num,axis) 
elseif strcmp(cen_type,'seg_high')
    vel_center=data.al_segmentation.segmented_vel_xyz.high_flux.mean; %file,pulse_num,axis) 
elseif strcmp(cen_type,'seg_low')
    vel_center=data.al_segmentation.segmented_vel_xyz.low_flux.mean;
elseif strcmp(cen_type,'seg_low_close')
    vel_center=data.al_segmentation.segmented_vel_xyz.close_low_flux.mean;
end    


mask_opts=[];
mask_opts.type='allow cyl';
mask_opts.radius=10e-3/anal_opts.global.fall_time;
mask_opts.num_norm_pts=5e7;

% save all the data that we get from each pulse
windowed_data=cell(size(data.mcp_tdc.al_pulses.pos.mean(:,:,1)));
num_shots=size(data.mcp_tdc.al_pulses.pos.mean,1);
num_pulses=size(data.mcp_tdc.al_pulses.pos.mean,2);

for shot=1:num_shots
    fprintf('%u\n',shot)
    for pulse=1:num_pulses
        %mean(window_data_this,1)
        cen_vel_this=squeeze(vel_center(shot,pulse,:));
        
        if strcmp(data_type,'all pulse')
            counts_vxyz_this=data.al_segmentation.unsegmented_vel_xyz.counts{shot,pulse};
        elseif strcmp(data_type,'seg low')
            counts_vxyz_this=data.al_segmentation.segmented_vel_xyz.low_flux.counts{shot,pulse};
        end
        
        num_counts_this=size(counts_vxyz_this,1);
        if num_counts_this>0
            cen_vel_this=repmat(permute(cen_vel_this,[2,1]),num_counts_this,1);
            if size(counts_vxyz_this)~=size(cen_vel_this)
                error('should be the same size')
            end
            counts_vxyz_this=counts_vxyz_this-cen_vel_this;
             %mean(window_data_this,1)
            %fprintf('%u\n',size(window_data_this,1))

            if strcmp(mask_opts.type,'allow cyl')
                radial_xz=rms(counts_vxyz_this(:,[1,3]),2);
                %min(radial_yt)
                mask=radial_xz<mask_opts.radius;
                vzxy_out=counts_vxyz_this(mask,:);
                windowed_data{shot,pulse}=vzxy_out;
            end
        end
    end
end


norm_counts=(rand(mask_opts.num_norm_pts,3)-0.5)*2;
random_range=40e-3/anal_opts.global.fall_time;
norm_counts=norm_counts*random_range;
if strcmp(mask_opts.type,'allow cyl')
    radial_xz=rms(norm_counts(:,[1,3]),2);
    %min(radial_yt)
    mask=radial_xz<mask_opts.radius;
    norm_windowed=norm_counts(mask,:);
end

   

%% combine all data
all_windows_comb=cat(1,windowed_data{1,:});
radius= vecnorm(all_windows_comb,2,2);
radius_norm= vecnorm(norm_windowed,2,2);



stfig('themal tails');
clf
sigma=5e-4;
sdat=smooth_hist(radius,'sigma',blur_radius);
s_ref=smooth_hist(radius_norm,'sigma',blur_radius,'edges',sdat.bin.edge);

% bin_volume=(4/3)*pi* (sdat.bin.edge(2:end).^3  - sdat.bin.edge(1:end-1).^3) ;
% flux=./bin_volume;

flux=sdat.count_rate.smooth./s_ref.count_rate.smooth_prob;
%flux=sdat.count_rate.smooth
%flux=s_ref.count_rate.smooth_prob


%x_mask_lims=[5e-3,15e-3]./anal_opts.global.fall_time;
x_mask_lims=[14e-3,40e-3]./anal_opts.global.fall_time;

xmask= sdat.bin.centers>x_mask_lims(1)  & sdat.bin.centers<x_mask_lims(2);
xdat=sdat.bin.centers(xmask);
ydat=flux(xmask) *1e-6;

plot(xdat,ydat )
%set(gca,'Yscale','log')
%set(gca,'Xscale','log')

opt = statset('TolFun',1e-10,'TolX',1e-10,...
    'MaxIter',1e4,... %1e4
    'UseParallel',1);
% cof_names={'sigma','amp'};
% modelfun=@(b,x) gaussian_function_1d(x,b(1),0,b(2),0);
% beta0=[0.02,1];

cof_names={'sigma','amp','offset'};
modelfun=@(b,x) gaussian_function_1d(x,b(1),0,b(2),b(3));
beta0=[0.02,1,0];

% cof_names={'sigma','amp','k4amp'};
% modelfun=@(b,x) gaussian_function_1d(x,b(1),0,b(2))+b(3)*x.^-4
% beta0=[0.016,0.01,1e-9];

% cof_names={'sigma','amp','k4amp','offset'};
% modelfun=@(b,x) gaussian_function_1d(x,b(1),0,b(2))+b(3)*x.^-4 +b(4)
% beta0=[0.016,0.01,1e-9,1e-3];

fitobject=fitnlm(xdat,ydat,modelfun,beta0,...
    'options',opt,...
    'CoefficientNames',cof_names);

xsamp=linspace(x_mask_lims(1),x_mask_lims(2),1e3)';
[prediction,ci]=predict(fitobject,xsamp,'Alpha',1-erf(1/sqrt(2))); %,'Prediction','observation'
  hold on
  
color_shaded=[1,1,1]*0.9;
colors_main=[0.5,0.6,0.2];
patch([xsamp', fliplr(xsamp')], [ci(:,1)', fliplr(ci(:,2)')], color_shaded,'EdgeColor','none')  %[1,1,1]*0.80
              
plot(xsamp(:,1),prediction,'-','LineWidth',1.0,'Color',colors_main)
hold off
set(gca,'Xscale','log')
set(gca,'Yscale','log')

fit_temperature=(fitobject.Coefficients.Estimate(1).^2) *const.mhe /(const.kb*3);
fit_temperature


%%
all_windows_comb=cat(1,windowed_data{:,1});
% % simple y axis
% radius= all_windows_comb(:,2);
% radius_norm= norm_windowed(:,2);
% signed radius
radius= vecnorm(all_windows_comb,2,2).*sign(all_windows_comb(:,2));
radius_norm= vecnorm(norm_windowed,2,2).*sign(norm_windowed(:,2));




stfig('themal tails');
clf
sdat=smooth_hist(radius,'sigma',blur_radius);
s_ref=smooth_hist(radius_norm,'sigma',blur_radius,'edges',sdat.bin.edge);

% bin_volume=(4/3)*pi* (sdat.bin.edge(2:end).^3  - sdat.bin.edge(1:end-1).^3) ;
% flux=./bin_volume;

flux=sdat.count_rate.smooth./s_ref.count_rate.smooth_prob;
%flux=sdat.count_rate.smooth
%flux=s_ref.count_rate.smooth_prob


%x_mask_lims=[5e-3,15e-3]./anal_opts.global.fall_time;
%x_mask_lims=[-30e-3,30e-3]./anal_opts.global.fall_time;
x_mask_lims=[-32e-3,32e-3]./anal_opts.global.fall_time;
xmask= sdat.bin.centers>x_mask_lims(1)  & sdat.bin.centers<x_mask_lims(2);

x_mask_lims_inner=[-1,1]*16e-3./anal_opts.global.fall_time;
xmask=xmask & ~(sdat.bin.centers>x_mask_lims_inner(1)  & sdat.bin.centers<x_mask_lims_inner(2) );
% sum(xmask)
xdat=sdat.bin.centers(xmask);
ydat=flux(xmask) *1e-6;

plot(xdat,ydat )
%set(gca,'Yscale','log')
%set(gca,'Xscale','log')

%


opt = statset('TolFun',1e-10,'TolX',1e-10,...
    'MaxIter',1e4,... %1e4
    'UseParallel',1);

cof_names={'sigma','amp'};
modelfun=@(b,x) gaussian_function_1d(x,b(1),0,b(2))
beta0=[0.02,1];

% cof_names={'sigma','amp','offset'};
% modelfun=@(b,x) gaussian_function_1d(x,b(1),0,b(2),b(3));
% beta0=[0.02,1,0];

% cof_names={'sigma','amp','k4amp','offset'};
% modelfun=@(b,x) gaussian_function_1d(x,b(1),0,b(2))+b(3)*x.^-4 +b(4)
% beta0=[0.016,0.01,1e-9,1e-3];



fitobject=fitnlm(xdat,ydat,modelfun,beta0,...
    'options',opt,...
    'CoefficientNames',cof_names);

xsamp=linspace(x_mask_lims(1),x_mask_lims(2),1e3)';
[prediction,ci]=predict(fitobject,xsamp,'Alpha',1-erf(1/sqrt(2))); %,'Prediction','observation'
  hold on
  
color_shaded=[1,1,1]*0.9;
colors_main=[0.5,0.6,0.2];
patch([xsamp', fliplr(xsamp')], [ci(:,1)', fliplr(ci(:,2)')], color_shaded,'EdgeColor','none')  %[1,1,1]*0.80
              
plot(xsamp(:,1),prediction,'-','LineWidth',1.0,'Color',colors_main)
hold off
%set(gca,'Xscale','log')
set(gca,'Yscale','log')

fit_temperature=[];
fit_temperature.val=(fitobject.Coefficients.Estimate(1).^2) *const.mhe /(const.kb*3);
fit_temperature.unc=fit_temperature.val*(fitobject.Coefficients.SE(1)./fitobject.Coefficients.Estimate(1))*2;
fit_temperature



bec_det=bec_properties([420,420,50],2e6);
condensate_fraction=[];
condensate_fraction.val=(1-(fit_temperature.val/bec_det.tc.finite_interacting)^3);
condensate_fraction.unc=(1-condensate_fraction.val)*(fit_temperature.unc/fit_temperature.val)*3;
%condensate_fraction.unc_range=(1-((fit_temperature.val+[-1,1]*fit_temperature.unc)./bec_det.tc.finite_interacting).^3)-condensate_fraction.val


end


%%
% % for reference what a gaussian looks like
% stfig('gauusian analytical')
% xsamp=linspace(1,5,1e3);
% plot(xsamp,gaussian_function_1d(xsamp,1,0))
% set(gca,'Yscale','log')
% set(gca,'Xscale','log')
