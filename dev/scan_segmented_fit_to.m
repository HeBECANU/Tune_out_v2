function to_seg_fits=scan_segmented_fit_to(anal_opts_fit_to,data)

anal_opts_fit_to.bootstrap = false;

delta_trap_freq_all = data.osc_fit.trap_freq_recons' - data.cal.freq_drift_model(data.mcp_tdc.time_create_write(:,1));

temp_cal=data.mcp_tdc.probe.calibration';
temp_cal(isnan(temp_cal))=1;    
delta_mask = ~isnan(delta_trap_freq_all);
probe_dat_mask = data.osc_fit.ok.rmse & ~isnan(data.wm_log.proc.probe.freq.act.mean)'...
& ~isnan(data.osc_fit.trap_freq_recons) & delta_mask' & temp_cal;

delta_max = max(delta_trap_freq_all(probe_dat_mask));
delta_min = min(delta_trap_freq_all(probe_dat_mask));

to_seg_fits=[];
to_seg_fits.num_shots=sum(probe_dat_mask);
to_seg_fits.fit_mask=probe_dat_mask;
probe_freq= data.wm_log.proc.probe.freq.act.mean(probe_dat_mask')*1e6; %convert to hertz
probe_freq_all= data.wm_log.proc.probe.freq.act.mean*1e6; %convert to hertz
trap_freq=data.osc_fit.trap_freq_recons(probe_dat_mask)';
cal_trap_freq=data.cal.freq_drift_model(data.mcp_tdc.time_create_write(probe_dat_mask,1));
shot_time_abs=data.mcp_tdc.time_create_write(probe_dat_mask,2);

%manual bootstrap rand(size(data.osc_fit.ok.rmse))>0.9

cdat=viridis(1000);
c_cord=linspace(0,1,size(cdat,1));

run_time_start=min(shot_time_abs);
to_seg_fits.time.run_start=run_time_start;
shot_time=shot_time_abs-run_time_start;
shot_time_scaled=shot_time/range(shot_time);
cdat_all=[interp1(c_cord,cdat(:,1),shot_time_scaled),...
    interp1(c_cord,cdat(:,2),shot_time_scaled),...
    interp1(c_cord,cdat(:,3),shot_time_scaled)];




delta_trap_freq=trap_freq-cal_trap_freq;
square_trap_freq= (trap_freq).^2-(cal_trap_freq).^2;
        

setpts_all=data.wm_log.proc.probe.freq.set;
diffs = diff(setpts_all);
[scan_edges,~] = find(diffs<-1e3);
to_idxs = scan_edges+round(0.5*gradient(scan_edges));
to_times = data.mcp_tdc.time_create_write(to_idxs,2)-run_time_start;

setpts=data.wm_log.proc.probe.freq.set(probe_dat_mask);
freq_offset=mean(probe_freq);
num_good_shots = length(setpts);
num_shots = length(data.mcp_tdc.shot_num);
if isempty(scan_edges)
    warning('No full scan completed')
    iimax = 0;
elseif length(scan_edges) ==1 
    warning('Only one full scan completed')
    scan_length = scan_edges(1);
    iimax = 2;
end
iimax = length(scan_edges)+1;
shot_idx_all = 1:num_shots;
%%
%now we do some fitting
%thresholds for CI
%sd         CI
%1          0.3174
%2          0.05
%3          2.699e-03
ci_size_disp=0.3174;%one sd %confidence interval to display
ci_size_cut_outliers=0.01; %confidence interval for cutting outliers

to_seg_fits.fit_all.freq.val=nan(iimax,1);
to_seg_fits.fit_all.freq.unc=nan(iimax,1);
to_seg_fits.fit_all.model=cell(iimax,1);

to_seg_fits.fit_trimmed.freq.val=nan(iimax,1);
to_seg_fits.fit_trimmed.freq.unc=nan(iimax,1);
to_seg_fits.fit_trimmed.model=cell(iimax,1);
to_seg_fits.fit_trimmed.to_unc_boot=nan(iimax,1);

figure()
fprintf('fitting tune out in segments %04u:%04u',iimax,0)
for ii=1:iimax
    %% Convert to time in hours?
    if ii == iimax
        seg_end = num_shots;
        seg_start = scan_edges(iimax-1)+1;
    else
        seg_end = scan_edges(ii);
        if ii==1
            seg_start = 1;
        else
            seg_start = scan_edges(ii-1)+1;
        end
    end
    % Plot three things. ONE: All setpoints by segment CHECK
    seg_mask_temp = [zeros(seg_start-1,1);ones(seg_end-seg_start+1,1);zeros(num_shots-seg_end,1)]'==1;
    
%     set_sel = setpts_all(seg_mask);
%     shot_idx = shot_idx_all(seg_mask_temp);
%     shot_idx_good = shot_idx_all(seg_mask);
    
    
%     to_seg_fits.time.start(ii) = shot_time_abs(seg_start);
%     to_seg_fits.time.stop(ii)  = shot_time_abs(seg_end);
   
    % Two: Set points of good shots by segment CHECK
    seg_mask = seg_mask_temp&probe_dat_mask;
    good_shot_idx = shot_idx_all(seg_mask);
    set_sel = setpts_all(seg_mask);
    seg_size = numel(set_sel);
    delta_sig = delta_trap_freq_all(seg_mask);
    
%     sfigure(60);
    if mod(ii,2) == 0
        marker='x';
    else
        marker ='o';
    end
    subplot(4,4,[1 2])
    pl= plot(shot_idx_all(seg_mask_temp),setpts_all(seg_mask_temp),marker);
    title('Setpoints by segment')
    pl.Color=cdat(ceil(800*ii/iimax),:);
    hold on
    if seg_size>0
        subplot(4,4,[3 4])
        pl=plot(good_shot_idx,set_sel,marker);
        pl.Color=cdat(ceil(800*ii/iimax),:);
        hold on
        title('Setpts of good shots')      
        
        subplot(4,4,[5 6])
        pl=plot(good_shot_idx,delta_sig,marker);
        pl.Color=cdat(ceil(800*ii/iimax),:);
        hold on
        title('Signal of good shots')
    
    
%     % Old code
    fprintf('\b\b\b\b%04u',ii)
%     if sum(seg_mask)>anal_opts_fit_to.min_pts %fits only work when there are a few points
%         set up the data input for the fit
        xdat=probe_freq_all(seg_mask)';
        ydat=delta_sig';
%         cdat=cdat_all(seg_mask);
% 
        xdat=xdat-freq_offset; %fits work better when they are scaled reasonably
        xdat=xdat*anal_opts_fit_to.scale_x;
        xydat=num2cell([xdat;ydat],1);
% 
        [to_freq,mdl_all]=fit_lin_data(xydat);
        to_seg_fits.fit_all.freq.val(ii)=to_freq/anal_opts_fit_to.scale_x+freq_offset;
        to_seg_fits.fit_all.model{ii}=mdl_all;
        mdl_coef=mdl_all.Coefficients;
        to_seg_fits.fit_all.freq.unc(ii)=abs((mdl_coef.Estimate(1)/mdl_coef.Estimate(2)))*...
            sqrt((mdl_coef.SE(1)/mdl_coef.Estimate(1))^2+...
            (mdl_coef.SE(2)/mdl_coef.Estimate(2))^2)/anal_opts_fit_to.scale_x;

        xsamp=linspace(min(xdat),max(xdat),1e3)'; %sample for the model curve
        [ysamp,yci]=predict(mdl_all,xsamp,'Prediction','observation','Alpha',ci_size_disp); %note the observation CI
% %         now make a new prediction with the model but with the CI to cut out outliers
        [~,yci_cull_lim]=predict(mdl_all,xdat','Prediction','observation','Alpha',ci_size_cut_outliers);
        is_outlier_idx=ydat>yci_cull_lim(:,1)' & ydat<yci_cull_lim(:,2)';

        xdat_culled=xdat(is_outlier_idx);
        ydat_culled=ydat(is_outlier_idx);

        culed_xydat=num2cell([xdat_culled;ydat_culled],1);

        [to_freq,mdl_culled]=fit_lin_data(culed_xydat);
        to_seg_fits.fit_trimmed.freq.val(ii)=to_freq/anal_opts_fit_to.scale_x+freq_offset;
        to_seg_fits.fit_trimmed.model{ii}=mdl_culled;
        mdl_coef=mdl_culled.Coefficients;
        to_seg_fits.fit_trimmed.freq.unc(ii)=abs((mdl_coef.Estimate(1)/mdl_coef.Estimate(2)))*...
            sqrt((mdl_coef.SE(1)/mdl_coef.Estimate(1))^2+...
            (mdl_coef.SE(2)/mdl_coef.Estimate(2))^2)/anal_opts_fit_to.scale_x;

        if anal_opts_fit_to.bootstrap
            fprintf('\n')
            boot=bootstrap_se(@fit_lin_data,culed_xydat,...
                'plots',false,...
                'replace',true,...
                'samp_frac_lims',[0.3,0.9],...%[0.005,0.9]
                'num_samp_frac',1,...
                'num_samp_rep',30,...
                'plot_fig_num',89,...
                'save_multi_out',0);
            fprintf('\n%04u',ii)
            to_seg_fits.fit_trimmed.boot{ii}=boot;
            to_seg_fits.fit_trimmed.to_unc_boot(ii)=boot.se_opp_unweighted/anal_opts_fit_to.scale_x;
            to_seg_fits.fit_trimmed.single_shot_uncert_boot(ii)=to_seg_fits.fit_trimmed.to_unc_boot(ii)*sqrt(numel(xdat_culled));
        else
            to_seg_fits.fit_trimmed.boot{ii}={};
            to_seg_fits.fit_trimmed.to_unc_boot(ii)=nan;
            to_seg_fits.fit_trimmed.single_shot_uncert_boot(ii)=nan;
        end
        xsamp_culled=linspace(min(xdat_culled),max(xdat_culled),1e3)';
        [ysamp_culled,yci_culled]=predict(mdl_culled,xsamp_culled,'Alpha',0.2); %'Prediction','observation'
        
        
% %         now plot the remaining data along with the fit model and the model CI
%         if anal_opts_fit_to.plots
%             now plot the data and the model together
            subplot(4,4,7)
            set(gcf,'color','w')
            plot(xsamp,ysamp,'k-')
            hold on
            plot(xsamp,yci,'r-')
            scatter(xdat,ydat,30,'square','filled')
            c =colorbar;
            c.Label.String = 'time (H)';
            caxis([0,range(shot_time)/(60*60)])
            plot(xdat,ydat,'bx')
            xlabel(sprintf('probe beam set freq - %.3f(GHz)',freq_offset*1e-9))
            ylabel('Response (Hz^2)')
            title('Good Data')
            first_plot_lims=[get(gca,'xlim');get(gca,'ylim')];
% %             color the ones that will be removed
            plot(xdat(~is_outlier_idx),ydat(~is_outlier_idx),'r.','markersize',15)
            hold off

            subplot(4,4,8)
            plot(xsamp_culled,ysamp_culled,'k-')
            hold on
            plot(xsamp_culled,yci_culled,'r-')
            scatter(xdat_culled,ydat_culled,30,'square','filled')
            colormap(viridis(1000))
            c =colorbar;
            c.Label.String = 'time (H)';
            caxis([0,range(shot_time)/(60*60)])
            hold off
            xlabel(sprintf('probe beam set freq - %.3f (GHz)',freq_offset*1e-9))
            ylabel('Response (Hz^2)')
            title('Fit Outliers Removed')
            set(gca,'xlim',first_plot_lims(1,:))
            set(gca,'ylim',first_plot_lims(2,:))
            pause(1e-6)
            set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])
            plot_name='TO_fits';
        end



            freq_plot_scaling=1e-9;
            if anal_opts_fit_to.bootstrap
                to_unc=freq_plot_scaling*to_seg_fits.fit_trimmed.to_unc_boot(ii);
            else
                to_unc=freq_plot_scaling*to_seg_fits.fit_trimmed.freq.unc(ii);
            end

            to_seg_freq=to_seg_fits.fit_trimmed.freq.val(ii);
            to_val=freq_plot_scaling*to_seg_freq;
            if ii==1
                to_val_ref = to_val;
            end
            subplot(4,4,[9 10])
            set(gcf,'color','w')
            to_time = to_times(ii)/3600;
            plot(to_time,to_val-to_val_ref,'kx')
            hold on
            plot(to_time,to_val-to_unc-to_val_ref,'b.')
            plot(to_time,to_val+to_unc-to_val_ref,'b.')
%             hold off
            xlabel('time (h)')
            ylabel('variation in TO fit (GHz)')

            subplot(4,4,[13 14])
            set(gcf,'color','w')
            plot(to_time,to_val,'kx')
            hold on
            plot(to_time,(to_val-to_unc),'b.')
            plot(to_time,(to_val+to_unc),'b.')
            xlabel('time (h)')
%             ylabel('variation in TO fit (GHz)')

            subplot(4,4,[11 12 15 16])
%             bins = linspace(-3,3,15);
            histogram((to_seg_fits.fit_trimmed.freq.val(1:ii))-nanmean(to_seg_fits.fit_trimmed.freq.val),15)
            xlabel(sprintf('Scan TO - %.3f (GHz)',1e-9*nanmean(to_seg_fits.fit_trimmed.freq.val)))
            ylabel('Number of scans')

end
    subplot(4,4,[1, 2])
hold off
subplot(4,4,[3,4])
hold off    
subplot(4,4,[5, 6])
hold off
subplot(4,4,7)
hold off
subplot(4,4,8)
hold off

suptitle('Inter-scan variation of TO')
end



% subplot(3,3,12)
% set(gcf,'color','w')
% mean_freq=nanmean(to_val);
% times=to_seg_fits.time.start-to_seg_fits.time.run_start;
% times=times/(60*60);
% plot(times,to_val-mean_freq,'k')
% hold on
% plot(times,to_val-to_unc-mean_freq,'b')
% plot(times,to_val+to_unc-mean_freq,'b')
% hold off
% xlabel('time (h)')
% ylabel('variation in TO fit (GHz)')

% 
% 
% 
% 
% 

% end


function [to_freq,fit_mdl]=fit_lin_data(in)
inmat=cell2mat(in);
xdat=inmat(1,:);
ydat=inmat(2,:);

modelfun = @(b,x) b(1)+ b(2).*x; %simple linear model
opts = statset('nlinfit');
%opts.RobustWgtFun = 'welsch' ; %a bit of robust fitting
%opts.Tune = 1;
beta0 = [1e-5,1e-2]; %intial guesses
fit_mdl = fitnlm(xdat,ydat,modelfun,beta0,'Options',opts);
to_freq=-fit_mdl.Coefficients.Estimate(1)/fit_mdl.Coefficients.Estimate(2);
end