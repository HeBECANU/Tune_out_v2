function to_seg_fits=scan_segmented_fit_to(anal_opts_fit_to,data)
to_seg_fits=[];
% Settings
ci_size_disp=0.3174;%one sd %confidence interval to display
ci_size_cut_outliers=0.1; %confidence interval for cutting outliers
cdat=viridis(1000);
%Import parameters & init
setpts_all=data.wm_log.proc.probe.freq.set;
temp_cal=data.mcp_tdc.probe.calibration';
anal_opts_fit_to.bootstrap = false;
delta_trap_freq_all = data.osc_fit.trap_freq_recons' - data.cal.freq_drift_model(data.mcp_tdc.time_create_write(:,1));
probe_freq= data.wm_log.proc.probe.freq.act.mean*1e6; %convert to hertz
probe_freq_all= data.wm_log.proc.probe.freq.act.mean*1e6; %convert to hertz
shot_time_abs=data.mcp_tdc.time_create_write(:,2);
shot_time=shot_time_abs-min(shot_time_abs);




% trap_freq=data.osc_fit.trap_freq_recons';
% cal_trap_freq=data.cal.freq_drift_model(data.mcp_tdc.time_create_write(:,1));
% square_trap_freq= (trap_freq).^2-(cal_trap_freq).^2;

% Compute useful objects & masks

diffs = diff(setpts_all);
temp_cal(isnan(temp_cal))=1;    
edge_mask = diffs < -1e3;
[to_seg_fits.scan_edges,~] = find(edge_mask);


delta_mask = ~isnan(delta_trap_freq_all);
probe_dat_mask = data.osc_fit.ok.rmse & ~isnan(data.wm_log.proc.probe.freq.act.mean)'...
& ~isnan(data.osc_fit.trap_freq_recons) & delta_mask' & ~temp_cal;
probe_dat_mask = data.osc_fit.ok.rmse & ~isnan(data.wm_log.proc.probe.freq.act.mean)'...
& ~isnan(data.osc_fit.trap_freq_recons) & delta_mask' & ~temp_cal;
last_good_shot = find(probe_dat_mask==1,1,'last');
to_seg_fits.scan_edges = to_seg_fits.scan_edges(to_seg_fits.scan_edges<last_good_shot);
run_time_start=min(shot_time_abs);
num_shots = length(data.mcp_tdc.shot_num);
shot_idx_all = 1:num_shots;
to_idxs = to_seg_fits.scan_edges+round(0.5*gradient(to_seg_fits.scan_edges));
to_times = data.mcp_tdc.time_create_write(to_idxs,2)-run_time_start;

% Mask lists
probe_freq= probe_freq(probe_dat_mask'); %convert to hertz
trap_freq=data.osc_fit.trap_freq_recons(:);
cal_trap_freq=data.cal.freq_drift_model(data.mcp_tdc.time_create_write(:,1));
square_trap_freq= (trap_freq).^2-(cal_trap_freq).^2;
freq_offset=mean(probe_freq);


if isempty(to_seg_fits.scan_edges)
    warning('No full scan completed')
    iimax = 0;
elseif length(to_seg_fits.scan_edges) ==1 
    warning('Only one full scan completed')
    iimax = 2;
else
    iimax = length(to_seg_fits.scan_edges); %Avoids partial scans at the end
end

% Init
to_seg_fits.fit_all.freq.val=nan(iimax,1);
to_seg_fits.fit_all.freq.unc=nan(iimax,1);
to_seg_fits.fit_all.model=cell(iimax,1);

to_seg_fits.fit_trimmed.freq.val=nan(iimax,1);
to_seg_fits.fit_trimmed.freq.unc=nan(iimax,1);
to_seg_fits.fit_trimmed.model=cell(iimax,1);
to_seg_fits.fit_trimmed.to_unc_boot=nan(iimax,1);
to_seg_fits.time.run_start=run_time_start;
to_seg_fits.num_shots=sum(probe_dat_mask);
to_seg_fits.fit_mask=probe_dat_mask;

to_seg_fits.set_sel = cell(iimax,1);
to_seg_fits.delta_sig = cell(iimax,1);
to_seg_fits.xdat=cell(iimax,1);
to_seg_fits.ydat=cell(iimax,1);
to_seg_fits.good_shot_idx=cell(iimax,1);
to_seg_fits.to_time = zeros(iimax,1);
to_seg_fits.seg_edges = zeros(iimax,2);
to_seg_fits.atom_num = zeros(iimax,2);
%% Process
fprintf('fitting tune out in segments %04u:%04u',iimax,0)
for ii=1:iimax 
    fprintf('\b\b\b\b%04u',ii)
    if ii == iimax
        seg_end = num_shots;
        seg_start = to_seg_fits.scan_edges(iimax-1)+1;
    else
        seg_end = to_seg_fits.scan_edges(ii);
        if ii==1
            seg_start = 1;
        else
            seg_start = to_seg_fits.scan_edges(ii-1)+1;
        end
    end
    
    % Select out data for scanwise fits
    seg_mask_temp = [zeros(seg_start-1,1);ones(seg_end-seg_start+1,1);zeros(num_shots-seg_end,1)]'==1;
    seg_mask = seg_mask_temp&probe_dat_mask;
    delta_sig = square_trap_freq(seg_mask);
    xdat=probe_freq_all(seg_mask)';
    ydat=delta_sig';
    xdat=anal_opts_fit_to.scale_x*(xdat-freq_offset); %fits work better when they are scaled reasonably
    xydat=num2cell([xdat;ydat],1);

    
    if size(ydat,2)>anal_opts_fit_to.min_pts 
        % Linear fit
        [to_freq,mdl_all]=fit_lin_data(xydat);
        % Other output stuff
        to_seg_fits.fit_all.model{ii}=mdl_all;
        mdl_coef=mdl_all.Coefficients;
        % The important stuff
        to_seg_fits.fit_all.freq.val(ii)=to_freq/anal_opts_fit_to.scale_x+freq_offset;
        to_seg_fits.fit_all.freq.unc(ii)=abs((mdl_coef.Estimate(1)/mdl_coef.Estimate(2)))*...
            sqrt((mdl_coef.SE(1)/mdl_coef.Estimate(1))^2+...
            (mdl_coef.SE(2)/mdl_coef.Estimate(2))^2)/anal_opts_fit_to.scale_x;

        
        %sample for the model curve
        xsamp=linspace(min(xdat),max(xdat),1e3)'; 
         %note the observation CI
        [ysamp,yci]=predict(mdl_all,xsamp,'Prediction','observation','Alpha',ci_size_disp);
        
        %%  Outlier removal
        % make a new prediction with the model but with the CI to cut out outliers
        [~,yci_cull_lim]=predict(mdl_all,xdat','Prediction','observation','Alpha',ci_size_cut_outliers);
        is_outlier_idx=ydat>yci_cull_lim(:,1)' & ydat<yci_cull_lim(:,2)';

        xdat_culled=xdat(is_outlier_idx);
        ydat_culled=ydat(is_outlier_idx);
        culled_ydat=num2cell([xdat_culled;ydat_culled],1);
        [to_freq,mdl_culled]=fit_lin_data(culled_ydat);
        
        % The refined fit results
        to_seg_fits.fit_trimmed.freq.val(ii)=to_freq/anal_opts_fit_to.scale_x+freq_offset;
        to_seg_fits.fit_trimmed.model{ii}=mdl_culled;
        % Other outputs
        mdl_coef=mdl_culled.Coefficients;
        to_seg_fits.fit_trimmed.freq.unc(ii)=abs((mdl_coef.Estimate(1)/mdl_coef.Estimate(2)))*...
            sqrt((mdl_coef.SE(1)/mdl_coef.Estimate(1))^2+...
            (mdl_coef.SE(2)/mdl_coef.Estimate(2))^2)/anal_opts_fit_to.scale_x;

        %% Bootstrapping
        if anal_opts_fit_to.bootstrap
            fprintf('\n')
            boot=bootstrap_se(@fit_lin_data,culled_ydat,...
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

        %% Plot the segmentwise fits
%             sfigure(601);
%             scat_clr = cdat(ceil(1000*ii/iimax),:);
%             subplot(1,2,1)
%             set(gcf,'color','w')
%             plot(xsamp,ysamp,'k-')
%             hold on
%             plot(xsamp,yci,'r-')
%             scatter(xdat,ydat,30,scat_clr,'square','filled');
%             
%             c =colorbar;
%             c.Label.String = 'time (H)';
%             caxis([0,range(shot_time)/(60*60)])
%             plot(xdat,ydat,'bx')
%             xlabel(sprintf('probe beam set freq - %.3f(GHz)',freq_offset*1e-9))
%             ylabel('Response (Hz^2)')
%             title('Good Data')
%             first_plot_lims=[get(gca,'xlim');get(gca,'ylim')];
% % %             color the ones that will be removed
%             plot(xdat(~is_outlier_idx),ydat(~is_outlier_idx),'r.','markersize',15)
% 
%             subplot(1,2,2)
%             plot(xsamp_culled,ysamp_culled,'k-')
%             hold on
%             plot(xsamp_culled,yci_culled,'r-')
%             scatter(xdat_culled,ydat_culled,30,scat_clr,'square','filled');
%             colormap(viridis(1000))
%             c =colorbar;
%             c.Label.String = 'time (H)';
%             caxis([0,range(shot_time)/(60*60)])
%             xlabel(sprintf('probe beam set freq - %.3f (GHz)',freq_offset*1e-9))
%             ylabel('Response (Hz^2)')
%             title('Fit Outliers Removed')
%             set(gca,'xlim',first_plot_lims(1,:))
%             set(gca,'ylim',first_plot_lims(2,:))
%             pause(1e-6)
%             set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])
%             plot_name='TO_fits';
        
        % Data to output/plot 
        to_seg_fits.seg_edges(ii,:) = [seg_start,seg_end];
        to_seg_fits.to_time(ii) = to_times(ii)/3600;
        to_seg_fits.set_sel{ii} = setpts_all(seg_mask);
        to_seg_fits.delta_sig{ii} = delta_sig;
        to_seg_fits.xdat{ii}=xdat*anal_opts_fit_to.scale_x;
        to_seg_fits.ydat{ii}=num2cell([xdat;ydat],1);
        to_seg_fits.atom_num(ii,:) = [nanmean(data.mcp_tdc.num_counts(seg_mask)),nanstd(data.mcp_tdc.num_counts(seg_mask))];
    end

end
subplot(1,2,1)
hold off
subplot(1,2,2)
hold off


%% Plot
% Could refactor so entire thing is in a loop over ii - needs a few things stored in to_seg_fits
fprintf('Plotting\n',iimax,0)
    for ii=1:iimax %Plots segmented data

        seg_mask_temp = [zeros(to_seg_fits.seg_edges(ii,1)-1,1);ones(to_seg_fits.seg_edges(ii,2)-to_seg_fits.seg_edges(ii,1)+1,1);...
                zeros(num_shots-to_seg_fits.seg_edges(ii,2),1)]'==1;
        seg_mask = seg_mask_temp&probe_dat_mask;

        sfigure(602);
        if mod(ii,2) == 0
            marker='x';
        else
            marker ='o';
        end

        if sum(seg_mask)>0
%             subplot(4,4,[1 2])
%             pl= plot(shot_idx_all(seg_mask_temp),setpts_all(seg_mask_temp),marker);
%             title('Setpoints by segment')
%             pl.Color=cdat(ceil(800*ii/iimax),:);
%             hold on
            
            subplot(4,4,[1 2])
            pl=plot(shot_time(seg_mask),to_seg_fits.set_sel{ii},marker);
            pl.Color=cdat(ceil(800*ii/iimax),:);
            hold on
            title('Setpts of good shots')
            xlabel('Time (h)')

            subplot(4,4,[3 4])
            pl=plot(shot_time(seg_mask),to_seg_fits.delta_sig{ii},marker);
            pl.Color=cdat(ceil(800*ii/iimax),:);
            hold on
            title('Signal of good shots')
            xlabel('Time (h)')
        end
    end
    
    % Plot everything that was stored
    
    to_time = to_seg_fits.to_time;
    to_val_scale = 1e-9;
    to_val = to_seg_fits.fit_all.freq.val*to_val_scale;
    to_val_trim = to_seg_fits.fit_trimmed.freq.val*to_val_scale;
    to_unc = to_seg_fits.fit_all.freq.unc*to_val_scale;
    to_unc_trim = to_seg_fits.fit_trimmed.freq.unc*to_val_scale;
    to_val_ref = to_val(1);
    to_val_ref_trim = to_val_trim(1);

    
    slopes = cellfun(@(x) x.Coefficients.Estimate(2),to_seg_fits.fit_all.model);
    slopes_err = cellfun(@(x) x.Coefficients.SE(2),to_seg_fits.fit_all.model);
    num_val = to_seg_fits.atom_num(:,1);
    num_unc = to_seg_fits.atom_num(:,2);
    
    to_diffs = abs(diff(to_val));
    to_diffs_trim = abs(diff(to_val_trim));
    
    [acf,lags,bound] = autocorr(to_val);
    [acf_trim,lags_trim,bound_trim] = autocorr(to_val_trim);
    
    subplot(4,4,[5 6])
    set(gcf,'color','w')
    plot(to_time,to_val-to_val_ref,'k')
    hold on
    plot(to_time,to_val-to_unc-to_val_ref,'b.-')
    plot(to_time,to_val+to_unc-to_val_ref,'b.-')
    xlabel('time (h)')
    ylabel('variation in TO fit (GHz)')
    title('TO drift')

    subplot(4,4,[7 8])
    plot(to_time,to_val_trim-to_val_ref_trim,'k')
    hold on
    plot(to_time,to_val_trim-to_unc_trim-to_val_ref_trim,'b.-')
    plot(to_time,to_val_trim+to_unc_trim-to_val_ref_trim,'b.-')
    xlabel('time (h)')
    ylabel('variation in TO fit (GHz)')
    title('Trimmed TO drift')
    
    
    subplot(4,4,[9])
    plot(lags*mean(diff(to_time)),acf,'ro')
    hold on
    plot([min(lags*mean(diff(to_time))),max(lags*mean(diff(to_time)))],[0,0],'k');
    title('TO autocorrelation')
    xlabel('Time (h)')
    ylabel('Autocorrelation')
   
    subplot(4,4,[10])
    plot(lags_trim*mean(diff(to_time)),acf_trim,'ro')
    hold on
    plot([min(lags*mean(diff(to_time))),max(lags*mean(diff(to_time)))],[0,0],'k');
    title('Trim TO autocorrelation')
    xlabel('Time (h)')
    ylabel('Autocorrelation')
    
    subplot(4,4,[11 12])
    plot(shot_time(1:last_good_shot),data.mcp_tdc.num_counts(1:last_good_shot))
    xlabel('Time (h)')
    ylabel('Num counts')
    title('Number trend')

    
    subplot(4,4,[13])
    histogram(to_diffs,10)
    title('\delta TO')
    xlabel('TO(t+1)-TO(t)')
    ylabel('counts')
    
    subplot(4,4,[14])
    errorbar(to_time,slopes,slopes_err)
    title('Fit gradient')
    xlabel('time (h)')
    ylabel('Gradient (a.u.)')

    subplot(4,4,[15])
    num_TO = normalize(num_val).*normalize(to_val);
    num_TO_mean = nanmean(num_TO);
    histogram(num_TO,8,'Normalization','pdf')
    hold on
    plot([num_TO_mean,num_TO_mean],[0,max(histcounts(num_TO,8,'Normalization','pdf'))],'r')
    ylabel('P(X)')
    title('Num counts * TO (normalized)')
    
    subplot(4,4,[16])
    slope_TO = normalize(slopes).*normalize(to_val);
    slope_TO_mean = nanmean(slope_TO);
    histogram(slope_TO,8,'Normalization','pdf')
    hold on
    plot([slope_TO_mean,slope_TO_mean],[0,max(histcounts(slope_TO,8,'Normalization','pdf'))],'r')
    ylabel('P(X)')
    title('fit slope * TO (normalized)')
    

end


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