function to_seg_fits=scan_segmented_fit_to(anal_opts_fit_to_seg,data)
%this function plots the tune out value for every scan of wavelength(~40 shots) across a data run to try to diagnose the
%drifts that we have observed.

% limitations
% wont work if do bi directional scaning
% does not use last segment

%impovements
% make time as time during probe not start shot

% Settings
ci_size_disp=1-erf(anal_opts_fit_to_seg.sigma_disp/sqrt(2));%one sd %confidence interval to display
ci_size_cut_outliers=1-erf(anal_opts_fit_to_seg.sigma_cut_outliers/sqrt(2)); %confidence interval for cutting outliers
cdat=viridis(1000);
anal_opts_fit_to_seg.bootstrap = false;

%scale the freq values to make the fit easier
scale_freq=anal_opts_fit_to_seg.scale_x;

%initalize
to_seg_fits=[];
%import the data to split up the scans
setpts_all=data.blue_probe.set;
mean_setpt=mean(setpts_all);
is_cal=col_vec(data.mcp_tdc.probe.calibration);
is_cal(isnan(is_cal))=0;    
probe_freq_all= data.blue_probe.act.mean;

shot_time_abs=data.mcp_tdc.time_create_write(:,2);
first_shot_time=shot_time_abs(1);
shot_time_rel=shot_time_abs-first_shot_time;

% split up the scans
set_pt_diffs = diff(setpts_all);
%find the majority gradient
temp_tol_diff=1; %1hz difference
step_up_mask=set_pt_diffs>temp_tol_diff;
step_down_mask=set_pt_diffs<-temp_tol_diff;
no_step_mask=set_pt_diffs==0;
sweep_down_or_up=sum(step_down_mask)<sum(step_up_mask);
%now use whatever is in the minority as the edge mask
if sweep_down_or_up
    edge_mask = step_down_mask;
else
    edge_mask = step_up_mask;
end
%this approach will only capture full scans of the TO
[to_seg_fits.scan_edges,~] = find(edge_mask);
% if to_seg_fits.scan_edges(1)>2 %why this???
%     to_seg_fits.scan_edges=cat(1,1,to_seg_fits.scan_edges);
% end
probe_dat_mask=col_vec(data.signal.dat_mask);
square_trap_freq_val= data.signal.square_probe_trap_freq.val;
square_trap_freq_unc=data.signal.square_probe_trap_freq.unc;
last_good_shot = find(probe_dat_mask==1,1,'last');

%exclude edges if they are past the last good shor 
to_seg_fits.scan_edges = to_seg_fits.scan_edges(to_seg_fits.scan_edges<=last_good_shot);
% move the fist and last (partial) scan data into the second and second last scan
to_seg_fits.scan_edges(1)=1;
to_seg_fits.scan_edges(end)=last_good_shot;

run_time_start=min(shot_time_abs);
num_shots = length(data.mcp_tdc.shot_num);
%shot_idx_all = 1:num_shots;

%get the mean time of the segment by averaging the edge times
time_seg_edges=data.mcp_tdc.time_create_write(to_seg_fits.scan_edges,2);
to_times = ((time_seg_edges(2:end)+time_seg_edges(1:end-1))/2); 
%alt. could use the time at the start of the seg...

% Mask lists
probe_freq= probe_freq_all(probe_dat_mask');
probe_freq_offset=mean(probe_freq);


if isempty(to_seg_fits.scan_edges)
    error('No full scan completed')
    iimax = 0;
elseif length(to_seg_fits.scan_edges) ==1 
    warning('Only one full scan completed')
    iimax = 2;
else
    iimax = length(to_seg_fits.scan_edges)-1; %Avoids partial scans at the end
end

% Init
to_seg_fits.fit_all.to_freq.val=nan(iimax,1);
to_seg_fits.fit_all.to_freq.unc=nan(iimax,1);
to_seg_fits.fit_all.model=cell(iimax,1);

to_seg_fits.fit_trimmed.to_freq.val=nan(iimax,1);
to_seg_fits.fit_trimmed.to_freq.unc=nan(iimax,1);
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
to_seg_fits.avg_osc_coefs_probe = cell(iimax,1);
to_seg_fits.avg_osc_coefs_cal = cell(iimax,1);

is_measument_not_outlier=nan(numel(probe_dat_mask),1);


%% Process
fprintf('fitting tune out in segments %04u:%04u',iimax,0)
for ii=1:iimax 
    fprintf('\b\b\b\b%04u',ii)
    seg_start= to_seg_fits.scan_edges(ii)+1;
    seg_end = to_seg_fits.scan_edges(ii+1);
    %this will miss partial scans at the end
    %if ii == iimax
    %    seg_end = num_shots;
    %    seg_start = to_seg_fits.scan_edges(iimax-1)+1;
    %else
  
    % Select out data for scanwise fits
    seg_mask_temp=col_vec(1:num_shots);
    seg_mask_temp=(seg_mask_temp>=seg_start & seg_mask_temp<seg_end);
    %%seg_mask_temp = [zeros(seg_start-1,1);ones(seg_end-seg_start+1,1);zeros(num_shots-seg_end,1)]'==1;
    
    seg_mask = seg_mask_temp&probe_dat_mask;
    seg_mask_cal = seg_mask_temp& is_cal;
    delta_sig = square_trap_freq_val(seg_mask);
    yunc=square_trap_freq_unc(seg_mask);
    wdat=1./yunc;
    wdat(isnan(wdat))=1e-20; %if nan then set to smallest value you can
    wdat=wdat./sum(wdat);
    xdat=probe_freq_all(seg_mask);
    ydat=delta_sig;
    xdat=scale_freq*(xdat-probe_freq_offset); %fits work better when they are scaled reasonably
    xywdat=num2cell([xdat,ydat,wdat],2);

    
    if numel(ydat)>anal_opts_fit_to_seg.min_pts 
        % Linear fit
        fit_out=fit_poly_with_int(xywdat,1,1,0);
        mdl_all=fit_out.fit_mdl;
        to_freq=fit_out.x_intercept.val;
        % Other output stuff
        to_seg_fits.fit_all.model{ii}=mdl_all;

        % The important stuff
        to_seg_fits.fit_all.model{ii}=fit_out.fit_mdl;
        to_seg_fits.fit_all.scale_x{ii}=scale_freq; %save the scale factor with the model
        to_seg_fits.fit_all.to_freq.val(ii)=fit_out.x_intercept.val/scale_freq+probe_freq_offset; %give the result back in hz
        to_seg_fits.fit_all.to_freq.unc(ii)=fit_out.x_intercept.unc.with_cov/scale_freq;%give the unc back in hz
        
        %sample for the model curve
        %xsamp=linspace(min(xdat),max(xdat),1e3)'; 
         %note the observation CI
        %[ysamp,yci]=predict(mdl_all,xsamp,'Prediction','observation','Alpha',ci_size_disp);
        
        %%  Outlier removal
        % make a new prediction with the model but with the CI to cut out outliers
        [~,yci_cull_lim]=predict(mdl_all,xdat,'Prediction','observation','Alpha',ci_size_cut_outliers);
        is_not_outlier_mask=logical(ydat>yci_cull_lim(:,1) & ydat<yci_cull_lim(:,2));
        culed_xywdat=xywdat(is_not_outlier_mask);
        cdat_culled=cdat(is_not_outlier_mask,:);
        xydatmat=cell2mat(culed_xywdat);
        xdat_culled= xydatmat(:,1);
        ydat_culled= xydatmat(:,2);
        wdat_culled= xydatmat(:,3);
        yunc_culled=yunc(is_not_outlier_mask);
        
        fit_out=fit_poly_with_int(culed_xywdat,1,0,1);
        mdl_culled=fit_out.fit_mdl;
        to_seg_fits.fit_trimmed.model{ii}=fit_out.fit_mdl;
        to_seg_fits.fit_trimmed.scale_x{ii}=scale_freq; %save the scale factor with the model
        to_seg_fits.fit_trimmed.to_freq.val(ii)=fit_out.x_intercept.val/scale_freq+probe_freq_offset; %give the result back in hz
        to_seg_fits.fit_trimmed.to_freq.unc(ii)=fit_out.x_intercept.unc.with_cov/scale_freq;%give the unc back in hz
        
        %to_freq=fit_out.x_intercept.val;
        
        % The refined fit results
        %to_seg_fits.fit_trimmed.to_freq.val(ii)=to_freq/scale_freq+opt_freq_offset;
        %to_seg_fits.fit_trimmed.model{ii}=mdl_culled;
        % Other outputs
        %mdl_coef=mdl_culled.Coefficients;
        %to_seg_fits.fit_trimmed.to_freq.unc(ii)=abs((mdl_coef.Estimate(1)/mdl_coef.Estimate(2)))*...
        %   sqrt((mdl_coef.SE(1)/mdl_coef.Estimate(1))^2+...
        %   (mdl_coef.SE(2)/mdl_coef.Estimate(2))^2)/scale_freq;

        %% Bootstrapping
        if anal_opts_fit_to_seg.bootstrap
            fprintf('\n')
            boot_fun=@(x) wrap_fit_for_boot(x,1,scale_freq);
            boot=bootstrap_se(boot_fun,culed_xywdat,...
                'plots',false,...
                'replace',true,...
                'samp_frac_lims',[0.3,0.9],...%[0.005,0.9]
                'num_samp_frac',1,...
                'num_samp_rep',30,...
                'plot_fig_name','segment bootstrap',...
                'save_multi_out',0);
            fprintf('\n%04u',ii)
            to_seg_fits.fit_trimmed.boot{ii}=boot;
            to_seg_fits.fit_trimmed.to_unc_boot(ii)=boot.se_opp_unweighted/scale_freq;
            to_seg_fits.fit_trimmed.single_shot_uncert_boot(ii)=to_seg_fits.fit_trimmed.to_unc_boot(ii)*sqrt(numel(xdat_culled));
        else
            to_seg_fits.fit_trimmed.boot{ii}={};
            to_seg_fits.fit_trimmed.to_unc_boot(ii)=nan;
            to_seg_fits.fit_trimmed.single_shot_uncert_boot(ii)=nan;
        end
        
        %xsamp_culled=linspace(min(xdat_culled),max(xdat_culled),1e3)';
        %[ysamp_culled,yci_culled]=predict(mdl_culled,xsamp_culled,'Alpha',0.2); %'Prediction','observation'

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
        is_measument_not_outlier(seg_mask)=is_not_outlier_mask;
        to_seg_fits.seg_edges(ii,:) = [seg_start,seg_end];
        to_seg_fits.to_time(ii) = to_times(ii);
        to_seg_fits.set_sel{ii} = setpts_all(seg_mask);
        to_seg_fits.delta_sig{ii} = delta_sig;
        to_seg_fits.xdat{ii}=xdat*scale_freq;
        to_seg_fits.ydat{ii}=num2cell([xdat;ydat],1); %TODO: FIX SHOULD JUST BE YDAT
        to_seg_fits.atom_num(ii,:) = [nanmean(data.mcp_tdc.num_counts(seg_mask)),nanstd(data.mcp_tdc.num_counts(seg_mask))];
        to_seg_fits.avg_osc_coefs_probe{ii} = squeeze(mean(data.osc_fit.model_coefs(seg_mask,:,:))); %Save the average fit coefficents for a scan segment also
        to_seg_fits.avg_osc_coefs_cal{ii} = squeeze(mean(data.osc_fit.model_coefs(seg_mask_cal,:,:)));
    end %if min number of points in scan
end %for each scan
fprintf('\n')
outlier_mask=is_measument_not_outlier;
outlier_mask(~isnan(outlier_mask))=~outlier_mask(~isnan(outlier_mask));
outlier_mask(isnan(outlier_mask))=false;
outlier_mask=logical(outlier_mask); %nans removed can be logical

%% mean TO val
% find the weighted mean TO freq(and the uncert) across scans

good_scan_mask=~isnan(to_seg_fits.fit_trimmed.to_freq.val) & ~isnan(to_seg_fits.fit_trimmed.to_freq.unc);
seg_unc=to_seg_fits.fit_trimmed.to_freq.unc(good_scan_mask);
seg_weights=seg_unc.^(-2);
to_seg_fits.fit_trimmed.to_freq.sewm=[];
[to_seg_fits.fit_trimmed.to_freq.sewm.unc,to_seg_fits.fit_trimmed.to_freq.sewm.val]=sewm(to_seg_fits.fit_trimmed.to_freq.val(good_scan_mask),seg_weights);


fprintf('total number of outliers %u\n',nansum(outlier_mask))
%% Plot
% Could refactor so entire thing is in a loop over ii - needs a few things stored in to_seg_fits
fprintf('Plotting\n',iimax,0)
stfig('TO segments','add_stack',1);

%clear the plot
if anal_opts_fit_to_seg.clear_plot
    clf
end
dead_rows = [];
for ii=1:iimax %Plots segmented data
    if ~(to_seg_fits.seg_edges(ii,2)>to_seg_fits.seg_edges(ii,1))
        %remove the dead row
        dead_rows = [dead_rows, ii];
        continue
    end
%   seg_mask_temp = [zeros(to_seg_fits.seg_edges(ii,1)-1,1);ones(to_seg_fits.seg_edges(ii,2)-to_seg_fits.seg_edges(ii,1)+1,1);...
%             zeros(num_shots-to_seg_fits.seg_edges(ii,2),1)]'==1;
%   isequal(col_vec(seg_mask_temp),seg_mask_temp2)
    seg_start= to_seg_fits.scan_edges(ii)+1;
    seg_end = to_seg_fits.scan_edges(ii+1);     
    
    seg_mask_temp=col_vec(1:num_shots);
    seg_mask_temp=(seg_mask_temp>=seg_start & seg_mask_temp<seg_end);
    
        
    seg_mask = seg_mask_temp& probe_dat_mask;

    if mod(ii,2) == 0
        marker='x';
    else
        marker ='+';
    end

    if sum(seg_mask)>0
%             subplot(4,4,[1 2])
%             pl= plot(shot_idx_all(seg_mask_temp),setpts_all(seg_mask_temp),marker);
%             title('Setpoints by segment')
%             pl.Color=cdat(ceil(800*ii/iimax),:);
%             hold on

        subplot(4,4,[1 2])

        pl=plot(shot_time_rel(seg_mask)/(60*60),(to_seg_fits.set_sel{ii}-mean_setpt)*1e-9,marker);
        pl.Color=cdat(ceil(800*ii/iimax),:);
        hold on
        title('Setpts of good shots (GHz)')
        xlabel('Time (h)')
        ylabel(sprintf('\\omega-%.3f',mean_setpt*1e-9))

        subplot(4,4,[3 4])
        pl=plot(shot_time_rel(seg_mask)/(60*60),to_seg_fits.delta_sig{ii},marker);
        pl.Color=cdat(ceil(800*ii/iimax),:);
        hold on
        title('Signal of good shots')
        xlabel('Time (h)') 
    end
end
subplot(4,4,[3 4])

if sum(outlier_mask)>0
    hold on
    plot(shot_time_rel(outlier_mask)/(60*60),square_trap_freq_val(outlier_mask),'or');
    hold off  
end
        %remove the dead row
to_seg_fits.fit_all.to_freq.val(dead_rows,:)=[];
to_seg_fits.fit_all.to_freq.unc(dead_rows,:)=[];
to_seg_fits.fit_all.model(dead_rows,:)=[];

to_seg_fits.fit_trimmed.to_freq.val(dead_rows,:)=[];
to_seg_fits.fit_trimmed.to_freq.unc(dead_rows,:)=[];
to_seg_fits.fit_trimmed.model(dead_rows,:)=[];
to_seg_fits.fit_trimmed.to_unc_boot(dead_rows,:)=[];

%             to_seg_fits.set_sel = cell(iimax,1);
%             to_seg_fits.delta_sig = cell(iimax,1);
to_seg_fits.xdat(dead_rows,:)=[];
to_seg_fits.ydat(dead_rows,:)=[];
to_seg_fits.good_shot_idx(dead_rows,:)=[];
to_seg_fits.to_time(dead_rows,:)=[];
to_seg_fits.atom_num(dead_rows,:)=[];
to_seg_fits.avg_osc_coefs_probe(dead_rows,:)=[];

% Plot everything that was stored

to_time = to_seg_fits.to_time;
to_val_scale = 1e-9;
to_val = to_seg_fits.fit_all.to_freq.val*to_val_scale;
to_val_trim = to_seg_fits.fit_trimmed.to_freq.val*to_val_scale;
to_unc = to_seg_fits.fit_all.to_freq.unc*to_val_scale;
to_unc_trim = to_seg_fits.fit_trimmed.to_freq.unc*to_val_scale;
to_val_ref = to_val(1);
to_val_ref_trim = to_val_trim(1);


% get the slope out
slopes = cellfun(@(x) x.Coefficients.Estimate(2),to_seg_fits.fit_all.model)*scale_freq;
slopes_err = cellfun(@(x) x.Coefficients.SE(2),to_seg_fits.fit_all.model)*scale_freq;

to_seg_fits.fit_trimmed.slope.val=slopes;
to_seg_fits.fit_trimmed.slope.unc=slopes_err;


num_val = to_seg_fits.atom_num(:,1);
num_unc = to_seg_fits.atom_num(:,2);

to_diffs = abs(diff(to_val));
to_diffs_trim = abs(diff(to_val_trim));

[acf,lags,bound] = autocorr(to_val);
[acf_trim,lags_trim,bound_trim] = autocorr(to_val_trim);
plot_name='single_run_drift_anal';
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
%   autocorr(to_val_trim)
plot(lags*mean(diff(to_time)),acf_trim,'ro')

hold on
plot([min(lags*mean(diff(to_time))),max(lags*mean(diff(to_time)))],[0,0],'k');
title('Trim TO autocorrelation')
xlabel('Time (h)')
ylabel('Autocorrelation')

subplot(4,4,[11 12])
plot(shot_time_rel(1:last_good_shot)/(60*60),data.mcp_tdc.num_counts(1:last_good_shot))
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
ylabel('Gradient (Hz(trap)^2/Hz(opt))')

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
saveas(gcf,fullfile(anal_opts_fit_to_seg.global.out_dir,[plot_name,'.png']))




    
end

function out=wrap_fit_for_boot(xywdat,ii,scale)
    fit_out=fit_poly_with_int(xywdat,ii,0,1);
    out=fit_out.x_intercept.val/scale;
end

%       mdl_coef=mdl_all.Coefficients;
%        to_seg_fits.fit_all.to_freq.val(ii)=to_freq/scale_freq+opt_freq_offset;
%         to_seg_fits.fit_all.to_freq.unc(ii)=abs((mdl_coef.Estimate(1)/mdl_coef.Estimate(2)))*...
%             sqrt((mdl_coef.SE(1)/mdl_coef.Estimate(1))^2+...
%             (mdl_coef.SE(2)/mdl_coef.Estimate(2))^2)/scale_freq;


% %apply some slightly over the top ckecks
% delta_mask = ~isnan(delta_trap_freq_all);
% probe_dat_mask = data.osc_fit.ok.rmse & data.mcp_tdc.all_ok'...
% & ~isnan(data.osc_fit.trap_freq_recons) & delta_mask' & ~temp_cal;
% cal_dat_mask = data.osc_fit.ok.rmse & data.mcp_tdc.all_ok'...
% & ~isnan(data.osc_fit.trap_freq_recons) & delta_mask' & temp_cal;
% trap_freq=data.osc_fit.trap_freq_recons';
% cal_trap_freq=data.cal.freq_drift_model(data.mcp_tdc.time_create_write(:,1));
% square_trap_freq_val= (trap_freq).^2-(cal_trap_freq).^2;