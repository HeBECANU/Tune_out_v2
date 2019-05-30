function corr_cancel_out=corr_ac_mains(data,anal_opts)
% Look for correlations of trap osc with mains
% try to remove the correlations

%function corr_ac_mains(data)
%end

if anal_opts.do
    
    terms=1; %number of terms to use in the mains waveform
    time_shifts=col_vec(linspace(anal_opts.tlim(1),anal_opts.tlim(2),1+round(range(anal_opts.tlim)/anal_opts.tsamp)));
    %construct a more convinent temp variable txyz_tmp wich is the position in mm for use in the fit
    x_tmp=data.mcp_tdc.al_pulses.vel_zxy.mean(:,:,3);
    x_tmp(~data.mcp_tdc.all_ok,:)=nan;
    t_tmp=col_vec(data.mcp_tdc.al_pulses.time_cen);    
    cancel_param=[0,0]; %[0.371686,0.00009]; %time,amp
    
    main_osc_corr=calc_corr(t_tmp,x_tmp,time_shifts,cancel_param, data.ai_log.ac_mains_fit,terms,4);
    fprintf('max norm correlation before cancelation %f\n',main_osc_corr.max_corr_norm_val)
    cancel_param=[main_osc_corr.max_corr_norm_time,main_osc_corr.max_corr_norm_val];
    stfig('mains_corr','add_stack',1);
    subplot(2,1,1)
    plot(time_shifts,main_osc_corr.xcorr_norm)
    xlabel('time shift')
    ylabel('corr amp')
    title('before cancelation')
    drawnow

    main_osc_cancel_corr=calc_corr(t_tmp,x_tmp,time_shifts,cancel_param, data.ai_log.ac_mains_fit,terms,4);
    fprintf('max norm correlation after cancelation  %f\n',main_osc_cancel_corr.max_corr_norm_val)
    fprintf('improvement in corr amp %f \n',main_osc_corr.max_corr_norm_val/main_osc_cancel_corr.max_corr_norm_val)

    stfig('mains_corr','add_stack',1);
    subplot(2,1,2)
    plot(time_shifts,main_osc_cancel_corr.xcorr_norm)
    xlabel('time shift')
    ylabel('corr amp')
    title('after cancelation')
    drawnow

    corr_cancel_out=[];
    corr_cancel_out.vel_zxy_corr_cancel=nan(size(x_tmp,1),size(x_tmp,2),3);
    corr_cancel_out.vel_zxy_corr_cancel(:,:,2)=main_osc_cancel_corr.x_canceled_array;
     
end

end

%    x_tmp=x_tmp-nanmean(x_tmp);
%     y_tmp=col_vec(squeeze(data.mcp_tdc.al_pulses.vel.mean(ii,:,3)));
%     y_tmp=y_tmp-nanmean(y_tmp);
%     z_tmp=col_vec(squeeze(data.mcp_tdc.al_pulses.vel.mean(ii,:,1)));
%     z_tmp=z_tmp-nanmean(z_tmp);

function out=calc_corr(tdat,xdat,time_shifts,cancel_param,models,terms,verbose)
    jjmax=numel(time_shifts);
    iimax=size(xdat,1);
    x_corr_osc_mains_norm=nan(iimax,jjmax);
    x_corr_osc_mains_unorm=nan(iimax,jjmax);
    x_canceled_array=xdat*nan;
    for ii=1:iimax
        ac_fit=models(ii);
        x_vec=col_vec(xdat(ii,:));
        if sum(isnan(x_vec))==0 %skip if nan in x_vec
            if isequal(cancel_param,[0,0])
                x_canceled_tmp=x_vec;
            else
                ac_tmp=ac_fit.models.var_terms(tdat+cancel_param(1),1)*cancel_param(2);
                ac_tmp=ac_tmp-mean(ac_tmp);
                x_canceled_tmp=x_vec-ac_tmp;
            end
            x_canceled_array(ii,:)=x_canceled_tmp';
            x_canceled_mean=x_canceled_tmp-mean(x_vec); %remove mean for xcorr
            for jj=1:jjmax
                ac_values=ac_fit.models.var_terms(tdat+time_shifts(jj),terms);
                ac_values=ac_values-mean(ac_values); %remove mean for xcorr
                x_corr_osc_mains_unorm(ii,jj)=sum(ac_values.*x_canceled_mean);
                x_corr_osc_mains_norm(ii,jj)=x_corr_osc_mains_unorm(ii,jj)/sum(ac_values.*ac_values);

            end
            if mod(ii,10)==0 && verbose>4
                stfig('mains_corr','add_stack',1);
                clf
                plot(time_shifts,nanmean(x_corr_osc_mains_norm,1))
                drawnow
                pause(1e-3)
            end
        end
    end

    out.xcorr_unorm=nanmean(x_corr_osc_mains_unorm,1);
    % find the largest correlation 
    [~,max_corr_idx]=max(abs(out.xcorr_unorm));
    out.max_corr_unorm_val=out.xcorr_unorm(max_corr_idx);
    out.max_corr_unorm_time=time_shifts(max_corr_idx);
    
    out.xcorr_norm=nanmean(x_corr_osc_mains_norm,1);
    [~,max_corr_idx]=max(abs(out.xcorr_norm));

    out.max_corr_norm_val=out.xcorr_norm(max_corr_idx);
    out.max_corr_norm_time=time_shifts(max_corr_idx);
    out.x_canceled_array=x_canceled_array;
    if verbose>2
        fprintf('max unorm correlation at time %f with amp %f\n',out.max_corr_unorm_time,out.max_corr_unorm_val)
        fprintf('max norm correlation at time %f with amp %f\n',out.max_corr_norm_time,out.max_corr_norm_val)
    end
    if verbose>3
        stfig('mains_corr','add_stack',1);
        clf
        subplot(2,1,1)
        plot(time_shifts,out.xcorr_norm)
        subplot(2,1,2)
        plot(time_shifts,out.xcorr_unorm)
    end

end