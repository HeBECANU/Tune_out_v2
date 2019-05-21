% Look for correlations of trap osc with mains
% try to remove the correlations

%function corr_ac_mains(data)
%end

iimax=size(data.mcp_tdc.counts_txy,2); 

terms=1;

tlim=[-1,1];
time_shifts=col_vec(linspace(tlim(1),tlim(2),range(tlim)/1e-3));
jjmax=numel(time_shifts);
x_corr_osc_mains_norm=nan(iimax,jjmax);
x_corr_osc_mains_unorm=nan(iimax,jjmax);
norm_corr=@(u,v) sum(u.*v)/sqrt(sum(u.^2)*sum(v.^2));
corr=@(u,v) sum(u.*v);

cancel_param=[0.371686,0.00009]; %[0.371686,0.00009];


for ii=1:iimax
    %position that data appears in data.mcp_tdc, not ness shot number
    %specified because we will remove any elements of osc_fit that did not
    %fit because of all_ok condition
    osc_fit.dld_shot_idx(ii)=ii;
    %shot number eg d123.txt as recorded by the tdc computer, not ness lv
    %number
    dld_shot_num=data.mcp_tdc.shot_num(ii);
    if data.mcp_tdc.all_ok(ii)
         osc_fit.dld_shot_num(ii)=dld_shot_num;
        %construct a more convinent temp variable txyz_tmp wich is the position in mm for use in the fit
        x_tmp=col_vec(squeeze(data.mcp_tdc.al_pulses.vel.mean(ii,:,2)));
        x_tmp=x_tmp-nanmean(x_tmp);
        y_tmp=col_vec(squeeze(data.mcp_tdc.al_pulses.vel.mean(ii,:,3)));
        y_tmp=y_tmp-nanmean(y_tmp);
        z_tmp=col_vec(squeeze(data.mcp_tdc.al_pulses.vel.mean(ii,:,1)));
        z_tmp=z_tmp-nanmean(z_tmp);
        txyz_tmp=cat(2,col_vec(data.mcp_tdc.al_pulses.time_cen),x_tmp,y_tmp,z_tmp);
        ac_data=data.ai_log.ac_mains_fit(ii);
        x_canceled=x_tmp-ac_data.models.var_terms(txyz_tmp(:,1)+cancel_param(1),1)*cancel_param(2);

        for jj=1:jjmax
            ac_values=ac_data.models.var_terms(txyz_tmp(:,1)+time_shifts(jj),terms);
            x_corr_osc_mains_norm(ii,jj)=norm_corr(ac_values,x_canceled);
            x_corr_osc_mains_unorm(ii,jj)=corr(ac_values,x_canceled);
        end
    end
    if mod(ii,10)==0
        stfig('mains_corr','add_stack',1);
        clf
        plot(time_shifts,nanmean(x_corr_osc_mains_norm,1))
        drawnow
        pause(1e-3)
    end
end

mean_xcorr_unorm=nanmean(x_corr_osc_mains_unorm,1);
[max_corr_val,max_corr_idx]=max(mean_xcorr_unorm);
fprintf('max unorm correlation at time %f with amp %f\n',time_shifts(max_corr_idx),max_corr_val)

mean_xcorr_norm=nanmean(x_corr_osc_mains_norm,1);
[max_corr_val,max_corr_idx]=max(mean_xcorr_norm);
fprintf('max norm correlation at time %f with amp %f\n',time_shifts(max_corr_idx),max_corr_val)
stfig('mains_corr','add_stack',1);
clf
subplot(2,1,1)
plot(time_shifts,mean_xcorr_norm)
subplot(2,1,2)
plot(time_shifts,mean_xcorr_unorm)