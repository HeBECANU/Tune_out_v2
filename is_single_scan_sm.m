function single_out=is_single_scan_sm(single_opts)




%%
%load('is_single_sm_start.mat');

pd_dat_uncmp=single_opts.pd_voltage(:,1);
pd_dat_cmp=single_opts.pd_voltage(:,2);
pzt_voltage=single_opts.pzt_voltage;

pd_dat_uncmp_filt=gaussfilt(single_opts.time,pd_dat_uncmp,...
    single_opts.pd_filt_factor*single_opts.pzt_scan_period);
pd_dat_cmp_filt=gaussfilt(single_opts.time,pd_dat_cmp,...
    single_opts.pd_filt_factor*single_opts.pzt_scan_period);
ptz_voltage_filt=gaussfilt(single_opts.time,pzt_voltage,...
    single_opts.pzt_filt_factor*single_opts.pzt_scan_period);


[pks_pd,pks_pzt] = findpeaks(pd_dat_uncmp_filt,ptz_voltage_filt,'MinPeakHeight', single_opts.peak_thresh(1));
single_out.pks.full.pd=pks_pd;
single_out.pks.full.pzt=pks_pzt;
single_out.pd_full_range=range(pd_dat_uncmp_filt);

[pks_pd,pks_pzt] = findpeaks(pd_dat_cmp_filt,ptz_voltage_filt,'MinPeakHeight',single_opts.peak_thresh(1));
single_out.pks.cmp.pd=pks_pd;
single_out.pks.cmp.pzt=pks_pzt;

%find the distance in pzt voltage to the nearest full peak
min_dist_cmp_to_full=arrayfun(@(x) min(abs(x-single_out.pks.full.pzt)),single_out.pks.cmp.pzt);
cmp_pks_mask=min_dist_cmp_to_full>single_opts.merge_cmp_full_pk_dist;
single_out.pks.cmp.pd=single_out.pks.cmp.pd(cmp_pks_mask);
single_out.pks.cmp.pzt=single_out.pks.cmp.pzt(cmp_pks_mask);

single_out.pks.all.pzt=[single_out.pks.cmp.pzt,single_out.pks.full.pzt];
single_out.pks.all.pd=[single_out.pks.cmp.pd,single_out.pks.full.pd];
single_out.pks.all.min_pzt_v_diff=min_diff(single_out.pks.all.pzt);
        
single_out.is_single_mode= single_out.pks.all.min_pzt_v_diff> single_opts.peak_distv_min_pass;
single_out.is_power_ok=single_out.pd_full_range>single_opts.pd_amp_min;

if single_opts.plot.all || (single_opts.plot.failed && ~single_out.is_single_mode)
    
    
    title_strs={'Failed','OK'};
    sfigure(1);
    clf;
    subplot(2,2,1)
    cla;
     set(gcf,'color','w')
    title('smoothing pzt v')
    plot(single_opts.time,pzt_voltage,'k')
    hold on
    plot(single_opts.time,ptz_voltage_filt,'b')
    hold off
    xlabel('time (ms)')
    ylabel('volts (v)')


    subplot(2,2,3)
    cla;
    plot(pzt_voltage,pd_dat_uncmp,'k')
    xlabel('pzt(v)')
    ylabel('pd (v)')
    hold on
    plot(ptz_voltage_filt,pd_dat_uncmp_filt,'b')
    plot(single_out.pks.full.pzt,single_out.pks.full.pd,'xk','markersize',20)
    plot(ptz_voltage_filt,pd_dat_cmp*single_opts.cmp_multiplier_disp,'r')
    plot(single_out.pks.cmp.pzt,single_out.pks.cmp.pd*single_opts.cmp_multiplier_disp,'rx','markersize',20);
    hold off
    title(title_strs(single_out.is_single_mode+1))

    pause(1e-3)
end %end plots

end