function single_out=is_single_scan_sm(single_opts)

%TODO
% -[x] remove peaks very close to the edge

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


%need to find peaks in time,pd space bc in pzt,pd space bc the pzt is not guaranteed to be strictly increasing
%to prevent having to do [v, idx] = closest_value(single_opts.time, pks_time) we just use the single input peak finder
% this only looses the option of using the peak widths in the future
[pks_pd,pks_idx] = findpeaks(pd_dat_uncmp_filt,'MinPeakHeight', single_opts.peak_thresh(1));
pks_pzt=ptz_voltage_filt(pks_idx);
%mask out the peaks that are too close to the upper or lower limits of the pzt scan
mask= pks_pzt<(max(ptz_voltage_filt)-single_opts.scan_clip_pzt)  & pks_pzt>(min(ptz_voltage_filt)+single_opts.scan_clip_pzt);
single_out.pks.full.pd=pks_pd(mask);
single_out.pks.full.pzt=pks_pzt(mask);
single_out.pd_full_range=range(pd_dat_uncmp_filt);

[pks_pd,pks_idx] = findpeaks(pd_dat_cmp_filt,'MinPeakHeight',single_opts.peak_thresh(1));
pks_pzt=ptz_voltage_filt(pks_idx);
mask= pks_pzt<(max(ptz_voltage_filt)-single_opts.scan_clip_pzt)  & pks_pzt>(min(ptz_voltage_filt)+single_opts.scan_clip_pzt);
single_out.pks.cmp.pd=pks_pd(mask);
single_out.pks.cmp.pzt=pks_pzt(mask);

%find the distance in pzt voltage to the nearest full peak
min_dist_cmp_to_full=arrayfun(@(x) min(abs(x-single_out.pks.full.pzt)),single_out.pks.cmp.pzt);
cmp_pks_mask=min_dist_cmp_to_full>single_opts.merge_cmp_full_pk_dist;
single_out.pks.cmp.pd=single_out.pks.cmp.pd(cmp_pks_mask);
single_out.pks.cmp.pzt=single_out.pks.cmp.pzt(cmp_pks_mask);

single_out.pks.all.pzt=[single_out.pks.cmp.pzt,single_out.pks.full.pzt];
single_out.pks.all.pd=[single_out.pks.cmp.pd,single_out.pks.full.pd];
single_out.pks.all.min_pzt_v_diff=min_diff_vec(single_out.pks.all.pzt);

if isnan(single_out.pks.all.min_pzt_v_diff)
    single_out.is_single_mode=false;
else
single_out.is_single_mode= single_out.pks.all.min_pzt_v_diff> single_opts.pzt_dist_sm;
end

if ~isempty(single_out.pd_full_range)
    single_out.is_power_ok=single_out.pd_full_range>single_opts.pd_amp_min;
else
    single_out.is_power_ok=false;
end

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
    plot(ptz_voltage_filt,pd_dat_uncmp,'k')
    xlabel('pzt (filt)(v)')
    ylabel('pd (v)')
    hold on
    plot(ptz_voltage_filt,pd_dat_uncmp_filt,'b')
    plot(single_out.pks.full.pzt,single_out.pks.full.pd,'xk','markersize',20)
    if isnan(single_opts.cmp_multiplier_disp)
        single_opts.cmp_multiplier_disp=max(pd_dat_uncmp)/max(pd_dat_cmp);
    end
    plot(ptz_voltage_filt,pd_dat_cmp*single_opts.cmp_multiplier_disp,'r')
    plot(single_out.pks.cmp.pzt,single_out.pks.cmp.pd*single_opts.cmp_multiplier_disp,'rx','markersize',20);
    hold off
    title(title_strs(single_out.is_single_mode+1))
    drawnow
    refresh(1)
end %end plots

end