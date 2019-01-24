function scan_data = process_scan(scan,config)
        scan = scan - config.zero_offset;
        T_cen = config.T_sweep(config.p_cent);
        line_lim = config.peak_width;
        line_lims = [T_cen - line_lim,T_cen + line_lim];
        T_sweep(idx_lims(1):idx_lims(2))
        peak_mask = T_win > line_lims(1) & T_win< line_lims(2);
        dT = mean(diff(T_win)); %Not correct but need a quick fix
        peak_power = sum(scan(peak_mask)*dT);
        back_power = sum(scan(~peak_mask)*dT);
        scan_data.T_fit = linspace(min(T_win),max(T_win),2000);
        scan_data.ratios(jj) = peak_power/back_power;
        scan_data.vals(jj,:) = [peak_power,back_power];
end