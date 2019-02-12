function light_on_data = light_on_analysis(config)
    dir_read=dir(fullfile(config.dir,config.log_name));
    file_names={dir_read.name};
    config.zero_offset = config.zero_offset;
    if ~isnan(config.test.num_files)
        num_files = config.test.num_files;
    else
        num_files = numel(file_names);
    end

    dir_data = cell(num_files,1);
    for pp = 1:num_files
        
        % Import data
        config.fname = file_names{pp};
        path=fullfile(config.dir,config.fname);
        fid = fopen(path,'r');
        raw_line = fread(fid,Inf,'*char')'; % this is the fastest string read method 3.5s for 20 files
        fclose(fid);
        ai_dat=jsondecode(raw_line);
    
        % get set up
        samples= size(ai_dat.Data,2);
        sr=ai_dat.sample_rate;
        acquire_time=samples/sr;
        T = linspace(0,acquire_time, samples);
        pd_raw = ai_dat.Data(2,:);
%         pd_raw = pd_raw;% - config.zero_offset;
        sfp_pzt_raw=ai_dat.Data(3,:);
        sub_ptz_raw=sfp_pzt_raw(1:end)/config.pzt_division;
        
        % Find the forward scans
        % Smoothing
        kernel=gausswin(ceil(4*config.pzt_volt_smothing_time*sr),4);
        kernel=kernel/sum(kernel(:));%normalize
        sub_dat_smooth=conv(sub_ptz_raw,kernel,'same');
        sub_dat_grad_smooth=diff(sub_dat_smooth)*sr;
        % Find positive slopes
        pos_slope_mask=sub_dat_grad_smooth>0;
        x=pos_slope_mask';
        % Find the transitions
        downs = strfind(x',[1 0]);
        ups = strfind(x',[0 1]);
        %Use only complete scans
        T_min = min(ups);
        T_max = max(downs);
        downs = downs(downs > T_min);
        ups = ups(ups < T_max);
        borders = [ups',downs'];
        if ~isnan(config.test.num_scans)
            num_scans = config.test.num_scans;
        else
            num_scans = size(borders,1);
        end        
        
        % Set up for scan iterator
        if config.pv_plot
            sfigure(68);
            clf;
            sfigure(69);
            clf;
        end
        scan_config = config;
        for ii=1:num_scans
            scan_config.idx_lims = borders(ii,:); % marks the edges of each scan
            sc_data.pd_raw = pd_raw;
            sc_data.T = T;
            sc_data.sub_dat_smooth = sub_dat_smooth;
            scan_data{ii} = process_scan(sc_data,scan_config);
        end %loop over scans 

        light_on_data{pp} = scan_data; 
%         light_on_data{pp}.raw = pd_raw(pos_slope_mask) - config.zero_offset;
    end


    
end

function scan_data = process_scan(sc_data,config)
    % Takes a single forward scan of the SFP and returns a list of the
    % peak/BG ratio (will change later)
    % Input: scan (array of photodiode readings)
    %        config.

    % INIT OUTPUT 
%     scan_data = cell(num_peaks,1);
    scan_data = [];
    scan_config = config;
    pd_raw = sc_data.pd_raw;
    T = sc_data.T;
    sub_dat_smooth = sc_data.sub_dat_smooth;
    idx_lims = config.idx_lims;
    % Extract the scan of interest & find peaks
    scan_data.sweep = pd_raw(idx_lims(1):idx_lims(2)) - config.zero_offset;
    scan_data.T_sweep = T(idx_lims(1):idx_lims(2)); % Timestamps in the scan
    scan_data.ptz_raw = sub_dat_smooth(idx_lims(1):idx_lims(2));
    [~,scan_config.locs] = findpeaks(scan_data.sweep,'MinPeakHeight',config.treshold);
    
%     if ~isempty(locs) % If there are peaks found
        
        % Whole-scan methods
        % Lebesgue
        if config.lebesgue_method
            config_lebesgue = scan_config;
            scan_data.lebesgue = lebesgue_method(scan_data,config_lebesgue);
        end
        if config.pv_method
            config_pv = scan_config;
            scan_data.pv_data = pv_method(scan_data,config_pv);
        end
        
end



function lebesgue_data = lebesgue_method(data,config)

    T = data.T_sweep;
    dT = mean(diff(T));
    X = data.sweep;
    X_s= sort(X);

    varfrac_ratios = zeros(100,1);
    varfrac_axis = zeros(100,1);
    varfrac_back = zeros(100,1);
    
    for ff=1:100
        frac = ff/100;
        Xs_mask = X_s<frac*config.lebesgue_thresh_demo;
        X_back = X_s(Xs_mask);
        X_peak = X_s(~Xs_mask);
        P_back = sum(X_back*dT);
        P_peak = sum(X_peak*dT);
        varfrac_back(ff) = P_back;
        varfrac_ratios(ff) = P_back/P_peak;
        varfrac_axis(ff) = frac*config.lebesgue_thresh_demo;
    end

    Xs_mask = X_s<config.lebesgue_thresh;
    X_back = X_s(Xs_mask);
    X_peak = X_s(~Xs_mask);
    P_back = sum(X_back*dT);
    P_peak = sum(X_peak*dT);
    lebesgue_data.P_back = P_back;
    lebesgue_data.P_ratio = P_back/P_peak;
    lebesgue_data.varf_ratios = varfrac_ratios;
    lebesgue_data.varf_axis = varfrac_axis;
    lebesgue_data.varf_back = varfrac_back;
    lebesgue_data.Xs = X_s;
    
end

function pv_data = pv_method(data,config)

    X = data.sweep;
    T = data.T_sweep;
    dT = mean(diff(T));
    num_peaks = length(config.locs);
    num_pairs = num_peaks - 1;
    width = config.peak_width;

    
    % Determine FWHM

    peak_widths = zeros(num_peaks,1);
    for ii=1:num_peaks
        idx_min= max(1,config.locs(ii)-floor(config.window_size/dT));
        idx_max= min(length(X),floor(config.window_size/dT)+config.locs(ii));
        T_cent = T(idx_min:idx_max)-T(config.locs(1));
        T_mid = median(T_cent);        
        coefs = [1,width,T_mid,0];
%         lorentz_fit = fitnlm(T_cent,X(idx_min:idx_max),@Lorentzian,coefs);
%         coefs = lorentz_fit.Coefficients.Estimate;
        
        if config.pv_plot
            sfigure(68);
            subplot(2,1,1)
            plot(T_cent,X(idx_min:idx_max),'k');
            hold on
            plot(T_cent,Lorentzian(coefs,T_cent),'r-.')
            subplot(2,1,2)
            plot(T_cent,X(idx_min:idx_max),'k');
            hold on
            plot(T_cent,Lorentzian(coefs,T_cent),'r-.')
            set(gca,'Yscale','log')
        end
    
        peak_widths(ii) = coefs(2);
    end


    peak_pows = zeros(num_pairs,1);
    back_pows = zeros(num_pairs,1);
    ratios = zeros(num_pairs,1);
    for ii=1:num_pairs
        width = config.pv_peak_factor*peak_widths(ii);
        idx_min= max(1,-floor(width/dT)+config.locs(ii));
        idx_max= min(floor(width/dT)+config.locs(ii+1),length(T));
        T_pp = T(idx_min:idx_max);
        T_mid = median(T_pp); 
        T_offset = T_pp-T_mid;
        


        %Integrate the peak & valley
        
        p1_idxs = [max(1,-floor(width/dT)+config.locs(ii)),min(floor(width/dT)+config.locs(ii),length(T))];
        p2_idxs = [max(1,-floor(width/dT)+config.locs(ii+1)),min(floor(width/dT)+config.locs(ii+1),length(T))];
        p1_vals = X(p1_idxs(1):p1_idxs(2));
        p2_vals = X(p2_idxs(1):p2_idxs(2));
        p1_pow = sum(p1_vals*dT);
        p2_pow = sum(p2_vals*dT);
        mn_pow = (p1_pow+p2_pow)/2;
        
        if config.pv_plot
            sfigure(69);
            subplot(2,2,1)
            plot(T_offset,X(idx_min:idx_max))
            hold on
            subplot(2,2,2)
            plot(T_offset,X(idx_min:idx_max))
            hold on
            set(gca,'Yscale','log')
            subplot(2,2,3)
            plot(T(p1_idxs(1):p1_idxs(2)),p1_vals)
            hold on
            subplot(2,2,4)
            plot(T(p1_idxs(1):p1_idxs(2)),p2_vals)
            hold on
        end
        
        % Integrate background power over same width
        
        bg_cent = floor(median(config.locs(ii:ii+1)));
        bg_idxs = bg_cent + floor(width/dT)*[-1,1];
        bg_idxs = [max(1,bg_idxs(1)),min(bg_idxs(2),length(X))];
        bg_vals = X(bg_idxs(1):bg_idxs(2));
        bg_pow  = sum(bg_vals*dT);
        
        peak_pows(ii) = mn_pow;
        back_pows(ii) = bg_pow;
        ratios(ii) = bg_pow/mn_pow;
        
            
    end
    pv_data.peaks = peak_pows;
    pv_data.backs = back_pows;
    pv_data.ratio = ratios;
end


function peak_data = peak_process(scan,config)
        scan = scan - config.zero_offset;
        T_cen = config.T_cen;
        T_win = config.T_win;
        line_lim = config.peak_width;
        dT = config.dT;
        
        line_lims = [T_cen - line_lim,T_cen + line_lim];
        peak_mask = T_win > line_lims(1) & T_win< line_lims(2);
       
        peak_power = sum(scan(peak_mask)*dT);
        back_power = sum(scan(~peak_mask)*dT);
        
        peak_data.powers = [peak_power,back_power];
        peak_data.raw = scan;
        peak_data.ratio_b2p = peak_power/back_power;
        peak_data.T_win = T_win;
end