function light_on_data = light_on_analysis(config)
    dir_read=dir(fullfile(config.dir,config.log_name));
    file_names={dir_read.name};
    config.zero_offset = config.zero_offset;
    if ~isnan(config.test.num_files)
        num_files = config.test.num_files;
    else
         num_files = numel(file_names);
    end

    all_vals = cell(num_files,1);
    raw_scans_all = cell(num_files,1);
    dir_ratios = cell(num_files,1);
    
    for pp = 1:num_files
        
        % Import data
        config.fname = file_names{pp};
        path=fullfile(config.dir,config.fname);
        fid = fopen(path,'r');
        raw_line = fread(fid,Inf,'*char')'; % this is the fastest string read method 3.5s for 20 files
        fclose(fid);
        ai_dat=jsondecode(raw_line);
    
        % Process
        samples= size(ai_dat.Data,2);
        sr=ai_dat.sample_rate;
        acquire_time=samples/sr;
        T = linspace(0,acquire_time, samples);
        pd_raw = ai_dat.Data(2,:);
        sfp_pzt_raw=ai_dat.Data(3,:);
        sub_ptz_raw=sfp_pzt_raw(1:end)/config.pzt_division;
        kernel=gausswin(ceil(4*config.pzt_volt_smothing_time*sr),4);
        kernel=kernel/sum(kernel(:));%normalize
        sub_dat_smooth=conv(sub_ptz_raw,kernel,'same');
        sub_dat_grad_smooth=diff(sub_dat_smooth)*sr;
        pos_slope_mask=sub_dat_grad_smooth>0;
        
        %now we want to select the first positive going sgement that is long enough
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
            
        % Init
        file_vals = cell(num_scans,1);
        raw_scan_peaks = cell(num_scans,1);
        file_ratios = cell(num_scans,1);

        if config.plot_out
            figure();
        end
 
        for ii=1:num_scans
            idx_lims = borders(ii,:);
            sweep = pd_raw(idx_lims(1):idx_lims(2));
            
            [pks_val,locs] = findpeaks(sweep,'MinPeakHeight',config.treshold);
            if ~isempty(locs)
                hwin_size = floor(0.5*median(diff(locs)));
                num_peaks = length(locs);
                scan_ratios = zeros(num_peaks, 1);
                scan_vals = zeros(num_peaks,2);
                raw_scan = cell(num_peaks,1);
                scan_config = config;
                scan_config.T_sweep = T(idx_lims(1):idx_lims(2));
                for jj = 1:num_peaks
%                     scan_data{jj} = process_scan(scan);
                    scan_config.p_cent = locs(jj);                 
                    scan_lims = [(max(1,p_cent-hwin_size)), min(length(sweep),p_cent + hwin_size)];
                    scan_config.T_win = T_sweep(scan_lims(1):scan_lims(2));
%                     p_edge = T_sweep(idx_lims);
                    if diff(scan_lims) > 1.5*hwin_size % ensure all or most of a scan is taken
                        scan = sweep(scan_lims(1):scan_lims(2));
                        scan_data{jj} = process_scan(scan,scan_config);
                    end 
                end %loop over peaks  
            end %isempty(locs)
            file_vals{ii} = scan_vals;
            raw_scan_peaks{ii} = raw_scan;
            file_ratios{ii} = scan_ratios;
        end %loop over scans 
        raw_scans_all{pp} = raw_scan_peaks;
%         dir_stats(pp,:) = [mean(scan_ratios),std(scan_ratios)];
        dir_ratios{pp} = file_ratios;
        all_vals{pp} = cell_vertcat(file_vals);        
    end
    a = cell_vertcat(dir_ratios);
    a = cell_vertcat(a{1});
    light_on_data.dir_ratios = a{1};
    a = cell_vertcat(all_vals);
    a = cell_vertcat(a{1});
    light_on_data.all_vals = a{1};
    b = cell_vertcat(raw_scans_all);
    b = cell(vertcat(b{1}));
    light_on_data.raw_data = b{1};
    
    %% ASSUMING PEAKS FIXED CORRECTLY
    light_on_data.stats = [mean(light_on_data.dir_ratios),std(light_on_data.dir_ratios)];
    
end

