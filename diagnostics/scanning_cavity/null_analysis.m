function null_data = null_analysis(config)
    % Set up for import

    dir_read=dir(fullfile(config.dir,config.log_name));
    file_names={dir_read.name};
    num_files = numel(file_names);

    all_pos_data = cell(num_files,1);
    for pp = 1:num_files
    % Load the data
        config.fname = file_names{pp};
        path=fullfile(config.dir,config.fname);
        fid = fopen(path,'r');
        raw_line = fread(fid,Inf,'*char')'; % this is the fastest string read method 3.5s for 20 files
        fclose(fid);
        ai_dat=jsondecode(raw_line);

        samples= size(ai_dat.Data,2);
        sr=ai_dat.sample_rate;
        acquire_time=samples/sr;
        T = linspace(0,acquire_time, samples);


        % Find the scan boundaries
        config.pzt_division=0.2498;
        sfp_pzt_raw=ai_dat.Data(3,:);
        sub_ptz_raw=sfp_pzt_raw(config.sampl_start:end)/config.pzt_division;
        kernel=gausswin(ceil(4*config.pzt_volt_smothing_time*sr),4);
        kernel=kernel/sum(kernel(:));%normalize
        sub_dat_smooth=conv(sub_ptz_raw,kernel,'same');
        sub_dat_grad_smooth=diff(sub_dat_smooth)*sr;
        pos_slope_mask=sub_dat_grad_smooth>0;
        pos_scan_data = ai_dat.Data(2,pos_slope_mask);

        % Store data
        all_pos_data{pp} = pos_scan_data';
    end
    % 


    null_zero = mean(pos_scan_data);
    null_std = std(pos_scan_data);
    fprintf(' - Average PD reading in %s: %.4f, std %.4f\n',config.dir,null_zero,null_std) 

    %Prepare output
    null_data.mean = null_zero;
    null_data.std = null_std;
    null_data.num_files = num_files;
    a = cell_vertcat(all_pos_data);
    null_data.data = a{1};
    %% Graphical output

    
    if config.plot_out
        
            sfigure(1);
            all_scan_data = cell2mat(cell_vertcat(all_pos_data));
            histogram(all_scan_data,100);

            title('Null scan voltage data')
            xlabel('Voltage')
            ylabel('Counts')
        
%         % % Inspect time series    
%         calib_(1) = sfigure();
%         clf;
%         subplot(2,2,1)
%         plot(pd_raw)
%         title('raw time series')
%         xlabel('time')
%         ylabel('Voltage')
% 
%         subplot(2,2,2)
%         histogram(pd_raw,150)
%         title('Histogram of values')
%         xlabel('Voltage')
%         ylabel('counts ')
% 
%         subplot(2,2,3)
%         histogram(pos_scan_data,75)
%         title('Histogram of slow scan')
%         xlabel('Voltage')
%         ylabel('counts ')
% 
%         subplot(2,2,4)
%         histogram(neg_scan_data,75)
%         title('Histogram of fast scan')
%         xlabel('Voltage')
%         ylabel('counts ')
    end

end