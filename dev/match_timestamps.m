<<<<<<< HEAD
function data_mcp_tdc = match_timestamps(data,anal_opts)
    if anal_opts.match_times
=======
function data_mcp_tdc = match_timestamps(data,anal_opts,match_times)
    if match_times
>>>>>>> master
        %try and match up the file with if it is a calibaration using the time
        %it is slightly overkill here to search each one, but just being extra
        %cautious/flexible
        time_thresh=4; %how close for the times to be considered the same shot
        %lets examine what the time difference does
        sfigure(45);
        set(gcf,'color','w')
        clf
        imax=min([size(data.labview.time,2),size(data.mcp_tdc.time_create_write,1)]);
        %imax=5000;
        time_diff=data.mcp_tdc.time_create_write(1:imax,2)'-anal_opts.dld_aquire-anal_opts.trig_dld-...
            data.labview.time(1:imax);
        mean_delay_labview_tdc=mean(time_diff);
        plot(time_diff)
        xlabel('shot number')
        ylabel('time between labview and mcp tdc')
        %to do include ai_log
        iimax=size(data.mcp_tdc.time_create_write(:,1),1);
        data_mcp_tdc.probe.calibration=nan(iimax,1);
        data_mcp_tdc.labview_shot_num=nan(iimax,1);
        %loop over all the tdc_files 
        for ii=1:iimax
            %predict the labview master trig time
            %use the write time to handle being unpacked from 7z
            est_labview_start=data.mcp_tdc.time_create_write(ii,2)...
                -anal_opts.trig_dld-anal_opts.dld_aquire-mean_delay_labview_tdc;
            [tval,nearest_idx]=closest_value(data.labview.time...
                ,est_labview_start);
            if abs(tval-est_labview_start)<time_thresh
                data_mcp_tdc.labview_shot_num(ii)=data.labview.shot_num(nearest_idx);
                data_mcp_tdc.probe.calibration(ii)=data.labview.calibration(nearest_idx);
            end 
        end
    else
        %just do it the boring way and hope that the tdc was set up right and
        %there were no false trigers
        %TO DO

    end
end