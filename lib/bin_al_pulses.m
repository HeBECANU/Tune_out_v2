function al_pulses=bin_al_pulses(anal_opts,data)


iimax=size(data.mcp_tdc.counts_txy,2);
al_pulses=[];
al_pulses.pulsedt=anal_opts.pulsedt;
al_pulses.window=nan(anal_opts.pulses,3,2); %initalize
al_pulses.num_counts=nan(iimax,anal_opts.pulses);
al_pulses.pos.mean=nan(iimax,anal_opts.pulses,3);
al_pulses.pos.std=nan(iimax,anal_opts.pulses,3);
global const %get gravity
fprintf('binning pulses in files %04u:%04u',size(data.mcp_tdc.counts_txy,2),0)
first_good_shot=true;
for shot=1:iimax
        if data.mcp_tdc.all_ok(shot)
            for pulse=1:anal_opts.pulses
                %set up time window centered arround t0
                t_pulse_cen=anal_opts.t0+anal_opts.pulsedt...
                    *(anal_opts.start_pulse+pulse-2);
                trange=t_pulse_cen+anal_opts.pulse_twindow*[-0.5,0.5];
                pulse_win_txy=[trange;anal_opts.xylim]; 
                counts_pulse=masktxy(data.mcp_tdc.counts_txy{shot},pulse_win_txy);
                if anal_opts.plot.all
                    stfig;
                    set(gcf,'Color',[1 1 1]);
                    subplot(3,1,1)
                    histogram(counts_pulse(:,1),100)
                    xlabel('t')
                    title('full')
                    subplot(3,1,2)
                    histogram(counts_pulse(:,2),100)
                    xlabel('x')
                    title('full')
                    subplot(3,1,3)
                    histogram(counts_pulse(:,3),100)
                    xlabel('y')
                    title('full')
                    pause(0.01)
                end
                if first_good_shot
                    %only need to store this on first shot becasue the same for
                    %all shots
                    al_pulses.window(pulse,:,:)=pulse_win_txy; 
                    al_pulses.time_cen(pulse,:)=t_pulse_cen;
                end
                al_pulses.num_counts(shot,pulse)=size(counts_pulse(:,3),1);
                al_pulses.pos.mean(shot,pulse,:)=mean(counts_pulse,1);
                al_pulses.pos.std(shot,pulse,:)=std(counts_pulse,1);
                %TODO: convert to velocity here and then take the mean,std
                vzxy_out=txy_to_vel(counts_pulse,...
                                    t_pulse_cen-anal_opts.global.fall_time,...
                                    -const.g0,...
                                    anal_opts.global.fall_dist);
                al_pulses.vel.mean(shot,pulse,:)=mean(vzxy_out,1);
                al_pulses.vel.std(shot,pulse,:)=std(vzxy_out,1);               
                
                
            end%pulse
            if first_good_shot,first_good_shot=false; end
        %check that the calculated t0 is close enough
        
        tol_t0_match=1e-2; %in factors of bin size
        tol_t0_match=tol_t0_match*anal_opts.pulsedt;
        mean_cen_time=mean(al_pulses.pos.mean(shot,:,1)-al_pulses.time_cen');
        if abs(mean_cen_time)>tol_t0_match
            est_t0=anal_opts.t0+mean_cen_time;
            warning('pulses are not centered in time pehaps t0 should be %.5f',est_t0)
        end
        end%is data.mcp_tdc.all_ok
        if mod(shot,10)==0,fprintf('\b\b\b\b%04u',shot),end     
%to set the pulse t0 right it can be handy to uncomment the next line
%
end%shots

%check if the average difference between the mean count position (vs the expected center) over all shots is outside a tolerance
tol_t0_match=1e-3; %in factors of bin size
tol_t0_match=tol_t0_match*anal_opts.pulsedt; 
mean_cen_time=nanmean(arrayfun(@(shotnum) mean(al_pulses.pos.mean(shotnum,:,1)-al_pulses.time_cen'),1:size(al_pulses.pos.mean,1)));
if abs(mean_cen_time)>tol_t0_match
    est_t0=anal_opts.t0+mean_cen_time;
    fprintf('t0 should be %.6f',est_t0)
end


fprintf('...Done\n') 


end



%% Could check that the pulse time is reasonable usign the fft of the count rate
% 
% tmp_dat=data.mcp_tdc.counts_txy{1};
% 
% opt_in.xdat=tmp_dat(:,1);
% opt_in.max=2.0;
% opt_in.min=0;
% opt_in.bins=1e5;
% opt_in.sigma=1e-4; %timescale that you wish to observe
% opt_in.bin_factor=100;
% out_struct=smooth_hist(opt_in);
% stfig('finding pulse freq','add_stack',1);
% clf
% subplot(2,1,1)
% plot(out_struct.bin.centers,out_struct.count_rate.smooth*1e-3)
% ylabel('count rate (khz)')
% xlabel('time (s)')
% xlim([0.4,0.7])
% subplot(2,1,2)
% 
% fft_out=fft_tx(out_struct.bin.centers,out_struct.count_rate.smooth);
% 
% 
% %%
% dom_opt=[];
% dom_opt.num_components=10;
% [components,details]=dominant_freq_components(out_struct.bin.centers,out_struct.count_rate.smooth,dom_opt);
% %%
% subplot(2,1,2)
% cla
% semilogy(details.fft_dat(1,:),abs(details.fft_dat(2,:)))
% hold on
% semilogy(fft_out(1,:),abs(fft_out(2,:)))
% xlim([0,1000])
% ylim([1e2,1e6])





%% approx way of finding converting mean position into velocity
% %convert the position data into velocity
% % initalize the output array
% al_pulses.vel_approx.mean=al_pulses.pos.mean*nan;
% % the z/time axis is slightly complicated because
% % z=z(0) + vz(0)t+1/2 g t^2
% % want to solve for v(0) given x(t_fall) the fall distance and gravity(in -ve x) , say bec pos is x=0 at t=0
% % -fall_dist=vz(0)t+1/2 g t^2
% % vz(0)=-fall_dist/t  -1/2 g t 
% % now we dont have the fall time directly, we can find the time difference to the expected pulse time 
% pulse_time_diff=al_pulses.pos.mean(:,:,1)-repmat(transpose(al_pulses.time_cen),size(al_pulses.vel_approx.mean,1),1);
% %isequal(col_vec(al_pulses.vel_approx.mean(1,:,1)),col_vec(al_pulses.pos.mean(1,:,1)'-al_pulses.time_cen))
% % then add add on the fall time
% pulse_fall_time=pulse_time_diff+anal_opts.global.fall_time;
% 
% al_pulses.vel_approx.mean(:,:,1)=-anal_opts.global.fall_dist./pulse_fall_time-(1/2)*-const.g0*pulse_fall_time;
% 
% % then we want to convert the x,y data into velocity using the fall time
% % bacause the TOF changes a little from the oscillation in z we can correct for this
% % lets see how big of a change this is
% % calc how much the TOF changes from the z osc
% %std(pulse_time_diff(1,:,1))/anal_opts.global.fall_time
% %this is a fraction of about 0.5 parts per thousand
% % the aprox method
% %al_pulses.vel_approx.mean(:,:,[2,3])=al_pulses.pos.mean(:,:,[2,3])/anal_opts.global.fall_time;
% % the precise method
% al_pulses.vel_approx.mean(:,:,[2,3])=al_pulses.pos.mean(:,:,[2,3])./repmat(pulse_fall_time,1,1,2);
% 
% % this is still aprox. bc. if you use the 'time in the air'(AL pulse to detector) for each count and convert to velocity in x,y and
% % take the mean velocity you will get a slightly different answer than if you take the spatial mean and
% % divide by the mean time to get the velocity 
% % this should be a difference in weighting of ~ 2* mean(al_pulses.pos.std(1,:,1))/anal_opts.global.fall_time 
% % or about 0.7% 
% % to do this all perfectly i need to move this calculation into the loop