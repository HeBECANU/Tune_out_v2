function al_pulses=segment_al_pulses_BEC_therm(anal_opts,data)

num_shots=size(data.mcp_tdc.counts_txy,2);
al_pulses=[];
al_pulses.pulsedt=anal_opts.pulsedt;
al_pulses.window=nan(anal_opts.pulses,3,2); %initalize
al_pulses.num_counts=nan(num_shots,anal_opts.pulses);
al_pulses.pos.mean=nan(num_shots,anal_opts.pulses,3);
al_pulses.pos.std=nan(num_shots,anal_opts.pulses,3);
al_pulses.vel.mean=nan(num_shots,anal_opts.pulses,3);
al_pulses.vel.std=nan(num_shots,anal_opts.pulses,3);
%al_pulses.seg_thresh=6e6; % in countes per cubic meter
seg_thresh=anal_opts.seg_thresh;

% internal options
filt_gaussian_factor=5; % how big to make the filter in number of sigma, recomend between 4 and 7 does not seem to have a big impact on speed

global const %get gravity

% set up the edges for the 3d histogram
edge_vecs=cell(1,3);
% first we have to find what pulse_twindows is in velocity
% we can get an approx expression with
% vel_range_z=anal_opts.global.fall_velocity*anal_opts.pulse_twindow*[-0.5,0.5]
% now we have the z axis converted to a distance, to convert to a velocity use
% vel_range_z = vel_range_z/anal_opts.global.fall_time

% for completeness ill use txy_to_vel
mock_pts=zeros(2,3);
mock_pts(:,1)=anal_opts.pulse_twindow*[-0.5,0.5];
mock_pts(:,1)=anal_opts.pulse_twindow*[-0.5,0.5];
anal_opts.xylim
vzxy_out=txy_to_vel(mock_pts,...
                            -anal_opts.global.fall_time,...
                            -const.g0,...
                            anal_opts.global.fall_dist);
% now strictly this isnt excatly right as the square in txy will not be a square in velocity
% but the difference is tiny for the small time range here "good enuff for straya"
xyrange_vel=anal_opts.xylim./repmat(anal_opts.global.fall_time,1,2);
range_vel_xyz=cat(1,xyrange_vel,vzxy_out(:,1)');

bins_zxy=ceil(diff(range_vel_xyz,1,2)/anal_opts.bin_size);
for ii=1:numel(edge_vecs)
   edge_vecs{ii}=linspace(range_vel_xyz(ii,1),range_vel_xyz(ii,2),bins_zxy(ii));
end


output_chars=fprintf('binning pulse %04u of %04u in file %04u of %04u',0,anal_opts.pulses,0,num_shots);
first_good_shot=true;
for shot=1:num_shots
        if data.mcp_tdc.all_ok(shot)
            for pulse=1:anal_opts.pulses
                %set up time window centered arround t0
                t_pulse_cen=anal_opts.t0+anal_opts.pulsedt...
                    *(anal_opts.start_pulse+pulse-2);
                trange=t_pulse_cen+anal_opts.pulse_twindow*[-0.5,0.5];
                pulse_win_txy=[trange;anal_opts.xylim]; 
                counts_pulse=masktxy_square(data.mcp_tdc.counts_txy{shot},pulse_win_txy);
                
                
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
                vxyz_out=vzxy_out;
                vxyz_out(:,1)=vzxy_out(:,2);
                vxyz_out(:,2)=vzxy_out(:,3);
                vxyz_out(:,3)=vzxy_out(:,1);
                %% find some stats on the vel
                al_pulses.unsegmented_vel_xyz.mean(shot,pulse,:)=mean(vxyz_out,1);
                al_pulses.unsegmented_vel_xyz.std(shot,pulse,:)=std(vxyz_out,1);
                al_pulses.unsegmented_vel_xyz.num_counts(shot,pulse)=size(vxyz_out,1);
                al_pulses.unsegmented_vel_xyz.ste(shot,pulse,:)=std(vxyz_out,1)/sqrt(size(vxyz_out,1)); 
                al_pulses.unsegmented_vel_xyz.counts{shot,pulse}=vxyz_out;
                
                
                %% make a 3d histogram
                [pulse_hist,~,vol_hist_midpts,count_bin_idx]=histcn(vxyz_out,edge_vecs{:});  
                pulse_hist=pulse_hist./(anal_opts.bin_size^3); %normalize to flux
                
                %now blur this histogram
                sigma_bins=anal_opts.blur_size/anal_opts.bin_size;
                filt_size=sigma_bins*filt_gaussian_factor; % to well represnt the gaussian
                % round to nearest odd number for imgaussfilt3
                filt_size=floor(filt_size/2)*2+1;
                pulse_hist_smooth=imgaussfilt3(pulse_hist,sigma_bins,...
                                                'FilterSize',filt_size); 
                mesh_cord=[];
                [mesh_cord.x,mesh_cord.y,mesh_cord.z]=meshgrid(vol_hist_midpts{1},...
                                                                vol_hist_midpts{2},...
                                                                vol_hist_midpts{3});                            
                
                %% make a slice plot of the 3d density distribution
                if anal_opts.plot.all
                    stfig('bin distribution');
                    sdat=smooth_hist(pulse_hist(:),'sigma',5e7);
                    plot(sdat.bin.centers,sdat.counts.smooth)
                    %set(gca,'Xscale','log')
                    set(gca,'Yscale','log') 
                    
    %                 zslice = linspace(zxyrange_vel(1,1),zxyrange_vel(1,2),3);
    %                 xslice = linspace(zxyrange_vel(2,1),zxyrange_vel(2,2),3);   
    %                 yslice = linspace(zxyrange_vel(3,1),zxyrange_vel(3,2),3); 
                    xslice = [-10]*1e-3 +al_pulses.unsegmented_vel_xyz.mean(shot,pulse,1);
                    yslice = [-10]*1e-3 +al_pulses.unsegmented_vel_xyz.mean(shot,pulse,2);
                    zslice = [-10]*1e-3 +al_pulses.unsegmented_vel_xyz.mean(shot,pulse,3);
                    stfig('3d slice plot');
                    clf
                    h=slice(mesh_cord.x,mesh_cord.y,mesh_cord.z,pulse_hist_smooth,zslice,xslice,yslice);
                    set(h,'edgecolor','none')
                    colormap(plasma)
                    alpha(0.5);       
                    colorbar
    %                 colormap(flipud(gray))
    %                 alpha('color');
                    xlabel('x')
                    ylabel('z')
                    zlabel('y')
                end
                %% 3d contor plot
                    if anal_opts.plot.all
                    stfig('3d contor plot');
                    clf
                    max_hist=max(pulse_hist_smooth(:));
                    flux_steps=max_hist*linspace(0.05,0.9,10);

                    hold on
                    kkmax=numel(flux_steps);
                    colors=viridis(kkmax);
                    face_alphas=linspace(0.1,0.9,kkmax)/(0.2*kkmax);
                    for kk=1:kkmax
                        patch(isosurface(mesh_cord.x,mesh_cord.y,mesh_cord.z,pulse_hist_smooth,flux_steps(kk)),'FaceColor',colors(kk,:),'EdgeColor','none','FaceAlpha',face_alphas(kk))
                    end
                    hold off
                    xlim(range_vel_xyz(1,:))
                    ylim(range_vel_xyz(2,:))
                    zlim(range_vel_xyz(3,:))
                    xlabel('x')
                    ylabel('y')
                    zlabel('z')
                    daspect([1,1,1])
                end
                
                %% segmentation
                
                if anal_opts.norm_seg_thresh
                    this_seg_thresh=seg_thresh*size(vxyz_out,1);
                else
                    this_seg_thresh=seg_thresh;
                end
                
                segmentation_mask=pulse_hist_smooth>this_seg_thresh;
                
                %% plot the segmentation
                if  anal_opts.plot.all
                    stfig('3d contor plot');
                    hold on
                    this_face_alpha=0.5;
                    patch(isosurface(mesh_cord.x,mesh_cord.y,mesh_cord.z,segmentation_mask,0.5),'FaceColor',[0.1,0.1,0.1],'EdgeColor','none','FaceAlpha',this_face_alpha)
                    hold off 
                end
                
                %% try filling the segmentation so its a disk not a donut
                
                % old way finding max radius from mean which does not give minimum bounding circle
%                 mean_seg_xyz=nan(1,3);
%                 segmentation_mask_nan=segmentation_mask*1;
%                 segmentation_mask_nan(segmentation_mask_nan==0)=nan;
%                 tmp_val_x=segmentation_mask_nan.*mesh_cord.x;
%                 mean_seg_xyz(1)=nanmean(tmp_val_x(:));
%                 tmp_val_y=segmentation_mask_nan.*mesh_cord.y;
%                 mean_seg_xyz(2)=nanmean(tmp_val_y(:));
%                 tmp_val_z=segmentation_mask_nan.*mesh_cord.z;
%                 mean_seg_xyz(3)=nanmean(tmp_val_z(:));
%                 % now that we have the mean position finx the max z y radius and 
%                 radius_xz= sqrt( (tmp_val_x-mean_seg_xyz(1)).^2 + (tmp_val_z-mean_seg_xyz(3)).^2 );
%                 max_radius_zy=max(radius_xz(:));
%                 tmp_y_delt=abs(tmp_val_y-mean_seg_xyz(2));
%                 max_y_delt= max(tmp_y_delt(:))
%                 segmentation_mask_cyl= sqrt( (mesh_cord.x-mean_seg_xyz(1)).^2 + ...
%                                                 (mesh_cord.z-mean_seg_xyz(3)).^2 ) < max_radius_zy  & ...
%                                         abs(abs(mesh_cord.y-mean_seg_xyz(2))) < max_y_delt;

                %% find the mean postion of the segmentation and then the max radius andy cord from the mean
                % this was used to find the disk radius before but it is not a minimum covering circle
                mean_seg_xyz=nan(1,3);
                segmentation_mask_nan=segmentation_mask*1;
                segmentation_mask_nan(segmentation_mask_nan==0)=nan;
                tmp_val_x=segmentation_mask_nan.*mesh_cord.x;
                mean_seg_xyz(1)=nanmean(tmp_val_x(:));
                tmp_val_y=segmentation_mask_nan.*mesh_cord.y;
                mean_seg_xyz(2)=nanmean(tmp_val_y(:));
                tmp_val_z=segmentation_mask_nan.*mesh_cord.z;
                mean_seg_xyz(3)=nanmean(tmp_val_z(:));
                %radius_xz= sqrt( (tmp_val_x-mean_seg_xyz(1)).^2 + (tmp_val_z-mean_seg_xyz(3)).^2 );
                %max_radius_zy=max(radius_xz(:));
                tmp_x_delt=abs(tmp_val_y-mean_seg_xyz(2));
                max_x_delt= max(tmp_x_delt(:));
                
                % output some of these details
                al_pulses.seg_details.flux_thresh.mean(shot,pulse,:)=mean_seg_xyz;
                
                %% turn the threshold segmentation which looks like a donut into a disk
                % we will find the minimum covering circle
                % fist collapse allong the y axis of the logical segmentation_mask
                mask_shadow_xz=squeeze(any(segmentation_mask,1));
                
                % find the cordinates of the edge positions
                mask_edges_xz=edge(mask_shadow_xz);
                [edges_idx_x, edges_idx_z] = find(mask_edges_xz ) ;
                edges_cord=[];
                edges_cord(:,1)=vol_hist_midpts{2}(edges_idx_x);
                edges_cord(:,2)=vol_hist_midpts{3}(edges_idx_z);
                % a diagnostic plot showing that the edge cord were found right
                if  anal_opts.plot.all
                    stfig('edges');
                    imagesc(mask_shadow_xz')
                    set(gca,'YDir','normal')
                    hold on
                    plot(edges_idx_x,edges_idx_z,'x')
                    hold off
                end
                % now we have a list of xy cord of the edges we need to find the minimum circle that covers these pts
                [center_xz,radius] = minboundcircle(edges_cord(:,1),edges_cord(:,2));
                % now remake the mask with the disk
                
                al_pulses.seg_details.min_cyl.cen(shot,pulse,:)=cat(1,center_xz(1),mean_seg_xyz(2),center_xz(2));
                al_pulses.seg_details.min_cyl.radius(shot,pulse)=radius;
                al_pulses.seg_details.min_cyl.ydelt(shot,pulse)=radius;
                
                segmentation_mask_cyl= sqrt( (mesh_cord.x-center_xz(1)).^2 + ...
                                                (mesh_cord.z-center_xz(2)).^2 ) < radius  & ...
                                        abs(abs(mesh_cord.y-mean_seg_xyz(2))) < max_x_delt;
                % add this to the 3d contor plot
                if  anal_opts.plot.all
                    stfig('3d contor plot');
                    hold on
                    plot3(mean_seg_xyz(1),mean_seg_xyz(2),mean_seg_xyz(3),'or')
                    patch(isosurface(mesh_cord.x,mesh_cord.y,mesh_cord.z,segmentation_mask_cyl,0.5),'FaceColor',[0.5,0.1,0.1],'EdgeColor','none','FaceAlpha',this_face_alpha)
                    hold off
                    view(20,20)
                    title('min span cyl sementation')
                end
                
                
                %% save the counts and the means
                % determine if a given count is in the segmentation mask using the bin index that histcn so nicely
                % returned
                count_in_high_flux=NaN(size(vzxy_out,1),1);
                for count=1:size(vzxy_out,1)
                    this_index=count_bin_idx(count,:);
                    count_in_high_flux(count)=segmentation_mask_cyl(this_index(1),this_index(2),this_index(3));
                end
                if any(isnan(count_in_high_flux)), error('did not set if count is high flux'), end
                
                
                vzxy_high_flux=vzxy_out(logical(count_in_high_flux),:);
                vzxy_low_flux=vzxy_out(~logical(count_in_high_flux),:);
                
                al_pulses.segmented_vel_xyz.high_flux.mean(shot,pulse,:)=mean(vzxy_high_flux,1);  
                al_pulses.segmented_vel_xyz.high_flux.std(shot,pulse,:)=std(vzxy_high_flux,1); 
                al_pulses.segmented_vel_xyz.high_flux.num_counts(shot,pulse)=size(vzxy_high_flux,1);
                al_pulses.segmented_vel_xyz.high_flux.ste(shot,pulse,:)=std(vzxy_high_flux,1)/sqrt(size(vzxy_high_flux,1)); 
                al_pulses.segmented_vel_xyz.high_flux.counts{shot,pulse}=vzxy_high_flux; 
               
                al_pulses.segmented_vel_xyz.low_flux.mean(shot,pulse,:)=mean(vzxy_low_flux,1);
                al_pulses.segmented_vel_xyz.low_flux.std(shot,pulse,:)=std(vzxy_low_flux,1);
                al_pulses.segmented_vel_xyz.low_flux.num_counts(shot,pulse)=size(vzxy_low_flux,1);
                al_pulses.segmented_vel_xyz.low_flux.ste(shot,pulse,:)=std(vzxy_low_flux,1)/sqrt(size(vzxy_low_flux,1));
                al_pulses.segmented_vel_xyz.low_flux.counts{shot,pulse}=vzxy_low_flux; 
                
                % lets find the counts that are in the low flux segmentation but are inside anal_opts.close_thresh
                %simple radius from origin
                %mask=vecnorm(vzxy_low_flux,2,2)<anal_opts.close_thresh;
                % radius from cen of high flux
                radius=vecnorm(vzxy_low_flux- repmat(mean(vzxy_high_flux,1),size(vzxy_low_flux,1),1) ...
                            ,2,2);
                mask=radius>anal_opts.close_thresh(1) &  radius<anal_opts.close_thresh(2);
                
                vxyz_close_low_flux=vzxy_low_flux(mask,:);
                al_pulses.segmented_vel_xyz.close_low_flux.mean(shot,pulse,:)=mean(vxyz_close_low_flux,1);
                al_pulses.segmented_vel_xyz.close_low_flux.std(shot,pulse,:)=std(vxyz_close_low_flux,1);
                al_pulses.segmented_vel_xyz.close_low_flux.num_counts=size(vxyz_close_low_flux,1);
                al_pulses.segmented_vel_xyz.close_low_flux.ste(shot,pulse,:)=std(vxyz_close_low_flux,1)/sqrt(size(vxyz_close_low_flux,1));
                al_pulses.segmented_vel_xyz.close_low_flux.counts{shot,pulse}=vxyz_close_low_flux; 
                
          
                if anal_opts.plot.all ||  mod(pulse,10)==0
                    fprintf(repmat('\b',[1,output_chars]))
                    output_chars=fprintf('binning pulse %04u of %04u in file %04u of %04u',pulse,anal_opts.pulses,shot,num_shots);
                end
                
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

% 
% 
% %function out=fit_temperature(anal_opts,data)
% %anal_opts.xylim=anal_opts.tdc_import.txylim(2:3,:);
% 
% % lets find the temperature of the BEC usign the AL pulse data
% % unfortunately the thermal cloud will not be in phase with the BEC part of the AL pulse
% % this is bc the thermal and BEC component have different trap frequencies
% % and bc the gaussian fit width using a simple masking procedure (fit_temperature_simple) is of a similar size ~10mm/s
% % to the oscillation amplitude ~10mm/s
% 
% %what we need to do then is to sgement each AL pulse into BEC component and thermal component
% % to do this we will use a 3d histogram, smooth it then threshold it
% 
% 
% % this will measure the termperature by taking a cylinder of counts in the weak axis from every atom laser pulse
% cyl_radius=101e-3;
% cen_type='mean';  %could be mean,window or thermal 
% if strcmp(cen_type,'mean')
%     vel_center=data.mcp_tdc.al_pulses.vel_zxy.mean; %file,pulse_num,axis) 
% 
% elseif strcmp(cen_type,'window')
%     vel_center=mean(data.mcp_tdc.al_pulses.window,3);
%     % make it the same size as above  file,pulse_num,axis
%     vel_center = permute(vel_center,[3 1 2]);
%     vel_center=repmat(vel_center,size(data.mcp_tdc.al_pulses.pos.mean,1),1,1);
% end
% 
% mask_opts=[];
% mask_opts.type='allow cyl';
% mask_opts.radius=5e-3/anal_opts.global.fall_time;
% 
% % save all the data that we get from each pulse
% windowed_data=cell(size(data.mcp_tdc.al_pulses.pos.mean(:,:,1)));
% iimax=size(data.mcp_tdc.al_pulses.pos.mean,1);
% jjmax=size(data.mcp_tdc.al_pulses.pos.mean,2);
% for ii=1:size(data.mcp_tdc.al_pulses.pos.mean,1)
%     fprintf('%u\n',ii)
%     file_counts=data.mcp_tdc.counts_txy{ii};
%     for jj=1:jjmax
%         t_window=data.mcp_tdc.al_pulses.window(jj,1,:);
%         t_pulse_cen= data.mcp_tdc.al_pulses.time_cen(jj);
%         mask_idx=fast_sorted_mask(file_counts(:,1),...
%             t_window(1),t_window(2));
%         window_data_this_txy=file_counts(mask_idx(1):mask_idx(2),:);
%         
%         vzxy_out=txy_to_vel(window_data_this_txy,...
%                                     t_pulse_cen-anal_opts.global.fall_time,...
%                                     -const.g0,...
%                                     anal_opts.global.fall_dist);
%                                 
%         
%         %mean(window_data_this,1)
%         cen_vel_this=squeeze(vel_center(ii,jj,:));
%         cen_vel_this=repmat(permute(cen_vel_this,[2,1]),size(vzxy_out,1),1);
%         vzxy_out=vzxy_out-cen_vel_this;
%          %mean(window_data_this,1)
%         %fprintf('%u\n',size(window_data_this,1))
%         
%         if strcmp(mask_opts.type,'allow cyl')
%             radial_yt=rms(vzxy_out(:,[1,3]),2);
%             %min(radial_yt)
%             mask=radial_yt<mask_opts.radius;
%             vzxy_out=vzxy_out(mask,:);
%             windowed_data{ii,jj}=vzxy_out;
%         end
%         
%     end
% end
% 
% num_norm_pts=2e8;
% norm_counts=(rand(num_norm_pts,3)-0.5)*2;
% random_range=32e-3/anal_opts.global.fall_time;
% norm_counts=norm_counts*random_range;
% if strcmp(mask_opts.type,'allow cyl')
%     radial_yt=rms(norm_counts(:,[1,3]),2);
%     %min(radial_yt)
%     mask=radial_yt<mask_opts.radius;
%     norm_windowed=norm_counts(mask,:);
% end
% 
%    
% 
% %% combine all data
% all_windows_comb=cat(1,windowed_data{:,1});
% radius= rms(all_windows_comb,2);
% radius_norm= rms(norm_windowed,2);
% 
% 
% 
% stfig('themal tails')
% clf
% sigma=5e-4;
% sdat=smooth_hist(radius,'sigma',sigma);
% s_ref=smooth_hist(radius_norm,'sigma',sigma,'edges',sdat.bin.edge);
% 
% % bin_volume=(4/3)*pi* (sdat.bin.edge(2:end).^3  - sdat.bin.edge(1:end-1).^3) ;
% % flux=./bin_volume;
% 
% flux=sdat.count_rate.smooth./s_ref.count_rate.smooth_prob;
% %flux=sdat.count_rate.smooth
% %flux=s_ref.count_rate.smooth_prob
% 
% 
% %x_mask_lims=[5e-3,15e-3]./anal_opts.global.fall_time;
% x_mask_lims=[5e-3,20e-3]./anal_opts.global.fall_time;
% 
% xmask= sdat.bin.centers>x_mask_lims(1)  & sdat.bin.centers<x_mask_lims(2);
% xdat=sdat.bin.centers(xmask);
% ydat=flux(xmask) *1e-6;
% 
% plot(xdat,ydat )
% %set(gca,'Yscale','log')
% %set(gca,'Xscale','log')
% 
% %
% 
% 
% opt = statset('TolFun',1e-10,'TolX',1e-10,...
%     'MaxIter',1e4,... %1e4
%     'UseParallel',1);
% cof_names={'sigma','amp','offset'};
% modelfun=@(b,x) gaussian_function_1d(x,b(1),0,b(2),b(3))
% beta0=[0.02,1,0];
% 
% % cof_names={'sigma','amp'};
% % modelfun=@(b,x) gaussian_function_1d(x,b(1),0,b(2))
% % beta0=[0.02,1];
% 
% fitobject=fitnlm(xdat,ydat,modelfun,beta0,...
%     'options',opt,...
%     'CoefficientNames',cof_names);
% 
% xsamp=linspace(x_mask_lims(1),x_mask_lims(2),1e3)';
% [prediction,ci]=predict(fitobject,xsamp,'Alpha',1-erf(1/sqrt(2))); %,'Prediction','observation'
%   hold on
%   
% color_shaded=[1,1,1]*0.9;
% colors_main=[0.5,0.6,0.2];
% patch([xsamp', fliplr(xsamp')], [ci(:,1)', fliplr(ci(:,2)')], color_shaded,'EdgeColor','none')  %[1,1,1]*0.80
%               
% plot(xsamp(:,1),prediction,'-','LineWidth',1.0,'Color',colors_main)
% hold off
% 
% fit_temperature=(fitobject.Coefficients.Estimate(1).^2) *const.mhe /(const.kb*3);
% fit_temperature
% 
% 
% %%
% all_windows_comb=cat(1,windowed_data{:,100});
% % radius= rms(all_windows_comb,2);
% % radius_norm= rms(norm_windowed,2);
% 
% radius= all_windows_comb(:,2);
% %radius_norm= norm_windowed(:,2);
% 
% 
% stfig('themal tails')
% clf
% sigma=5e-4;
% sdat=smooth_hist(radius,'sigma',sigma);
% %s_ref=smooth_hist(radius_norm,'sigma',sigma,'edges',sdat.bin.edge);
% 
% % bin_volume=(4/3)*pi* (sdat.bin.edge(2:end).^3  - sdat.bin.edge(1:end-1).^3) ;
% % flux=./bin_volume;
% 
% %flux=sdat.count_rate.smooth./s_ref.count_rate.smooth_prob;
% flux=sdat.count_rate.smooth
% %flux=s_ref.count_rate.smooth_prob
% 
% 
% %x_mask_lims=[5e-3,15e-3]./anal_opts.global.fall_time;
% x_mask_lims=[-30e-3,30e-3]./anal_opts.global.fall_time;
% xmask= sdat.bin.centers>x_mask_lims(1)  & sdat.bin.centers<x_mask_lims(2);
% 
% x_mask_lims_inner=[-1,1]*8e-3./anal_opts.global.fall_time;
% xmask=xmask & ~(sdat.bin.centers>x_mask_lims_inner(1)  & sdat.bin.centers<x_mask_lims_inner(2) );
% sum(xmask)
% xdat=sdat.bin.centers(xmask);
% ydat=flux(xmask) *1e-6;
% 
% plot(xdat,ydat )
% %set(gca,'Yscale','log')
% %set(gca,'Xscale','log')
% 
% %
% 
% 
% opt = statset('TolFun',1e-10,'TolX',1e-10,...
%     'MaxIter',1e4,... %1e4
%     'UseParallel',1);
% cof_names={'sigma','amp','offset'};
% modelfun=@(b,x) gaussian_function_1d(x,b(1),0,b(2),b(3))
% beta0=[0.02,1,0];
% 
% % cof_names={'sigma','amp'};
% % modelfun=@(b,x) gaussian_function_1d(x,b(1),0,b(2))
% % beta0=[0.02,1];
% 
% fitobject=fitnlm(xdat,ydat,modelfun,beta0,...
%     'options',opt,...
%     'CoefficientNames',cof_names);
% 
% xsamp=linspace(x_mask_lims(1),x_mask_lims(2),1e3)';
% [prediction,ci]=predict(fitobject,xsamp,'Alpha',1-erf(1/sqrt(2))); %,'Prediction','observation'
%   hold on
%   
% color_shaded=[1,1,1]*0.9;
% colors_main=[0.5,0.6,0.2];
% patch([xsamp', fliplr(xsamp')], [ci(:,1)', fliplr(ci(:,2)')], color_shaded,'EdgeColor','none')  %[1,1,1]*0.80
%               
% plot(xsamp(:,1),prediction,'-','LineWidth',1.0,'Color',colors_main)
% hold off
% 
% fit_temperature=(fitobject.Coefficients.Estimate(1).^2) *const.mhe /(const.kb*3);
% fit_temperature
% 
% 
% 
% %%
% % for reference what a gaussian looks like
% stfig('gauusian analytical')
% xsamp=linspace(1,5,1e3);
% plot(xsamp,gaussian_function_1d(xsamp,1,0))
% set(gca,'Yscale','log')
% set(gca,'Xscale','log')
% %%
% 

%end