function [laser_ok,laser_sm_details]=is_laser_single_mode(input_struct)
%is_laser_single_mode - uses data from a scanning  Fabry-Perot Interferometer (such as https://www.thorlabs.com/newgrouppage9.cfm?objectgroup_id=859)
% to deterime if a laser is single mode. 
% The code suports a two channel measurment of the photodide voltage to enable better dynamic range.
%
% Syntax:  [laser_ok,laser_sm_details]=is_laser_single_mode(import_opts)
%
%inputs
%       input_struct.pd_voltage [matrix N x [1,2]] photodiode voltages, [uncompressed,compressed]
%       input_struct.times [matrix N x 1] time for each of the sample
%       points, must be monotonicaly increasing
%       input_struct.pzt_voltage [matrix N x 1] peizo voltage for the sfp
%       input_struct.scan_type [string] type of pzt scan 'sawtooth','triangle','sine'
%       input_struct.num_checks=10; %how many places to check that the laser is single mode
%       input_struct.thresh_cmp_peak=20e-3; %theshold on the compressed signal to be considered a peak
%       input_struct.peak_distv_min_pass=4.5;%minimum pzt voltage between peaks for the laser to be considered single mode
%       input_struct.scan_time=14e-3;  %estimate of the sfp scan time,used to set the window and the smoothing
%       input_struct.peak_distance_thresh_cmp_full_min [scalar] distance (in pzt voltage) between peaks in the
%       compressed voltage and uncompressed voltage to be considered new peaks
%       input_struct.pd_amp_min %minimum range of the pd signal to indicate the laser has sufficient power
%       input_struct.plot.all=false;
%       input_struct.plot.failed=false;

%output
%      laser_ok  [logical] is the laser single mode
%      laser_sm_details [struct]
%      laser_sm_details.failure_mode [string] 

    
% Known BUGS/ Possible Improvements
%   -return error if no peaks are found
%   -adaptive determination of scan time
%   -handle sawtooth and triangle wave
%   -adaptively handle if the peak int the pd_cmp vs pd_full is the same peak
%   -adaptive smoothing and thresholding
%   -

%   #design dections/devlog
%   - should use input struct or string-val args



% Author: Bryce Henson
% email: Bryce.Henson[a circle]live.com  %YOU MUST INCLUDE '[matlab][is_laser_single_mode]' in the subject line OR I WILL NOT REPLY
% Last revision:2019-05-02

%------------- BEGIN USER VAR --------------
load('is_sm_start.mat')  %DEV DEV DEV REMOVE FOR USE


pzt_scan_freq_lims=[1,2e3];
scan_num_exta_facor=2;%what factor does the number of scans in the data have to exceed the requested query scans in order to just sample a subset of the scans
peak_distance_thresh_cmp_full_min=0.3;%volts difference full-cmp peak in order to be considered new peak
ptz_filt_factor_deriv=1e-3; %fraction of a scan period, for deteriming turning pt
ptz_filt_factor_pks=1e-3;
pd_filt_factor=1e-3;
peak_thresh=[0,0];
cmp_multiplier_disp=50;


input_struct.plot.all=false ;

%------------- END USER VAR --------------
%------------- BEGIN CODE --------------


if ~isfield(input_struct,'peak_distv_min_pass')
    input_struct.peak_distv_min_pass=1;
end
if ~isfield(input_struct,'pd_amp_min')
    input_struct.pd_amp_min=1;
end


%% parse input options

valid_pzt_scan_types={'sawtooth','triangle','sine'};
if ~ismember( input_struct.scan_type,valid_pzt_scan_types)
    warning('valid .scan_types are %s',strjoin(valid_pzt_scan_types,', '))
elseif isequal(input_struct.scan_type,'sine')
    warning('not yet supported')
end


% try to guess the scan freq

% plot(input_struct.times,input_struct.pzt_voltage)
fft_out=fft_tx(input_struct.times,input_struct.pzt_voltage,'window','hanning');

%mask out the lowest and highest frequencies
idxs_freq_range=fast_sorted_mask(fft_out(1,:),pzt_scan_freq_lims(1),pzt_scan_freq_lims(2));
fft_out=fft_out(:,idxs_freq_range(1):idxs_freq_range(2));
%figure(5)
%semilogy(fft_out(1,:),abs(fft_out(2,:)))

%get the frequency with the max powe
[~,idx_max_power_freq]=max(abs(fft_out(2,:)));
pzt_scan_freq=fft_out(1,idx_max_power_freq);
pzt_scan_period=1/pzt_scan_freq;


% use every scan or just sample
% determine roughly how many pzt scans are in the data
% if there are may more than input_struct.num_checks then just sample a
% subset of the scans
%otherwise split up the scans and interogate all of them

time_range=input_struct.times(1)- input_struct.times(end);
scans_in_data=time_range/pzt_scan_period;

% if the scan type is triangle or sine then there are twice the number of
% scans if positive and negative polarity is taken
if ismember(input_struct.scan_type,{'triangle','sine'})
    scans_in_data=scans_in_data*2;
end

use_subset_sample=scans_in_data>input_struct.num_checks*scan_num_exta_facor;

%% split up the scans and feed it to the is_single_scan_sm

if use_subset_sample
    if isinf(input_struct.num_checks)
        error('runtime error infinite number of checks but use_subset_sample=true')
    end
    test_times=linspace(input_struct.times(1),input_struct.times(end),input_struct.num_checks+2); 
    test_times=test_times(2:end-1);
    iimax=numel(test_times);
    %probe a scan centered arround test_time(ii)
    for ii=1:iimax
        time_cen=test_times(ii);
        % to guarantee one good full scan we need to select 2 (and a little bit) scan periods
        time_lims=time_cen+[-1,1]*pzt_scan_period; 
        
        
    end
    
else
    %smooth the pzt data
    pzt_data_filt=gaussfilt(input_struct.times,input_struct.pzt_voltage,ptz_filt_factor_deriv*pzt_scan_period);
    mean_pzt_data_filt=mean(pzt_data_filt);
    pzt_data_deriv=gradient(pzt_data_filt,input_struct.times); %calculate the single sided derivative
    %find the zero crossings of the derivative
    [xing_idx,xing_t]=crossing(pzt_data_deriv,input_struct.times,0,[]);
    xing_idx=xing_idx';
    xing_t=xing_t';
    %determine if these crossings are above or below the mean value with the assumption that if the zero derivative pt
    %is above the mean then the voltage is turning arround to be a negative slope
    xing_positive=pzt_data_filt(xing_idx)<mean_pzt_data_filt;
    
    % alternate method for determing sign of slope
    %pzt_data_2nd_deriv=gradient(pzt_data_deriv,input_struct.times);
    %determine if these crossings they are turning from negative  increasing or decreasing
    %xing_positive=pzt_data_2nd_deriv(xing_idx)>0;
    %
%     clf
%     plot(input_struct.times,pzt_data_filt)
%     hold on
%     plot(xing_t(xing_positive),xing_t(xing_positive)*0+6,'x')
%     plot(xing_t(~xing_positive),xing_t(~xing_positive)*0+19,'x')

    %now determine the phase of the sawtooth from the difference of the timings
    xing_tdiff=diff(xing_t);
    xing_odd_time_diff=mean(xing_tdiff(1:2:end));
    xing_even_time_diff=mean(xing_tdiff(2:2:end));
    if max(xing_odd_time_diff,xing_even_time_diff)<min(xing_odd_time_diff,xing_even_time_diff)*10
        error('based on the timing differences between derivaive=0 points you dont have a sawtooth (or the xing detection didnt work)')
    end
    %this finds if the even or odd diferences between crossings is longer
    sawtooth_phase=xing_odd_time_diff>xing_even_time_diff;
    %calulate if the long part of the sawtooth is positive(true) or negative(false) slope
    sawtooth_polaity=pzt_data_filt(xing_idx(sawtooth_phase+2))<pzt_data_filt(xing_idx(sawtooth_phase+3));
    % this is equivelent to this
    %pzt_val_interp(xing_t(sawtooth_phase+2)) < pzt_val_interp(xing_t(sawtooth_phase+3))
    
    %make a matrix with start and stop indicies for all scans (for the long slope)
    seg_start_idx=xing_idx(sawtooth_polaity+1:2:end);
    seg_stop_idx=xing_idx(sawtooth_polaity+2:2:end);
    %trim extra components in the start index
    seg_start_idx=seg_start_idx(1:numel(seg_stop_idx)); 
    if ~isequal(size(seg_start_idx),size(seg_stop_idx)), error('runtime error start stop idx vec for segment not the same size'), end
    
    %TODO:now adjust these segments to trim off the edges
    plot(input_struct.times(seg_start_idx(1):seg_stop_idx(1)),...
        pzt_data_filt(seg_start_idx(1):seg_stop_idx(1)))
    hold on
    plot(input_struct.times(seg_start_idx(1):seg_stop_idx(1)),...
        input_struct.pd_voltage(seg_start_idx(1):seg_stop_idx(1),2))
    hold off
    
    scans=[];
    single_opts=[];
    single_opts.peak_thresh=peak_thresh;
    single_opts.pd_filt_factor=pd_filt_factor;
    single_opts.pzt_scan_period=pzt_scan_period;
    single_opts.pzt_filt_factor=ptz_filt_factor_pks;
    single_opts.merge_cmp_full_pk_dist=peak_distance_thresh_cmp_full_min;
    single_opts.peak_distv_min_pass=input_struct.peak_distv_min_pass;
    single_opts.pd_amp_min=input_struct.pd_amp_min;
    single_opts.plot=input_struct.plot;

    single_opts.cmp_multiplier_disp=cmp_multiplier_disp;
    pk_spacing_list=[];
    iimax=numel(seg_start_idx);
    for ii=1:iimax
        
        single_opts.pd_voltage=input_struct.pd_voltage(seg_start_idx(ii):seg_stop_idx(ii),:);
        single_opts.time=input_struct.times(seg_start_idx(ii):seg_stop_idx(ii));
        single_opts.pzt_voltage=input_struct.pzt_voltage(seg_start_idx(ii):seg_stop_idx(ii));
        scans{ii}=is_single_scan_sm(single_opts);
        scans{ii}.start_idx=seg_start_idx(ii);
        %find the next falling edge
        scans{ii}.stop_idx=seg_stop_idx(ii);
        pk_spacing_list=[pk_spacing_list;diff(scans{ii}.pks.full.pzt)];
    end


end



%%

%%If the check if the laser is single mode



% test_sm_while=true; %intialize while loop flag, give up with one bad detection
% jj=0;
% 
% sweep=[];
% single_mode_vec=false(1,jjmax);
% 
% jj=0;
% while test_sm_while
%     jj=jj+1;
%     time_start=test_times(jj);
%     probe_sampl_start=max(1,1+floor(time_start*sr));
%     probe_sampl_stop=min(samples,probe_sampl_start+ceil(args_single.window_time*sr));
%     if probe_sampl_stop-probe_sampl_start>=args_single.window_time*sr %check that we have enough points to work with
%         
%         sub_dat_grad_smooth=diff(sub_dat_smooth)*sr;
%         pos_slope_mask=sub_dat_grad_smooth>0;
%         %now we want to select the first positive going sgement that is long enough
%         x=pos_slope_mask';
%         % Find the transitions
%         down= strfind(x',[1 0]);
%         ups= strfind(x',[0 1]);
%         %we take the ramp to investigate as the between the first rising edge of the mask and the next
%         %falling edge
%         sweep{jj}.start_idx=ups(1);
%         %find the next falling edge
%         sweep{jj}.stop_idx=down(find(down>sweep{jj}.start_idx,1)); %fuck so elegant, fuck this is so much beter than LabView
%         sweep{jj}.pzt_raw=sub_ptz_raw(sweep{jj}.start_idx:sweep{jj}.stop_idx);
%         sweep{jj}.pzt_smooth=sub_dat_smooth(sweep{jj}.start_idx:sweep{jj}.stop_idx);
%         sweep{jj}.pd_full_raw=sfp_pd_raw(probe_sampl_start+sweep{jj}.start_idx:probe_sampl_start+sweep{jj}.stop_idx);
%         sweep{jj}.pd_cmp_raw=ai_dat.Data(4,probe_sampl_start+sweep{jj}.start_idx:probe_sampl_start+sweep{jj}.stop_idx);
%         sweep{jj}.pd_cmp_raw=sweep{jj}.pd_cmp_raw-median(sweep{jj}.pd_cmp_raw);
%         %threshold out all data that is 7% the max value for the peak finding thresh
%         thesh=0.07*max(sweep{jj}.pd_full_raw);
%         [pks_val,locs] = findpeaks(sweep{jj}.pd_full_raw,'MinPeakHeight',thesh);
%         sweep{jj}.pks.full.pd=pks_val;
%         sweep{jj}.pks.full.pzt=sweep{jj}.pzt_smooth(locs);
%         [pks_val,locs] = findpeaks(sweep{jj}.pd_cmp_raw,'MinPeakHeight',args_single.sfp.thresh_cmp_peak);
%         sweep{jj}.pks.cmp.pd=pks_val;
%         sweep{jj}.pks.cmp.pzt=sweep{jj}.pzt_smooth(locs);
%         %find the distance in pzt voltage to the nearest full peak
%         min_dist_cmp_to_full=arrayfun(@(x) min(abs(x-sweep{jj}.pks.full.pzt)),sweep{jj}.pks.cmp.pzt);
%         %peak_distance_thresh_cmp_full_min
%         sweep{jj}.pks.cmp.pzt=sweep{jj}.pks.cmp.pzt(min_dist_cmp_to_full>peak_distance_thresh_cmp_full_min);
%         sweep{jj}.pks.cmp.pd=sweep{jj}.pks.cmp.pd(min_dist_cmp_to_full>peak_distance_thresh_cmp_full_min);
%         sweep{jj}.pks.all.pzt=[sweep{jj}.pks.cmp.pzt,sweep{jj}.pks.full.pzt];
%         sweep{jj}.pks.all.pd=[sweep{jj}.pks.cmp.pd,sweep{jj}.pks.full.pd];
%         sweep{jj}.pks.all.min_pzt_v_diff=min_diff(sweep{jj}.pks.all.pzt);
% 
%         if sweep{jj}.pks.all.min_pzt_v_diff<args_single.sfp.peak_dist_min_pass
%             fprintf('\nLASER IS NOT SINGLE MODE!!!\n%04i',0')
%             single_mode_vec(jj)=false;%give up if anything looks bad
%         else
%            single_mode_vec(jj)=true;
% 
%         end
% 
%         if args_single.plot.all || (args_single.plot.failed && ~single_mode_vec(jj))
%             title_strs={'Failed','OK'};
%             sfigure(1);
%             subplot(2,2,1)
%             cla;
%              set(gcf,'color','w')
%             title('smoothing pzt v')
%             time=(1:numel(sub_ptz_raw))/sr;
%             plot(time*1e3,sub_ptz_raw,'r')
%             hold on
%             plot(time*1e3,sub_dat_smooth,'k')
%             hold off
%             xlabel('time (ms)')
%             ylabel('volts (v)')
% 
%             subplot(2,2,2)
%             cla;
%             %select a single scan
%             %to do this we will mask out the positive slope
%             diff_sub_raw=diff(sub_ptz_raw)*sr;
%             diff_time=time(1:end-1)+0.5*(time(2)-time(1));
%             plot(diff_time*1e3,diff_sub_raw)
%             hold on
%             plot(diff_time*1e3,sub_dat_grad_smooth)
%             plot(diff_time*1e3,pos_slope_mask*2e3)
%             hold off
%             ylim([min(diff_sub_raw),max(diff_sub_raw)])
%             xlabel('time (ms)')
%             ylabel('grad (v/s)')
% 
% 
%             subplot(2,2,3)
%             cla;
%             plot(sweep{jj}.pzt_smooth,sweep{jj}.pd_full_raw,'k')
%             xlabel('pzt(v)')
%             ylabel('pd (v)')
%             hold on
%             plot(sweep{jj}.pks.full.pzt,sweep{jj}.pks.full.pd,'xk','markersize',20)
%             plot(sweep{jj}.pzt_smooth,sweep{jj}.pd_cmp_raw*args_single.cmp_multiplier_disp,'r')
%             plot(sweep{jj}.pks.cmp.pzt,sweep{jj}.pks.cmp.pd*args_single.cmp_multiplier_disp,'rx','markersize',20);
%             hold off
%             title(title_strs(single_mode_vec(jj)+1))
% 
%             pause(1e-6)
%         end %end plots
%     end %enough points to work with
% 
% if jj>=jjmax || ~single_mode_vec(jj) %if something is not single mode do not continue
%     test_sm_while=false;
% end 
% end%end while loop over checking mode

end



