function [laser_ok,laser_sm_details]=is_laser_single_mode(in_struct)
%is_laser_single_mode - uses data from a scanning  Fabry-Perot Interferometer (such as https://www.thorlabs.com/newgrouppage9.cfm?objectgroup_id=859)
% to deterime if a laser is single mode. 
% The code suports a two channel measurment of the photodide voltage to enable better dynamic range.
%
% Syntax:  [laser_ok,laser_sm_details]=is_laser_single_mode(import_opts)
%
%inputs
%       in_struct
%                   .pd_voltage [matrix N x [1,2]] photodiode voltages, [uncompressed,compressed]
%                   .times [matrix N x 1] time for each of the sample
%       points, must be monotonicaly increasing
%                   .pzt_voltage [matrix N x 1] peizo voltage for the sfp
%                   .scan_type [string] type of pzt scan 'sawtooth','triangle','sine'
%                   .num_checks=10; %how many places to check that the laser is single mode
%                   .peak_thresh=[0,0]; %theshold on the [uncompressed,compressed] signal to be considered a peak
%                   .pzt_dist_sm %minimum pzt voltage between peaks for the laser to be considered single mode
%                   .pzt_dist_pd_cmp [scalar] distance (in pzt voltage) between peaks in the
%       compressed voltage and uncompressed voltage to be considered new peaks
%                   .scan_time=14e-3;  %estimate of the sfp scan time,used to set the window and the smoothing
%                   .pd_filt_factor    %fraction of a scan to smooth the pd data by for peak detection
%                   .ptz_filt_factor_pks  %fraction of a scan to smooth the pzt data by for peak detection
%                   .pzt_filt_factor_deriv %fraction of a scan to smooth the data by for derivative detection
%                   .pd_amp_min %minimum range of the pd signal to indicate the laser has sufficient power
%                   .plot.all=false;
%                   .plot.failed=false;
%                   .verbose=1 [int] zero no output bigger numbers more verbose

%output
%      laser_ok  [logical] is the laser single mode
%      laser_sm_details [struct]
%      laser_sm_details.failure_mode [string] 

    
% Known BUGS/ Possible Improvements
%   -return error if no peaks are found
%   -[x] adaptive determination of scan time
%   -handle sawtooth and triangle wave
%   -[x] adaptively handle if the peak int the pd_cmp vs pd_full is the same peak
%   -[x] adaptive smoothing and thresholding
%   -


% Author: Bryce Henson
% email: Bryce.Henson[a circle]live.com  %YOU MUST INCLUDE '[matlab][is_laser_single_mode]' in the subject line OR I WILL NOT REPLY
% Last revision:2019-05-07

%------------- BEGIN USER VAR --------------
%%load('is_sm_start.mat')  %DEV DEV DEV REMOVE FOR USE


pzt_scan_freq_lims=[0.1,10e3]; %this defines the expected range of the pzt scan freqs, not very important
scan_num_exta_facor=2;%what factor does the number of scans in the data have to exceed the requested query scans in order to just sample a subset of the scans
cmp_multiplier_disp=NaN;
st_pt_min_pzt_factor=1; %how big the difference between stationary points needs to be, in factors of the std of the pzt voltage
scan_clip_factor=1e-2; %remove peaks that are not further than this away from the edges, in factors of the std of the pzt voltage

%in_struct.pzt_scan_period=50e-3; %DEV DEV DEV REMOVE FOR USE


in_struct.plot.all=false ;

%------------- END USER VAR --------------
%------------- BEGIN CODE --------------

%% parse input options

if ~isfield(in_struct,'pzt_dist_sm') || isnan(in_struct.pzt_dist_sm) || isinf(in_struct.pzt_dist_sm) || isempty(in_struct.pzt_dist_sm)
    warning('you have not specified pzt_dist_sm which sets the mimumum distance (in pzt voltage) between peaks for laser to pass this single mode test')
    in_struct.pzt_dist_sm=inf;
    warning('this code will retun the distribution of peak distances so you can set one')
    
end
if ~isfield(in_struct,'pd_amp_min')
    in_struct.pd_amp_min=1;
end

if ~isfield(in_struct,'pd_filt_factor')
    warning('setting pd_filt_factor to default')
    in_struct.pd_filt_factor=1e-3;
end
if ~isfield(in_struct,'ptz_filt_factor_pks')
    warning('setting ptz_filt_factor_pks to default')
    in_struct.ptz_filt_factor_pks=1e-3;
end
if ~isfield(in_struct,'pzt_filt_factor_deriv')
    warning('setting pzt_filt_factor_deriv to default')
    in_struct.pzt_filt_factor_deriv=1e-3;
end

if ~isfield(in_struct,'pzt_dist_pd_cmp')
    warning('setting pzt_dist_pd_cmp to default')
    in_struct.pzt_dist_pd_cmp=0.1*in_struct.pzt_dist_sm;
end


if ~isfield(in_struct,'peak_thresh') || isempty(in_struct.peak_thresh) || sum(isnan(in_struct.peak_thresh))~=0
    warning('you did not specify a peak threshould we will try and guess what is should be based on the distribution of values in the pd input')
     in_struct.peak_thresh=[];
    for ii=1:size(in_struct.pd_voltage,2)
        [f,x]=ecdf(in_struct.pd_voltage(:,ii));
        cum_fun=@(xq) interp1(x(2:end),f(2:end),xq)-0.96;
        root=fzero(cum_fun,[min(x),max(x)]); %median(in_struct.pd_voltage(:,ii))
        in_struct.peak_thresh(ii) =root;
    end
end

if ~isfield(in_struct,'verbose') || isempty(in_struct.verbose)
    in_struct.verbose=1;
end

if ~isfield(in_struct,'plot') || isempty(in_struct.plot)
    in_struct.plot.all=false;
    in_struct.plot.failed=false;
else 
    if ~isfield(in_struct.plot,'all')
        in_struct.plot.all=false;
    end
    if ~isfield(in_struct.plot,'failed')
         in_struct.plot.failed=false;
    end
end

    
    
    

valid_pzt_scan_types={'sawtooth','triangle','sine'};
if ~ismember( in_struct.scan_type,valid_pzt_scan_types)
    warning('valid .scan_types are %s',strjoin(valid_pzt_scan_types,', '))
elseif isequal(in_struct.scan_type,'sine')
    warning('not yet supported')
end

%% measure the scan freq
% plot(in_struct.times,in_struct.pzt_voltage)
% if the pzt_scan_period has been specified take it as a guess and then we can take the precise answer using
% a fft of the first 10 oscillations
if isfield(in_struct,'pzt_scan_period')
    tmax=min(in_struct.times(1)+in_struct.pzt_scan_period*10,in_struct.times(end));
    [~,idx_max]=closest_value(in_struct.times,tmax);
    pztfreq_subset_data=[in_struct.times(1:idx_max),in_struct.pzt_voltage(1:idx_max)];
    %clear('tmax','idx_max')
else %otherwise we do the fft on the whole thing
    pztfreq_subset_data=[in_struct.times,in_struct.pzt_voltage];
end
% do a padded FFt to get some better freq resolution when just using a simple max fft bin
%subtract the mean to remove the DC peak (so the max bin operation does not get confused)
fft_out=fft_tx(pztfreq_subset_data(:,1),pztfreq_subset_data(:,2)-mean(pztfreq_subset_data(:,2)),'window','hanning','padding',10);

pzt_scan_std=std(pztfreq_subset_data(:,2));

%mask out the lowest and highest frequencies
idxs_freq_range=fast_sorted_mask(fft_out(1,:),pzt_scan_freq_lims(1),pzt_scan_freq_lims(2));
fft_out=fft_out(:,idxs_freq_range(1):idxs_freq_range(2));
%figure(5)
%semilogy(fft_out(1,:),abs(fft_out(2,:)))

%get the frequency with the max power
[~,idx_max_power_freq]=max(abs(fft_out(2,:)));
pzt_scan_freq=fft_out(1,idx_max_power_freq);
pzt_scan_period=1/pzt_scan_freq;
laser_sm_details.pzt_scan_period=pzt_scan_period;

%% checking pzt waveform
% for the sawtooth waveform we need to know the waveform polarity (long +ve gradient, or long -ve gradient) in order to
% properly split up the subsets.
% it is also prudent  to check that the pzt waveform has only two zero derivative points per period

% we again select the first 10 pzt scan periods
tmax=min(in_struct.times(1)+pzt_scan_period*10,in_struct.times(end));
[~,idx_max]=closest_value(in_struct.times,tmax);
wfcheck_subset_data=[in_struct.times(1:idx_max),in_struct.pzt_voltage(1:idx_max)];

[st_pts,level_xing]=stpt_and_level_xing(wfcheck_subset_data(:,1),wfcheck_subset_data(:,2),...
                                            in_struct.pzt_filt_factor_deriv*pzt_scan_period,nan);
                                     
%check that the level xing alternates in gradient
if ~is_alternating_logical_vec(level_xing.positive_deriv)
    error('problem in pzt voltage data, mean crossings of pzt signal does not alternate derivative')
end
%and the same in curvature
if ~is_alternating_logical_vec(st_pts.positive_curvature)
    error('problem in pzt voltage data, stpts of pzt signal does not alternate curvature')
end
%check that when there is a stpt above the mean that this implies that there is a negative curvature
if ~isequal(st_pts.above_mean,~st_pts.positive_curvature)
     error('problem in pzt voltage data, the pzt value at a stationay point is above the mean but does not have negative curvature ')
end

%now lets check that there are 2 st pts for each mean xing
%if the number of elements is odd trim by one element to take an even number
iimax=numel(level_xing.idx);
iimax=iimax-mod(iimax,2);
for ii=1:iimax
    idx_start=ii;
    idx_end=ii+1;
    number_of_st_pts=numel(fast_sorted_mask(st_pts.idx,idx_start,idx_end));
    if number_of_st_pts~=2
        error('pzt waveform has %u st pts between mean crossings',number_of_st_pts)
    end
end

%find how the odd and even st.pt. spacing compate

stpt_time_diff=diff(st_pts.time);
stpt_val_diff=diff(wfcheck_subset_data(st_pts.idx,2));
stpt_even_odd_time_diff=[mean(stpt_time_diff(2:2:end)),mean(stpt_time_diff(1:2:end))];
stpt_even_odd_val_diff=[mean(stpt_val_diff(2:2:end)),mean(stpt_val_diff(1:2:end))];
stpt_diff_parity_ratio=max(stpt_even_odd_time_diff)/min(stpt_even_odd_time_diff);

if strcmp(in_struct.scan_type,'sawtooth') && stpt_diff_parity_ratio<10
    warning('based on the timing differences between st. points of the pzt waveform you dont have a sawtooth (or the xing detection didnt work)')
elseif ismember(in_struct.scan_type,{'triangle','sine'}) && stpt_diff_parity_ratio>5
    warning('based on the timing differences between st. points of the pzt waveform you dont have a triangle or sine waveform (or the xing detection didnt work)')
end

%calculate the type of sawtooth long positive gradient or long negative gradient
if strcmp(in_struct.scan_type,'sawtooth') 
    [~,idx]=max(stpt_even_odd_time_diff);
    sawtooth_polarity_positive=stpt_even_odd_val_diff(idx)>0; %is the sawtooth a positive type
end
  

%% sample or full check
% determine roughly how many pzt scans are in the data
% if there are may more than in_struct.num_checks then just sample a
% subset of the scans
%otherwise split up the scans and interogate all of them

time_range=in_struct.times(end)- in_struct.times(1);
scans_in_data=time_range/pzt_scan_period;

% if the scan type is triangle or sine then there are twice the number of
% scans if positive and negative polarity is taken
if ismember(in_struct.scan_type,{'triangle','sine'})
    scans_in_data=scans_in_data*2;
end

use_subset_sample=scans_in_data>in_struct.num_checks*scan_num_exta_facor;


%% split up the scans and feed it to the is_single_scan_sm

scans=[];
single_opts=[];
single_opts.peak_thresh=in_struct.peak_thresh;
single_opts.pd_filt_factor=in_struct.pd_filt_factor;
single_opts.pzt_scan_period=pzt_scan_period;
single_opts.pzt_filt_factor=in_struct.ptz_filt_factor_pks;
single_opts.merge_cmp_full_pk_dist=in_struct.pzt_dist_pd_cmp;
single_opts.pzt_dist_sm=in_struct.pzt_dist_sm;
single_opts.pd_amp_min=in_struct.pd_amp_min;
single_opts.plot=in_struct.plot;
single_opts.scan_clip_pzt=scan_clip_factor*pzt_scan_std;
single_opts.cmp_multiplier_disp=cmp_multiplier_disp;
pk_spacing_list=[];


laser_ok=true;
if use_subset_sample
    if isinf(in_struct.num_checks)
        error('runtime error infinite number of checks but use_subset_sample=true')
    end
    test_times=linspace(in_struct.times(1),in_struct.times(end),in_struct.num_checks+2); 
    test_times=test_times(2:end-1);
    iimax=numel(test_times);
    %probe a scan centered arround test_time(ii)
    subset=[];
    for ii=1:iimax
        time_cen=test_times(ii);
        % to guarantee one good full scan we need to select 2 (and a little bit) scan periods
        time_lims=time_cen+[-1,1]*pzt_scan_period*1.1; 
        subset_idxs=fast_sorted_mask(in_struct.times,time_lims(1),time_lims(2));
        subset.time=in_struct.times(subset_idxs(1):subset_idxs(2));
        subset.pzt=in_struct.pzt_voltage(subset_idxs(1):subset_idxs(2));
        subset.pd=in_struct.pd_voltage(subset_idxs(1):subset_idxs(2),:);
%         figure(3)
%         subplot(2,1,1)
%         plot(subset.time,subset.pzt)
%         subplot(2,1,2)
%         plot(subset.time,subset.pd(:,1))
        st_pts=stpt_and_level_xing(subset.time,subset.pzt,...
                                            in_struct.pzt_filt_factor_deriv*pzt_scan_period,nan);
        st_pts.value=subset.pzt(st_pts.idx);
        pzt_diff_between_st_pts=diff(st_pts.value);
        %select pairs that have a large enough difference
        scan_mask=pzt_diff_between_st_pts>pzt_scan_std*st_pt_min_pzt_factor;
        if strcmp(in_struct.scan_type,'sawtooth') 
            %select the pairs of stationary points that give a differenc that is the right sign(positive if sawtooth_polarity_positive==true)
            if sawtooth_polarity_positive
                scan_mask=scan_mask & pzt_diff_between_st_pts>0;
            else
                scan_mask=scan_mask & pzt_diff_between_st_pts<0;
            end
        end
        scan_stpt_idxs=find(scan_mask,1); %find the first scan that satisfies the conditions
        %make a matrix with start and stop indicies for all scans (for the long slope)
        seg_start_idx=st_pts.idx(scan_stpt_idxs);
        seg_stop_idx=st_pts.idx(scan_stpt_idxs+1);

        
        single_opts.pd_voltage=subset.pd(seg_start_idx:seg_stop_idx,:);
        single_opts.time=subset.time(seg_start_idx:seg_stop_idx);
        single_opts.pzt_voltage=subset.pzt(seg_start_idx:seg_stop_idx);
        scans{ii}=is_single_scan_sm(single_opts);
        scans{ii}.start_idx=seg_start_idx;
        scans{ii}.stop_idx=seg_stop_idx;
        if ~scans{1}.is_power_ok || ~scans{1}.is_single_mode
            laser_ok=false;
        end
        
        
    end
    
else
    %get the stationary points of the pzt waveform
    st_pts=stpt_and_level_xing(in_struct.times,in_struct.pzt_voltage,...
                                            in_struct.pzt_filt_factor_deriv*pzt_scan_period,nan);
                                        
    st_pts.value=in_struct.pzt_voltage(st_pts.idx);
    pzt_diff_between_st_pts=diff(st_pts.value);
    %select pairs that have a large enough difference
    scan_mask=pzt_diff_between_st_pts>pzt_scan_std*st_pt_min_pzt_factor;
    if strcmp(in_struct.scan_type,'sawtooth') 
        %select the pairs of stationary points that give a differenc that is the right sign(positive if sawtooth_polarity_positive==true)
        if sawtooth_polarity_positive
            scan_mask=scan_mask & pzt_diff_between_st_pts>0;
        else
            scan_mask=scan_mask & pzt_diff_between_st_pts<0;
        end
    end
    scan_stpt_idxs=find(scan_mask); %find the first scan that satisfies the conditions
    %make a matrix with start and stop indicies for all scans (for the long slope)
    seg_start_idx=st_pts.idx(scan_stpt_idxs);
    seg_stop_idx=st_pts.idx(scan_stpt_idxs+1);                                   

    %trim extra components in the start index
    seg_start_idx=seg_start_idx(1:numel(seg_stop_idx)); 
    if ~isequal(size(seg_start_idx),size(seg_stop_idx)), error('runtime error start stop idx vec for segment not the same size'), end
    
    %some plots for debuging
%     figure(1)
%     subplot(2,1,1)
%     plot(in_struct.times(seg_start_idx(1):seg_stop_idx(1)),...
%         in_struct.pzt_voltage(seg_start_idx(1):seg_stop_idx(1)))
%     hold on
%     plot(in_struct.times(seg_start_idx(1):seg_stop_idx(1)),...
%         in_struct.pd_voltage(seg_start_idx(1):seg_stop_idx(1),2))
%     hold off
%     subplot(2,1,2)
%      plot(in_struct.times(seg_start_idx(end):seg_stop_idx(end)),...
%         in_struct.pzt_voltage(seg_start_idx(end):seg_stop_idx(end)))
%     hold on
%     plot(in_struct.times(seg_start_idx(end):seg_stop_idx(end)),...
%         in_struct.pd_voltage(seg_start_idx(end):seg_stop_idx(end),2))
%     hold off
    

    iimax=numel(seg_start_idx);
    for ii=1:iimax
        single_opts.pd_voltage=in_struct.pd_voltage(seg_start_idx(ii):seg_stop_idx(ii),:);
        single_opts.time=in_struct.times(seg_start_idx(ii):seg_stop_idx(ii));
        single_opts.pzt_voltage=in_struct.pzt_voltage(seg_start_idx(ii):seg_stop_idx(ii));
        scans{ii}=is_single_scan_sm(single_opts);
        scans{ii}.start_idx=seg_start_idx(ii);
        scans{ii}.stop_idx=seg_stop_idx(ii);
        if ~scans{1}.is_power_ok || ~scans{1}.is_single_mode
            laser_ok=false;
        end
    end
end %use subset of scans

laser_sm_details.scans=scans;

%%

tmp_cell=arrayfun(@(x) diff(x{1}.pks.all.pzt),scans,'UniformOutput',0);
pk_spacing_list=cat(1,tmp_cell{:});


tmp_cell=arrayfun(@(x) (x{1}.pks.full.pd),scans,'UniformOutput',0);
pk_amp_list=cat(1,tmp_cell{:});

%%

%TODO give more stats stats

if in_struct.verbose>1 || isinf(in_struct.pzt_dist_sm)
    fprintf('photodiode peak spacing in pzt voltage\n')
    fprintf('mean %.3f , std %.3f , median %.3f',mean(pk_spacing_list),std(pk_spacing_list),median(pk_spacing_list))
    figure('Name','is_laser_single_mode: pk spacing dist','NumberTitle','off')
    histogram(pk_spacing_list,round(numel(pk_spacing_list)/2))
end
if isinf(in_struct.pzt_dist_sm)
    laser_ok=false;
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

