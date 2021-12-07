%make_nice_sas_plot

data_dir='Z:\EXPERIMENT-DATA\2019_Forbidden_Transition\supplementary_data\wm_calibrations\20190715T0941_wm_cal';
%file_name='2p_cs_log_20181231T014329.txt'
files=dir(fullfile(data_dir,'2p_cs_log_*.txt'));
files={files.name};
scans_data=cell(1,0);
data_idx=1;
iimax=numel(files);
file_lengths=nan(1,iimax);
for ii=1:numel(files)
    file_name=files{ii}
    path=fullfile(data_dir,file_name);
    fid = fopen(path,'r');
    raw_lines = textscan(fid,'%s','Delimiter','\n'); % this is the fastest string read method 3.5s for 20 files
    fclose(fid);
    raw_lines=raw_lines{1};
    jjmax=numel(raw_lines)
    file_lengths(ii)=jjmax;
    for jj=1:jjmax
        tmp_st=jsondecode(raw_lines{jj});
        tmp_st.log_fname=file_name;
        scans_data{data_idx}=tmp_st;
        data_idx=data_idx+1;
    end
end

%%
f_table=[];
f_table.cs_6SF4_6PF2=351725718.50-4021.776399375-339.64; %MHZ
f_table.cs_6SF4_6PF3=351725718.50-4021.776399375-188.44; %MHZ
f_table.cs_6SF4_6PF4=351725718.50-4021.776399375+12.815; %MHZ
f_table.cs_6SF4_6PF5=351725718.50-4021.776399375+263.81; %MHZ
f_table.cs_6SF4_6PF3co4=(f_table.cs_6SF4_6PF3+f_table.cs_6SF4_6PF4)/2;
f_table.cs_6SF4_6PF4co5=(f_table.cs_6SF4_6PF4+f_table.cs_6SF4_6PF5)/2;

%% fit a single scan with a voigt
pmt_gain=20e6;
rescale_curr_mult=1e9;
pmt_gain=1;
rescale_curr_mult=1;
ii=1; %14   
outlier_thresh=3;

single_scan_dat=scans_data{ii};
rel_samp_time=single_scan_dat.parameters.sample_time_posix;
start_samp_time=rel_samp_time(1); %posix time is first entry
rel_samp_time(1)=0;
aq_cen_time=mean(rel_samp_time)+start_samp_time;
set_freq=single_scan_dat.parameters.set_freq;

num_samples_guess=500*(mean(diff(single_scan_dat.parameters.sample_time_posix(2:end)))-0.5);

pmt_curr_mean=single_scan_dat.parameters.pmt_voltage_mean/pmt_gain;
this_curr_std=single_scan_dat.parameters.pmt_voltage_std/pmt_gain;
this_curr_std=abs(this_curr_std);
pmt_curr_ste=this_curr_std/sqrt(num_samples_guess);
mean_freq=mean(set_freq);

%predictor=set_freq-mean_freq;
predictor=set_freq-f_table.cs_6SF4_6PF4co5;
response=pmt_curr_mean*rescale_curr_mult;
opts=[];
opts.predictor=predictor;
opts.response=response;
opts.response_err=pmt_curr_ste*rescale_curr_mult;
%opts.response_noise=this_curr_std*rescale_curr_mult;
opts.num_samples_per_pt=num_samples_guess;
opts.do_plots=true;
opts.sigma_outlier=outlier_thresh;
tmp_out_st=fit_1p_data(opts);
%tmp_out_st.fit_noise_no_cull
tmp_out_st.fit_cull.voigt