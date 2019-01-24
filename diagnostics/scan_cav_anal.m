clear all
close all

% Let's not worry about the peak integration issues for now. Complete the
% workflow;
% First, import the null data. Select the forward scan and find the mean.
    % This is the zero offset (but one still finds negative values in the
    % scans, yes? At least they integrate to positive results. 
% For a given directory (filter configuration);
%   For each scan around a peak, find the average peak/bg ratio.
%   Return a list of all of these for each file.
% Output distribution of these ratios.
% Alright, that's done.
% Tomorrow I'll take some data with slower scans and hopefully that
% resolves the peaks better.
% For now, attention is required elsewhere.

% Methods for peak fixing:
    % Lebesgue method (sort, cut, integrate)
    % Height fixing (by max); width estimation
    % fit/local smoothing

% TODO: Correct for PD offset
% TODO: Check zero-peak files
% TODO: Using calibration offset, compute mean power ratio per file
% Then compile and interpret (i.e. physics)

    %ai_dat format:
    %   1: probe beam PD; low, because attenuated by SFP?
    % sfp_pzt_raw=ai_dat.Data(3,:);
    % sfp_pd_raw=ai_dat.Data(2,:);
    % compressed =ai_dat.Data(4,:);

fwtext('CAVITY SCAN ANALYSIS')

config.log_name='log_analog_in_*';
config.datadir = 'C:\Data\blue_sfp';
config.savedir = 'C:\Users\jacob\Documents\Projects\Tune_out_v2\figs\cavity_analysis';
config.scan_time=14e-3;  %estimate of the sfp scan time,used to set the window and the smoothin
config.pzt_volt_smothing_time=config.scan_time/100;
config.pzt_division=0.2498;
config.sampl_start = 1;
config.demo_max = 4e3;
config.treshold = 1e-1;
config.plot_out = false;
config.test.num_files = NaN;
config.test.num_scans = 10;

fwtext('DARK SCANS')
light_off_dirnames = {'probe_off','red_off'};
all_null = cell(numel(light_off_dirnames),1);
config_null = config;
for ii=1:numel(light_off_dirnames)
    config_null.dir = fullfile(config.datadir,light_off_dirnames{ii});    
    null_data = null_analysis(config_null);
    all_null{ii} = null_data.data;
end
fwtext('Null scan complete')
all_null = cell_vertcat(all_null);

%%

fwtext('LIGHT ON')
light_on_dirnames = {'probe_on_no_filt','probe_on_one_filt','probe_on_two_filt','probe_on_all_filt'};
config_light = config;
config_light.plot_out = false;
config_light.lorentz_fit = false;
config_light.zero_offset =mean(all_null{1});
config_light.peak_width = 0.5e-4;
config_light.logscale = false;

for ii=1:numel(light_on_dirnames)
    config_light.dirname= light_on_dirnames{ii};
    config_light.dir = fullfile(config.datadir,config_light.dirname);
    light_on_data = light_on_analysis(config_light);
    inspect_light_on(light_on_data,config_light)
    fprintf(['Completed ',config_light.dirname,', peak/back ratio %.4f +/- %.4f\n'],light_on_data.stats)
end
fwtext('All dirs done')



