%% bryce scanning cav background measurement

% this is my rewrite of the scaningin cavity analysis
% a lot of code will be grabbed from jacobs code
% one of the primary reasons is that i will try to block and unblock the scanning cavity from the light in order to do a
% 'chopper' measurment of the background
% the code will have to automaticaly detect when the light is blocked and not
% also a lot of development has been done in the is_laser_single_mode project which should make all of this a lot easier
%

%%
fprintf('setting up path\n') %cli_format_text may not be in path yet

% find this .m file's path, this must be in the project root dir
project_root_folder = fileparts(fileparts(fileparts(which(mfilename))));
cd(project_root_folder)
% Add that folder plus all subfolders to the path.
addpath(genpath(project_root_folder));%add all subfolders to the path to find genpath_exclude
path_to_genpath=fileparts(which('genpath_exclude'));
path(pathdef) %clean up the path back to the default state to remove all the .git that were added
addpath(project_root_folder)
addpath(path_to_genpath)
addpath(genpath_exclude(fullfile(project_root_folder,'lib'),'\.')) %dont add hidden folders
addpath(genpath_exclude(fullfile(project_root_folder,'dev'),'\.'))
addpath(genpath_exclude(fullfile(project_root_folder,'bin'),'\.'))
addpath(genpath_exclude(fullfile(project_root_folder,'diagnostics'),'\.'))
hebec_constants %call the constants function that makes some globals


cli_format_text('Starting','cen',1)
%%
config=[]
config.dir='..\20190722_sfp_measurments\blue_4_filt_chop';
config.log_name='log_analog_in_';

%%
dir_read=dir(fullfile(config.dir,[config.log_name,'*.txt']));
file_names={dir_read.name};
args_single=[];
args_single.plot.all=false;
args_single.dir=config.dir;
%set the widths to integerate over
args_single.int_width.peak=0.1;
args_single.int_width.trough=0.5;
args_single.sfp_finesse=234;


anal_out=[]
for ii=1:2%numel(file_names)
    args_single.fname=file_names{ii};
    anal_out{ii}=background_single_file(args_single);
end




% find the mean backgroudn power
measured_background_ratios=cellfun(@(x) x.pd_light_on.background_pow,anal_out,'UniformOutput',false);
measured_background_ratios=cat(1,measured_background_ratios{:});

fprintf('mean background ratio %s\n',string_value_with_unc(mean(measured_background_ratios),std(measured_background_ratios)/numel(measured_background_ratios),'b'))

backgrounds_hist=smooth_hist(measured_background_ratios,'sigma',3e-5);
stfig('backgroudn measurment');
area( backgrounds_hist.bin.centers,backgrounds_hist.count_rate.smooth)

%%

