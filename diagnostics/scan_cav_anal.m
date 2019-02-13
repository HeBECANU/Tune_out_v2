% clear all
% close all

fwtext('CAVITY SCAN ANALYSIS')

config.log_name='log_analog_in_*';

% config.datadir = 'C:\Data\blue_sfp';
% light_on_dirnames = {'probe_on_one_filt'};

config.datadir = 'C:\Data\20190125_sfp';
light_on_dirnames = {'sfp_0_filter','sfp_1_filter','sfp_2_filter','sfp_3_filter'};
config.savedir = 'C:\Users\jacob\Documents\Projects\Tune_out_v2\figs\cavity_analysis';
config.scan_time=14e-3;  %estimate of the sfp scan time,used to set the window and the smoothin
config.pzt_volt_smothing_time=config.scan_time/100;
config.pzt_division=0.2498;
config.sampl_start = 1;
config.demo_max = 4e3;
config.treshold = 1e-1; %min peak height for find_peaks
config.plot_out = false;
config.R = 0.9875;
config.pv_method = true;
config.lebesgue_method = true;
config.num_dirs = 4;
config.test.num_files = 4;
config.test.num_scans = 100;


%% Calibration data
% fwtext('DARK SCANS')
% null_dirs = 'C:\Data\blue_sfp';
% light_off_dirnames = {'probe_off','red_off'};
% all_null = cell(numel(light_off_dirnames),1);
% config_null = config;
% for ii=1:numel(light_off_dirnames)
%     config_null.dir = fullfile(null_dirs,light_off_dirnames{ii});    
%     null_data = null_analysis(config_null);
%     all_null{ii} = null_data.data;
% end
% fwtext('Null scan complete')
% all_null = cell_vertcat(all_null);

%% Import light-on data
% 
fwtext('LIGHT ON')
config_light = config;
config_light.plot_out = false;
config_light.lorentz_fit = false;
config_light.zero_offset =-0.017;%mean(all_null{1});
config_light.peak_width = 1e-4;
config_light.window_size = 2e-3;
config_light.valley_width = 0.7;
config_light.logscale = false;
config_light.lebesgue_thresh = 1e-2*5; %Empirically chosen
config_light.lebesgue_thresh_demo = .5; %~1% of peak height
config_light.pv_peak_factor = 5;
config_light.pv_plot = false;
light_on_data = cell(numel(light_on_dirnames),1);
if ~isnan(config.num_dirs)
    num_dirs = config.num_dirs;
else
    num_dirs = numel(light_on_dirnames);
end
for ii=1:num_dirs
    fwtext({'Starting on %s',light_on_dirnames{ii}})
    config_light.dirname= light_on_dirnames{ii};
    config_light.dir = fullfile(config.datadir,config_light.dirname);
    light_on_data{ii} = light_on_analysis(config_light);
end
fwtext('All dirs analyzed')

%%
% close all
fwtext('Presenting results')
config_disp = config_light;
config_disp.lbg = false;
config_disp.insp_hists= true;


if config_disp.insp_hists
    sfigure(1337);
    clf;
end
insp_data = cell(num_dirs,1);
for ii=1:num_dirs

    insp_data{ii} = inspect_light_on(light_on_data{ii},config_disp);
    if config_disp.lbg
        present_lebesgue(insp_data{ii},config)        
    end
end
if config.pv_method
    sfigure(200);
    X = 0:num_dirs-1;
    Y = cellfun(@(x) x.pv_stats(1),insp_data);
    Y_err = cellfun(@(x) x.pv_stats(2),insp_data);
    Y_num = cellfun(@(x) x.pv_stats(3),insp_data);
    Y_stderr = Y_err./Y_num;
    errorbar(X,Y,Y_stderr,'bo')
    title('Background dependence on filters: peak/valley method')
    xlabel('Number of filters')
    ylabel('Background/peak ratio')
    xlim([-1,num_dirs])
    xticks([0,1,2,3])
end


if config.lebesgue_method
    fprintf('Lebesgue cutoff used: %f\n', config_light.lebesgue_thresh)

    sfigure(100)
    % subplot(1,2,1)
    X = 0:num_dirs-1;
    Y = cellfun(@(x) x.L_mean,insp_data);
    Y_err = cellfun(@(x) x.L_std,insp_data);
    Y_num = cellfun(@(x) x.L_num,insp_data);
    Y_stderr = Y_err./Y_num;
    errorbar(X,Y,Y_stderr,'ro')
    title('Background dependence on filters: Lebesgue method')
    xlabel('Number of filters')
    ylabel('Background/peak ratio')
    xlim([-1,num_dirs])
    xticks([0,1,2,3])
end

%% So yep, the worst place to put a delta-function perturbation is at 1sd.
% Why? If the filters are Gaussian, the shift is linear in detuning and
% power -> shift = power*detuning, power going like a Gaussian

% So: Given the ratio, put a Lorentzian of the same width with a given
%power at 1SD? 
% 

% subplot(1,2,2)
fwtext('Done!')
fwtext('')


% allstats = cell_vertcat(cellfun(@(x) x.stats, light_on_data,'UniformOutput',false)');
% allstats=allstats{1};
% errorbar(allstats(:,1),allstats(:,2))
% title('back/peak ratio vs num filters')



function present_lebesgue(insp_data,config)

%     figure()
%     for jj = 1:numel(light_on_data{ii}) %loop over files
%        for kk = 1:numel(light_on_data{ii}{jj}) %loop over scans
%             plot(sort(light_on_data{ii}{jj}{kk}.sweep))
%             hold on
%        end
%     end
%     set(gca,'Yscale','log')
%     title('Sorted scan data')
    
    LB = insp_data.L_back;
    fprintf('Directory ratio %.5f +- %.5f, integrated BG %e +- %e \n',...
        nanmean(insp_data.all_ratios),nanstd(insp_data.all_ratios),mean(LB),std(LB))
    
    if config.disp.lbg
        figure()
    %     suptitle(sprintf('Directory %u', ii))
        subplot(2,1,1)

        for jj=1:numel(insp_data.L_ratios)
            plot(insp_data.L_ratios{jj},'x')
            hold on
        end
        title('Back/peak ratios by scan')
        xlabel('Scan number')
        ylabel('Back/peak ratio: Lebesgue method')
        ylim([0,0.08])


        subplot(4,1,3)
        for jj = 1:numel(insp_data.varf_ratios)
            Y = insp_data.varf_ratios{jj};
            X = insp_data.varf_axis{jj};
            for kk = 1:numel(X)
                plot(X{kk},Y{kk},'k')
                hold on
            end
        end
        title(sprintf('Ratio vs cutoff for dir'))
        xlabel('Integration cutoff')
        ylabel('Background/peak ratio')
        ylim([-0.01,0.15])

        subplot(4,1,4)
        for jj = 1:numel(insp_data.varf_back)
            Y = insp_data.varf_back{jj};
            X = insp_data.varf_axis{jj};
            for kk = 1:numel(X)
                plot(X{kk},Y{kk},'k')
                hold on
            end
        end
        title('Back values')
    end
end