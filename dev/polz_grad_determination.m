
clear all
hebec_constants();
cli_header('Getting polz data');
root_data_dir='E:\Data\Tuneout\to_main_data';
files = dir(root_data_dir);
files=files(3:end);
% Get a logical vector that tells which is a directory.
dir_mask = [files.isdir];
folders=files(dir_mask);
folders=arrayfun(@(x) fullfile(root_data_dir,x.name),folders,'UniformOutput' ,false);
loop_config.dir=folders;


% % Retrieve the PD data
% profile on
cache_opts.verbose = 0;
cache_opts.force_recalc = 0;
% root_data_dir='E:\Data\Tuneout\to_main_data';
% files = dir(root_data_dir);
% files=files(3:end);
% Get a logical vector that tells which is a directory.
% dir_mask = [files.isdir];
% folders=files(dir_mask);
% %convert to the full path
% folders=arrayfun(@(x) fullfile(root_data_dir,x.name),folders,'UniformOutput' ,false);
% %folders=folders(randperm(numel(folders)));%randomize the ordering
% folders;


% get pol data
scandata = load_pocessed_to_data(loop_config);

pol_opts.location = 'post';%post, pre_cen, pre_left, pre_right
pol_opts.predict_method = 'only_data';%'full_fit_pref_fit','full_fit_pref_data','full_fit_only','only_data'; %obs (obsovation) fit (pertial fit) full_fit (fit with all parameters free)
                                        %'interp_only','interp_pref_data','interp_pref_interp'
                                        %'gauss_only','gauss_pref_data','gauss_pref_interp'
pol_opts.smoothing=3; %deg
pol_opts.wrap_hwp=0;

pol_opts.hwp=scandata.drift.wp.hwp;
pol_opts.qwp=scandata.drift.wp.qwp;
pol_model=pol_data_query(pol_opts);

polz_theta=pol_model.theta.val;
polz_v=pol_model.v.val;
polz_cont=pol_model.cont.val;

% get PD etc
max_file_num = length(folders);
run_colours = viridis(max_file_num);
for folder_idx = 1:max_file_num
    cli_header(1,'Trying folder %u/%u',folder_idx,max_file_num);
    this_folder = folders{folder_idx};
    % is there the output we want?
    if ~exist(fullfile(this_folder,'out'),'dir')
       cli_header(2,'Output folder not detected, skipping analysis and plot');
    else
        try 
            dir_files = dir(fullfile(this_folder,'out'));
            dir_files =arrayfun(@(x) fullfile(x.name),dir_files,'UniformOutput' ,false);
            dir_files = dir_files(3:end);
            dfile_use = dir_files{end}; % hopefully
            % need: 
            load_folder = fullfile(this_folder,'out',dfile_use);
    
            cache_opts.dir = fullfile(load_folder);
            alpha_info{folder_idx} = simple_function_cache(cache_opts,@get_alpha_info,{load_folder});
            
            
            % recalculate the shotwise d_alpha_d_f error as it's the wrong
            % dimensionality
            % % We have the slope (omega_probe^2) in terms of some constants:
            % % omega_probe.^2 = 4*A*laser_power*real_part_of_alpha./beam_waist.^4;
            % % where
            A = (const.mhe*const.c*const.epsilon0)^-1;
            % % first attempt will just be guessing the power and using the relative
            % % stability of the PD
            PD_rel_unc = alpha_info{folder_idx}.scan_pd_stds./alpha_info{folder_idx}.scan_pd_means;
            beam_waist.val = 10e-6;
            beam_waist.unc = 0.05*beam_waist.val; % generous error margin
            % % hence
            % % real_part_of_alpha = (omega_probe.^2.*beam_waist.^4/(4*A*beam_power)
            % % and thus we have the gradient
            %     d_alpha_d_f.unc = (4*A)*sqrt( (PD_rel_unc.^4).^2 + ...
            %                                    (4*beam_power.val*beam_waist.unc*beam_waist.val^-5)^2).*data.to_fit_seg.fit_trimmed.slope.val;
            %         d_alpha_d_f.unc = PD_rel_unc.*data.to_fit_seg.fit_trimmed.d_alpha_d_f.val;
            %         d_alpha_d_f.val = fit_trimmed.d_alpha_d_f.val;

            % If this is done, then save; or wrap in cache function
            catch e
                fprintf(2,'caught error:\n%s\n',getReport(e,'extended'))
                stack_trace_report=arrayfun(@(x) sprintf('line %u: %s \n(%s)\n',x.line,x.name,x.file),e.stack,'UniformOutput',false);
                fprintf(2,'stack:\n%s',[stack_trace_report{:}])
                diary off
                %error('stop')
        end
    end
end

cli_header('plot done');
% title(sprintf('Results from 
% do next: Extract zero crossing and polarization from folders also
% Fit grad vs PD power?
alpha_info = alpha_info(~cellfun(@(x) isempty(x), alpha_info));
% profile off
% profile viewer

%% make a table out of the results

waist_est = [15,2]*1e-6;

same_size_filter = cell2mat(cellfun(@(x) isequal([0,0],size(x.scan_pd_means)-size(x.slope.val)), alpha_info,'uni',0));
pd_mean_cat = cell_vertcat(cellfun(@(x) x.scan_pd_means,alpha_info(same_size_filter),'uni',0)');
pd_unc_cat = cell_vertcat(cellfun(@(x) x.scan_pd_stds,alpha_info(same_size_filter),'uni',0)');
fit_slope_cat = cell_vertcat(cellfun(@(x) x.slope.val,alpha_info(same_size_filter),'uni',0)');%Hz^2 / Hz;
fit_slope_unc_cat = cell_vertcat(cellfun(@(x) x.slope.unc,alpha_info(same_size_filter),'uni',0)'); %Hz^2 / Hz;
to_time_cat = cell_vertcat(cellfun(@(x) x.to_time,alpha_info(same_size_filter),'uni',0)');
to_val_cat = cell_vertcat(cellfun(@(x) x.to.val,alpha_info(same_size_filter),'uni',0)');
to_val_unc_cat = cell_vertcat(cellfun(@(x) x.to.unc,alpha_info(same_size_filter),'uni',0)');

% Need a power conversion
% Taken from spetroscopy work, v similar wavelengths, hoping the PD
% settings are the same... But not looking likely so see below
        % set_pow_pairs = [0.35,0.0144;
        %                 0.15,0.00519;
        %                 0.1,0.0034;
        %                 0.10,3.31e-3;			
        %                 0.15,5.14e-3;			
        %                 0.20,6.5e-3; 			
        %                 0.25,8.9e-3;
        %                 0.011,0;
        %                 0.25, 0.00845];
        % set_pow_pairs = fliplr(set_pow_pairs);
        % pd_mdl = fitlm(set_pow_pairs(:,1),set_pow_pairs(:,2));
        % pow_conversion = pd_mdl.Coefficients.Estimate(2); %W / V

% from TO data we have more limited picture
% 20190205T1416_to_hwp_310_polmin_121_nuller_reconfig
% 						power previously 155mw at 1.25v
% 						after adj 138mw at 1.25v
% 20190208T2202_to_hwp_299_polmin_99.5_nuller_reconfig
% 						set 1.00 130mw after power samp bc,  86 at cammera
% 						1.25v = 134mw after power samp bc, 86mw at cammera
zero_val = 0*0.011;
pd_cal.dir = {'20190205T1416_to_hwp_310_polmin_121_nuller_reconfig','20190208T2202_to_hwp_299_polmin_99.5_nuller_reconfig'};
pd_cal.val = [1.25,1.25];
pd_cal.pow = [0.138,0.134];
pd_cal.factor = pd_cal.pow./(pd_cal.val);
pd_cal.run_vals = [nan,nan];
for ii=1:length(pd_cal.dir)
    pd_cal.run_vals(ii) = find(cellfun(@(x) contains(x,pd_cal.dir{ii}),folders));
end

%make run idx bit

num_runs = max(size(alpha_info(same_size_filter)));
num_scans_all = cellfun(@(x) size(x.scan_pd_means,1),alpha_info(same_size_filter));
run_labels = (1:sum(num_scans_all))';
cc = 0;
for ii = 1:num_runs
    for jj = 1:num_scans_all(ii)
       cc = cc + 1;
       run_labels(cc) = ii;
    end
end

data_table = table(to_time_cat{1},pd_mean_cat{1},pd_unc_cat{1},fit_slope_cat{1},fit_slope_unc_cat{1},...
    to_val_cat{1},to_val_unc_cat{1},pd_mean_cat{1}*mean(pd_cal.factor),run_labels,...
    'VariableNames',{'to_time','pd_mean','pd_unc','slope','slope_unc','to_val','to_unc','pd_power','run_idx'});
% size(data_table)

pol_table = table(scandata.drift.to_time,scandata.drift.wp.hwp,scandata.drift.wp.qwp,...
    pol_model.v.val,pol_model.v.unc,pol_model.theta.val,pol_model.theta.unc,...
    'VariableNames',...
    {'to_time','hwp','qwp','v','v_unc','theta','theta_unc'});

% %
% the data table is prob smaller
dtab_mask = nan(size(data_table,1),1);
for ii=1:size(data_table,1)
   dtab_mask(ii) = find(pol_table.to_time == data_table.to_time(ii)); 
end
pol_tab_take = pol_table(dtab_mask,:);

join_table = [data_table,pol_tab_take(:,2:end)];

join_table.conversion_factor = (8*join_table.pd_power)/(8*pi^3*const.mhe*const.c*const.epsilon0*waist_est(1)^4);
join_table.polz_grad =  join_table.slope ./ join_table.conversion_factor; 



% % Analyse and modelling
% showing that only the PD voltage explains the data
predictor_vars = table(join_table.pd_mean,join_table.v,join_table.theta);
wp_nan_mask = isnan(join_table.v) | isnan(join_table.theta);
wp_run_tab = join_table(~wp_nan_mask,:);

run_pd_means = nan(num_runs,1);
run_pd_uncs = nan(num_runs,1);
run_slope_means = nan(num_runs,1);
run_slope_uncs = nan(num_runs,1);
run_mid_time = nan(num_runs,1);
run_v = nan(num_runs,1);
run_theta = nan(num_runs,1);
run_to = nan(num_runs,1);
run_power = nan(num_runs,1);
run_grad = nan(num_runs,1);

for ii = 1:num_runs
    sub_table = join_table(join_table.run_idx==ii,:);
    run_pd_means(ii) = mean(sub_table.pd_mean);
    run_pd_uncs(ii) = sqrt(mean(sub_table.pd_unc.^2));
    run_slope_means(ii) = mean(sub_table.slope);
    run_slope_uncs(ii) = sqrt(mean(sub_table.slope_unc.^2));
    run_mid_time(ii) = mean(sub_table.to_time);
    run_v(ii) = mean(sub_table.v);
    run_theta(ii) = mean(sub_table.theta);
    run_to(ii) = mean(sub_table.to_val);
    run_power(ii) = mean(sub_table.pd_power);
    run_grad(ii) = mean(sub_table.polz_grad);
    
end

run_table = table(run_mid_time,run_pd_means,run_pd_uncs,run_slope_means,...
    run_slope_uncs,run_v,run_theta,run_to,run_power,run_grad);
run_mask = isnan(run_table.run_v) | isnan(run_table.run_theta);
pow_cal_table = join_table(any(join_table.run_idx == pd_cal.run_vals,2),:);

predictor_table = table(wp_run_tab.v,wp_run_tab.theta,wp_run_tab.pd_mean,'VariableNames',...
    {'V','Theta','PD'});
predictor_vals = (table2array(predictor_table));

[pd_sole_fit,pd_sole_gof,pd_sole_out] = fit(join_table.pd_mean,join_table.slope,'poly1');

linmdl = fitlm(predictor_vals,(wp_run_tab.slope),'PredictorVars',...
    {'V','Theta','PD'})
res_tab = anova(linmdl)


% now try using all the runs (noting we've filtred out those without WP
% data)
trunc_table = join_table(~isnan(join_table.qwp),:);




% linmdl = fitlm(

pow_mdl_cal = fitlm(pow_cal_table.pd_power,pow_cal_table.slope);
pow_mdl_wp = fitlm(wp_run_tab.pd_power,wp_run_tab.slope);
pow_mdl_full = fitlm(join_table.pd_power,join_table.slope);



conversion_factor = (8*pow_cal_table.pd_power)/(8*pi^3*const.mhe*const.c*const.epsilon0*waist_est(1)^4);
% slope = d Omega^2 / d f
%d Omega^2 /d f =  conv_factor * d alpha / d f thus
polz_grad =  pow_cal_table.slope ./ conversion_factor; % the slope is in GHz so need to convert to Hz
% compare to theory value


% Well that's not so awesome - factor 100 disagreement should be traceable
% to some error in the calculation... 

% using the full final data set
conversion_factor = (8*join_table.pd_power)/(8*pi^3*const.mhe*const.c*const.epsilon0*waist_est(1)^4);
% slope = d Omega^2 / d f
%d Omega^2 /d f =  conv_factor * d alpha / d f thus
polz_grad =  join_table.slope ./ conversion_factor; % the slope is in GHz so need to convert to Hz
% compare to theory value
trunc_table.conversion_factor = (8*trunc_table.pd_power)/(8*pi^3*const.mhe*const.c*const.epsilon0*waist_est(1)^4);
trunc_table.polz_grad = trunc_table.slope ./ trunc_table.conversion_factor;



dAlpha_dOmega_theory = const.a0^3 * .108e-5; % per MHz
cli_header('Results:');
cli_header(1,'Mean conversion factor: %.2e',mean(conversion_factor));
cli_header(1,'Avg polz grad %.3e(%.3e)',mean(trunc_table.polz_grad),std(trunc_table.polz_grad));
cli_header(1,'Vs theory value %.3e',dAlpha_dOmega_theory*1e-6); %convert to per Hz
cli_header(1,'Using full dataset:');
cli_header(2,'Avg polz grad %.3e(%.3e)',mean(polz_grad),std(polz_grad));
% Well that's not so awesome - factor 100 disagreement should be traceable
% to some error in the calculation... 
% Could the slope be underestimated? Or is the theory value wrong? 
% from the plot inset we ahve 1e-41 / THz = tiny! 1e-55/Hz? way smaller
% than theory seems to indicate?
% using the full final data set
%%
P = [mean(trunc_table.pd_power),std(trunc_table.pd_power)];
w = [15e-6,10e-6];
A = (P(1).*w(1).^-4)/(const.mhe*const.epsilon0*const.c*pi^3)
A_unc = sqrt((-4*(P(1).*w(1).^-5)/(const.mhe*const.epsilon0*const.c*pi^3)* w(2)).^2 + ...
         ((1.*w(1).^-4)/(const.mhe*const.epsilon0*const.c*pi^3)*P(2)).^2)
Om_grad = 30e-9;
polzgrad=Om_grad / A
polzgradunc=Om_grad * A_unc / A^2
thry = 1.78e-53;
[polzgrad,polzgradunc]/thry
%%

% run_pow_mdl = fitlm(run_predictor_table(:,3),

% % Do plots

h = 2;
w = 3;
stfig('out 2');
clf

subplot(h,w,1)

yyaxis left
hold on
% errorbar(join_table.to_time,join_table.slope,join_table.slope_unc,...
%     'r.','MarkerSize',1)
errorbar(run_table.run_mid_time-min(run_table.run_mid_time),run_table.run_slope_means,...
    run_table.run_slope_uncs,run_table.run_slope_uncs,...
    0*run_table.run_pd_uncs,0*run_table.run_pd_uncs,...
    '.','MarkerSize',15)

ylabel('Fit slope')
yyaxis right

hold on
% errorbar(join_table.to_time,join_table.pd_mean,join_table.pd_unc,...
%     'k.','MarkerSize',1)
errorbar(run_table.run_mid_time-min(run_table.run_mid_time),run_table.run_pd_means,run_table.run_pd_uncs,...
    '.','MarkerSize',5)
ylabel('PD voltage')
xlabel('scan timestamp')



set(gca,'FontSize',14)

subplot(h,w,2)
cla
hold on

% errorbar(run_table.run_slope_means,run_table.run_theta,...
%     0*run_table.run_pd_uncs,0*run_table.run_pd_uncs,...
%         run_table.run_slope_uncs,run_table.run_slope_uncs,...
%     'k.','MarkerSize',15)
% errorbar(run_table.run_slope_means,run_table.run_v,...
%     0*run_table.run_pd_uncs,0*run_table.run_pd_uncs,...
%     run_table.run_slope_uncs,run_table.run_slope_uncs,...
%     'r.','MarkerSize',15)
errorbar(run_table.run_theta,run_table.run_slope_means,...
        run_table.run_slope_uncs,run_table.run_slope_uncs,...
        0*run_table.run_pd_uncs,0*run_table.run_pd_uncs,...
    'k.','MarkerSize',15)
errorbar(run_table.run_v,run_table.run_slope_means,...
    run_table.run_slope_uncs,run_table.run_slope_uncs,...
        0*run_table.run_pd_uncs,0*run_table.run_pd_uncs,...
    'r.','MarkerSize',15)
ylabel('Fit slope')
xlabel('polz parameter')
legend('theta','v')
set(gca,'FontSize',14)



subplot(h,w,3)
hold on
errorbar(join_table.pd_mean,join_table.slope,...
        join_table.slope_unc,join_table.slope_unc,...
            join_table.pd_unc,join_table.pd_unc,...
    'k.','MarkerSize',5)

errorbar(wp_run_tab.pd_mean,wp_run_tab.slope,...
        wp_run_tab.slope_unc,wp_run_tab.slope_unc,...
            wp_run_tab.pd_unc,wp_run_tab.pd_unc,...
    'b.','MarkerSize',5)

errorbar(run_table.run_pd_means,run_table.run_slope_means,...
    run_table.run_slope_uncs,run_table.run_slope_uncs,...
    run_table.run_pd_uncs,run_table.run_pd_uncs,...
    'g.','MarkerSize',5)

errorbar(pow_cal_table.pd_mean,pow_cal_table.slope,...
        pow_cal_table.slope_unc,pow_cal_table.slope_unc,...
            pow_cal_table.pd_unc,pow_cal_table.pd_unc,...
    'r.','MarkerSize',5)

legend({'Scans','with wp','run avgs','With pow'},'Location','SouthEast')
xlabel('PD voltage (V)')
ylabel('Fit slope (Hz$^2$/Hz)')
set(gca,'FontSize',14)



subplot(h,w,4)
hold on
plot(join_table.pd_power,join_table.slope,...
    'k.','MarkerSize',10)
plot(wp_run_tab.pd_power,wp_run_tab.slope,'g.','MarkerSize',10)
plot(pow_cal_table.pd_power,pow_cal_table.slope,'r.','MarkerSize',10)
ylabel('Fit slope (Hz$^2$/Hz)')
xlabel('Pow estimate (W)')
legend({'All Scans','scans with wp','scans With pow'},'Location','SouthEast')
set(gca,'FontSize',14)

subplot(2,3,5)
hold on
plot(join_table.pd_power,polz_grad,'k.','MarkerSize',10)
plot(trunc_table.pd_power,trunc_table.polz_grad,'g.','MarkerSize',10)
plot(run_table.run_power,run_table.run_grad,'bx','MarkerSize',10)
legend('All scans','scans with WP data','run averages')
xlabel('Probe power')
ylabel('d$\alpha$/d$\omega$')
set(gca,'FontSize',14)

hedges = 1.1e54*linspace(min(join_table.polz_grad),max(join_table.polz_grad),41);
subplot(2,3,6)
hold on
histogram(1e54*join_table.polz_grad,hedges,'Normalization','PDF')
histogram(1e54*trunc_table.polz_grad,hedges,'Normalization','PDF')
legend('All scans','WP run avg')
xlabel('d $\alpha$ /d $\omega$ run avg ($\times10^{-54}$)')
ylabel('Number of runs')
set(gca,'FontSize',14)


%%
cli_header('Theory values');
d_alpha_df_theory_CGS = 0.108e-5 * const.a0^3 %CGS, converting m/MHz to cm/Hz
d_alpha_df_theory_SI = d_alpha_df_theory_CGS*(4*pi*const.epsilon0)/1e6;
alpha_T_cgs = 2*0.1841e-2 * const.a0^3; %in CGS bohr-radii
alpha_T_SI  = 1e-6*alpha_T_cgs*(4*pi*const.epsilon0)
beta_T_pred = alpha_T_SI/d_alpha_df_theory_SI

bT= [5430.1 -1911.0 4467.2
    2635.0 -927.3 2167.7 ];
bT(:,2:end) = bT(:,2:end) + bT(:,1);
mean_bT = mean(bT(:,1))
range_bT = [min(bT,[],'all'),max(bT,[],'all')]
%%

vals.dalpha.exp = [4,2]*1e-54;
vals.dalpha.theory = 1.78e-53;
vals.betaT.exp = [6,2.5]*1e3;
vals.betaT.theory = 3.4e3;
vals.betaV.exp = [1.5,.1]*1e4;
vals.betaV.theory = nan;

vals.ratio.dalpha = vals.dalpha.exp./vals.dalpha.theory;
vals.ratio.betaT= vals.betaT.exp./vals.betaT.theory;
vals.ratio.betaV = vals.betaV.exp./vals.betaV.theory;
vals.ratio



%%

function alpha_info = get_alpha_info(target_dir)

%          run_alpha_info{folder_idx} = simple_function_cache(get_alpha_info{folder_idx},target_path,cache_opts);
        
%         cli_header(2,'Output folder last entry:');
%         dir_files = dir(fullfile(this_folder,'out'));
%         dir_files =arrayfun(@(x) fullfile(x.name),dir_files,'UniformOutput' ,false);
%         dir_files = dir_files(3:end);
%         dfile_use = dir_files{end}; % hopefully
%         % need: 
%         cli_header(1,'loading data_anal_full from %s',this_folder);
%         cli_header(2,'Using output %s...',dfile_use);    
        if ~exist(fullfile(target_dir,'data_anal_full.mat'))
            cli_header('No data file found in out folder, skipping');
            alpha_info = [];
        else
            
            loadfile = fullfile(target_dir,'data_anal_full');
            data_full = load(loadfile);
            cli_header('Loading done');    



            %% Process and plot
            data = data_full.data;
            pd_mean_all = data.ai_log.pd.mean(data.to_fit_seg.fit_mask);
            pd_std_all = data.ai_log.pd.std(data.to_fit_seg.fit_mask);
            pd_median_all = data.ai_log.pd.median(data.to_fit_seg.fit_mask);
            % data.ai_log.ok.reg_pd - should be included in fit_mask?


            shots_used = data.mcp_tdc.shot_num(data.to_fit_seg.fit_mask);
            % segmented_shots = which_seg_are_you_in(shot_used);
            num_shots = length(data.mcp_tdc.shot_num);
            fprintf('Getting probe behaviour in scan segments');
            num_scans = length(data.to_fit_seg.scan_edges);
            iimax = num_scans-1;

            scan_pd_means = zeros(iimax,1);
            scan_pd_stds = zeros(iimax,1);
            scan_pd_medians = zeros(iimax,1);
            for ii=1:iimax 
                seg_start= data.to_fit_seg.scan_edges(ii)+1;
                seg_end = data.to_fit_seg.scan_edges(ii+1);
    %             cli_header(2,'scan %u Shots %u:%u', ii,seg_start,seg_end);
                seg_mask_temp=col_vec(1:num_shots);
                seg_mask_temp=(seg_mask_temp>=seg_start & seg_mask_temp<seg_end);
                seg_mask = seg_mask_temp & data.to_fit_seg.fit_mask;
    %             cli_header(2,'%u/%u shots kept in scan %u', sum(seg_mask),sum(seg_mask_temp),ii);
                scan_pd_means(ii) = mean(data.ai_log.pd.mean(seg_mask));
                scan_pd_medians(ii) = mean(data.ai_log.pd.median(seg_mask));
                scan_pd_stds(ii) = sqrt(sum(data.ai_log.pd.std(seg_mask).^2));
            end


                alpha_info.to.val = data.to_fit_seg.fit_trimmed.to_freq.val;
               alpha_info.to.unc = data.to_fit_seg.fit_trimmed.to_freq.unc;
               alpha_info.slope.val = data.to_fit_seg.fit_trimmed.slope.val;
               alpha_info.slope.unc = data.to_fit_seg.fit_trimmed.slope.unc;
               alpha_info.to_time = data.to_fit_seg.to_time;
    % Why are the above and below different sizes sometimes?
                alpha_info.pd_mean_all = pd_mean_all;
                alpha_info.pd_std_all = pd_std_all;
                alpha_info.pd_median_all = pd_median_all;
                alpha_info.scan_pd_means = scan_pd_means;
                alpha_info.scan_pd_medians = scan_pd_medians;
                alpha_info.scan_pd_stds = scan_pd_stds;
                
        end
end