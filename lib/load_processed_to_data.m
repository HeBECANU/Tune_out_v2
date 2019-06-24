function data = load_pocessed_to_data(loop_config)
selected_dirs = 1:numel(loop_config.dir); %which files to loop over (currently all)

drift_data_compiled.to_val{:,1}=[];
drift_data_compiled.to_val{:,2}=[];
drift_data_compiled.to_time=[];
drift_data_compiled.atom_num=[];
drift_data_compiled.grad{:,1}=[];
drift_data_compiled.grad{:,2}=[];
drift_data_compiled.avg_coef = [];
drift_data_compiled.avg_coef_cal = [];
drift_data_compiled.avg_coef_unc = [];
drift_data_compiled.avg_coef_cal_unc = [];
main_data_compiled.lin_fit{1} = [];
main_data_compiled.lin_fit{2} = [];
main_data_compiled.quad_fit{1} = [];
main_data_compiled.quad_fit{2} = [];
main_data_compiled.shots = [];
main_data_compiled.time = [];
main_data_compiled.grad = [];
main_data_compiled.freq = [];
main_data_compiled.scan_num = [];
%break_idx = [];
main_data_compiled.set_pt = [];
for loop_idx=selected_dirs
    current_dir = loop_config.dir{loop_idx};
    if ~strcmp(current_dir(end),'\')
        current_dir = [current_dir,'\'];
    end
    if exist('Z:\EXPERIMENT-DATA\2018_Tune_Out_V2','file')~=7, error(sprintf('dir \n %s \n does not exist',current_dir)) ,end
    out_dirs=dir([current_dir,'out\']);
    if size(out_dirs,1)==0, error(sprintf('dir \n %s \n does not contain any out dirs',current_dir)), end 
    offset=0; %to make sure file contains the mat files
    most_recent_dir=out_dirs(end-offset,1);
    %runs through all the out put dirs for a given run and looks for saved data, if none is there
    %skips that data dir
    try
        while ~and(isfile([current_dir,'out\',most_recent_dir.name,'\main_data.mat']),isfile([current_dir,'out\',most_recent_dir.name,'\drift_data.mat']))
            offset = offset + 1;
            most_recent_dir=out_dirs(end-offset,1);
            %check = drift_data.avg_coef; %check if it has the avg coefs update
        end
    catch e
        fprintf('\n dir: %s didnt work \n',current_dir)
        msgText = getReport(e)
        continue
    end
    load([current_dir,'out\',most_recent_dir.name,'\main_data.mat'])
    load([current_dir,'out\',most_recent_dir.name,'\drift_data.mat'])
    
    %scrapes set points if you want to redo some analysis
    try
        %fprintf('\n dir: %s , st pt = %u \n',current_dir,main_data.set_pt)
        main_data_compiled.set_pt = [main_data_compiled.set_pt,main_data.set_pt];
    catch
        fprintf('\n dir: %s , no st pt record \n',current_dir)
    end
    
    %append to main structure
    drift_data_compiled.to_val{:,1}=[drift_data_compiled.to_val{:,1};drift_data.to_val{:,1}];
    drift_data_compiled.to_val{:,2}=[drift_data_compiled.to_val{:,2};drift_data.to_val{:,2};];
    drift_data_compiled.to_time=[drift_data_compiled.to_time;drift_data.to_time];
    drift_data_compiled.atom_num=[drift_data_compiled.atom_num;drift_data.atom_num];
    grad_temp=cell2mat(drift_data.grad);
    drift_data_compiled.grad{:,1}=[drift_data_compiled.grad{:,1};grad_temp(:,1)];
    drift_data_compiled.grad{:,2}=[drift_data_compiled.grad{:,2};grad_temp(:,2)];
    
    temp = cell2mat(drift_data.avg_coefs);
    c_uncs = reshape(temp(:,2),8,numel(drift_data.avg_coefs))';
    c_vals = reshape(temp(:,1),8,numel(drift_data.avg_coefs))';
    
    drift_data_compiled.avg_coef = [drift_data_compiled.avg_coef;c_vals];
    drift_data_compiled.avg_coef_unc = [drift_data_compiled.avg_coef_unc;c_uncs];
    
    try
    R = drift_data.avg_coefs_cal;
    R = R(~cellfun('isempty',R));
    temp = cell2mat(R);
    c_uncs = reshape(temp(:,2),8,numel(R))';
    c_vals = reshape(temp(:,1),8,numel(R))';
    
    drift_data_compiled.avg_coef_cal = [drift_data_compiled.avg_coef_cal;c_vals];
    drift_data_compiled.avg_coef_cal_unc = [drift_data_compiled.avg_coef_cal_unc;c_uncs];
    catch
    end
    
    main_data_compiled.lin_fit{1} = [main_data_compiled.lin_fit{1};main_data.lin_fit{1}];
    main_data_compiled.lin_fit{2} = [main_data_compiled.lin_fit{2};main_data.lin_fit{2}];
    main_data_compiled.quad_fit{1} = [main_data_compiled.quad_fit{1};main_data.quad_fit{1}];
    main_data_compiled.quad_fit{2} = [main_data_compiled.quad_fit{2};main_data.quad_fit{2}];
    main_data_compiled.shots = [main_data_compiled.shots;main_data.shots];
    main_data_compiled.time = [main_data_compiled.time;nanmean(drift_data.to_time)];
    main_data_compiled.scan_num = [main_data_compiled.scan_num;numel(drift_data.atom_num)];
    main_data_compiled.grad = [main_data_compiled.grad;nanmean(grad_temp(:,1))];
    main_data_compiled.freq = [main_data_compiled.freq;nanmean(c_vals(:,2))];
    
    %break_idx = [break_idx;numel(drift_data.to_time)];
end
data.drift = drift_data_compiled;
data.main = main_data_compiled;
end



% processing steps

% 
% clear drift_data main_data %To make sure there isn't any bleed through between directories
% 
% %Scan segmented data
% drift_data.to_val{:,1}=data.to_fit_seg.fit_trimmed.freq.val;
% drift_data.to_val{:,2}=data.to_fit_seg.fit_trimmed.freq.unc;
% drift_data.to_time=data.to_fit_seg.to_time;
% drift_data.atom_num=data.to_fit_seg.atom_num(:,1);
% for kk=1:numel(data.to_fit_seg.fit_trimmed.model)
%     drift_data.grad{kk,1}=data.to_fit_seg.fit_trimmed.model{kk,1}.Coefficients{2,1};
%     drift_data.grad{kk,2}=data.to_fit_seg.fit_trimmed.model{kk,1}.Coefficients{2,2};
% end
% drift_data.model=data.to_fit_seg.fit_trimmed.model;
% drift_data.avg_coefs = data.to_fit_seg.avg_coefs;
% drift_data.avg_coefs_cal = data.to_fit_seg.avg_coefs_cal;
% 
% %The analysis of the whole run
% main_data.set_pt = anal_opts.probe_set_pt;
% main_data.lin_fit{1} = to_freq_val_lin;
% main_data.lin_fit{2} = to_freq_unc_lin;
% main_data.quad_fit{1} = to_freq_val_quad;
% main_data.quad_fit{2} = to_freq_unc_quad;
% main_data.shots = tot_num_shots;
% save([anal_opts.global.out_dir,'main_data.mat'],'main_data')
% save([anal_opts.global.out_dir,'drift_data.mat'],'drift_data')
