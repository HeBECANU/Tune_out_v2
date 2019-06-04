function data = load_pocessed_to_data(loop_config)
selected_dirs = 1:numel(loop_config.dir); %which files to loop over (currently all)

% this script will combine the analysed data from multiple folders 
drift_data_compiled.to.val=[];
drift_data_compiled.to.unc=[];
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
    %if exist('Z:\EXPERIMENT-DATA\2018_Tune_Out_V2','file')~=7, error(sprintf('dir \n %s \n does not exist',current_dir)) ,end
    out_dirs=dir(fullfile(current_dir,'out'));
    out_dirs=out_dirs(3:end);
    out_dirs=out_dirs(cat(1,out_dirs.isdir));
    if size(out_dirs,1)==0, error(sprintf('dir \n %s \n does not contain any out dirs',current_dir)), end 
    % convert the folder name (iso time) to posix time
    time_posix=cellfun(@(x) posixtime(datetime(datenum(x,'yyyymmddTHHMMSS'),'ConvertFrom','datenum')),{out_dirs.name});
    [~,sort_idx]=sort(time_posix,'descend');
    out_dirs=out_dirs(sort_idx);
    looking_for_data_dir=1;
    folder_index=1;
    %runs through all the out put dirs for a given run and looks for saved data, if none is there
    %skips that data dir
    while looking_for_data_dir
        try
            out_instance_folder_path=fullfile(current_dir,'out',out_dirs(folder_index).name,'done.txt')
            if (exist(out_instance_folder_path,'file') || ...
                exist(out_instance_folder_path,'file')) && ...
                exist(out_instance_folder_path,'file')
                looking_for_data_dir=0;
            else
                folder_index=folder_index+1;
                if folder_index>numel(out_dirs) %if beyon the end of the folder list return nan;
                     looking_for_data_dir=0;
                     folder_index=nan;
                     warning('did not find a valid output direcory for folder %s',current_dir)
                end
            end
            %~and(isfile([current_dir,'out\',most_recent_dir.name,'\main_data.mat']),isfile([current_dir,'out\',most_recent_dir.name,'\drift_data.mat']))
            %offset = offset + 1;
            %most_recent_dir=out_dirs(end-offset,1);
            %check = drift_data.avg_coef; %check if it has the avg coefs update
        
        catch e
            fprintf('\n dir: %s didnt work \n',current_dir)
            msgText = getReport(e)
            continue
        end
    end
    if ~isnan(folder_index)
        load(fullfile(current_dir,'out',out_dirs(folder_index).name,'data_results.mat'),'to_fit_seg','to_fit_all','anal_opts')
        % now do some serious data plumbing
        %scrapes set points if you want to redo some analysis
        try
            %fprintf('\n dir: %s , st pt = %u \n',current_dir,main_data.set_pt)
            main_data_compiled.set_pt = cat(1,main_data_compiled.set_pt,anal_opts.ai_log.pd.set_probe);
        catch
            fprintf('\n dir: %s , no st pt record \n',current_dir)
        end

        %append to main structure
        drift_data_compiled.to.val{end+1}=to_fit_seg.fit_trimmed.to_freq.val;
        drift_data_compiled.to.unc{end+1}=to_fit_seg.fit_trimmed.to_freq.unc;
        drift_data_compiled.to_time(end+1)=to_fit_seg.to_time;
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

        main_data_compiled.lin_fit{1} = [main_data_compiled.lin_fit{1};to_fit_seg.fit_trimmed.to_freq.val];
        main_data_compiled.lin_fit{2} = [main_data_compiled.lin_fit{2};to_fit_seg.fit_trimmed.to_freq.unc];
        main_data_compiled.quad_fit{1} = [main_data_compiled.quad_fit{1};main_data.quad_fit{1}];
        main_data_compiled.quad_fit{2} = [main_data_compiled.quad_fit{2};main_data.quad_fit{2}];
        main_data_compiled.shots = [main_data_compiled.shots;main_data.shots];
        main_data_compiled.time = [main_data_compiled.time;nanmean(drift_data.to_time)];
        main_data_compiled.scan_num = [main_data_compiled.scan_num;numel(drift_data.atom_num)];
        main_data_compiled.grad = [main_data_compiled.grad;nanmean(grad_temp(:,1))];
        main_data_compiled.freq = [main_data_compiled.freq;nanmean(c_vals(:,2))];

        %break_idx = [break_idx;numel(drift_data.to_time)];
    end
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
