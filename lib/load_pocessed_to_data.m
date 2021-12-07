function data = load_pocessed_to_data(loop_config)
selected_dirs = 1:numel(loop_config.dir); %which files to loop over (currently all)

% this script will combine the analysed data from multiple folders 


drift_data_compiled.wp.qwp=[];
drift_data_compiled.wp.hwp=[];
drift_data_compiled.to.val=[];
drift_data_compiled.to.unc=[];
drift_data_compiled.to_time=[];
drift_data_compiled.atom_num_probe.val=[];
drift_data_compiled.atom_num_probe.std=[];
drift_data_compiled.grad.val=[];
drift_data_compiled.grad.unc=[];
drift_data_compiled.avg_coef = [];
drift_data_compiled.avg_coef_cal = [];
drift_data_compiled.avg_coef_unc = [];
drift_data_compiled.avg_coef_cal_unc = [];
drift_data_compiled.probe_shots=[];

main_data_compiled.wp.qwp=[];
main_data_compiled.wp.hwp=[];

main_data_compiled.lin.to.val=[];
main_data_compiled.lin.to.unc=[];
main_data_compiled.lin.grad.val=[];
main_data_compiled.lin.grad.unc=[];

main_data_compiled.quad.to.val=[];
main_data_compiled.quad.to.unc=[];
main_data_compiled.quad.grad.val=[];
main_data_compiled.quad.grad.unc=[];

main_data_compiled.shots_probe = [];
main_data_compiled.start_time = [];
main_data_compiled.atom_num.mean =[];
main_data_compiled.atom_num.std=[];

main_data_compiled.grad = [];
main_data_compiled.freq = [];
main_data_compiled.scan_num = [];
%break_idx = [];
main_data_compiled.set_pt = [];
for loop_idx=selected_dirs
    current_dir = loop_config.dir{loop_idx};
    fprintf('importing data from \n %s \n',current_dir)
    if ~strcmp(current_dir(end),'\')
        current_dir = [current_dir,'\'];
    end
    %if exist('Z:\EXPERIMENT-DATA\2018_Tune_Out_V2','file')~=7, error(sprintf('dir \n %s \n does not exist',current_dir)) ,end
    out_dirs=dir(fullfile(current_dir,'out'));
    out_dirs=out_dirs(3:end);
    out_dirs=out_dirs(cat(1,out_dirs.isdir));
    if size(out_dirs,1)==0
        warning(sprintf('dir \n %s \n does not contain any out dirs',current_dir)), 
    else 
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
                out_instance_folder_path=fullfile(current_dir,'out',out_dirs(folder_index).name,'done.txt');
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
            %parse the folder name into qwp hwp angles
            wp_angles=parse_folder_name_into_wp_angles(current_dir);

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

            drift_data_compiled.to.val=cat(1,drift_data_compiled.to.val,to_fit_seg.fit_trimmed.to_freq.val);
            drift_data_compiled.to.unc=cat(1,drift_data_compiled.to.unc,to_fit_seg.fit_trimmed.to_freq.unc);
            
            %hacky way to get the same size in the qwp,hwp data
            drift_data_compiled.wp.qwp=cat(1,drift_data_compiled.wp.qwp,wp_angles.qwp_ang.*ones(size(to_fit_seg.fit_trimmed.to_freq.val)));
            drift_data_compiled.wp.hwp=cat(1,drift_data_compiled.wp.hwp,wp_angles.hwp_ang.*ones(size(to_fit_seg.fit_trimmed.to_freq.val)));

            drift_data_compiled.to_time=cat(1,drift_data_compiled.to_time,to_fit_seg.to_time); %mean of segment start, end time
            drift_data_compiled.atom_num_probe.val=cat(1,drift_data_compiled.atom_num_probe.val,to_fit_seg.atom_num(:,1));
            drift_data_compiled.atom_num_probe.std=cat(1,drift_data_compiled.atom_num_probe.std,to_fit_seg.atom_num(:,2));
            % annoyingly i didnt output how many shots in each segment, but this does the trick
            drift_data_compiled.probe_shots=cat(1,drift_data_compiled.probe_shots,cellfun(@(x) size(x,1),to_fit_seg.xdat));
            
            if ~isequal(size(drift_data_compiled.atom_num_probe.val),size( drift_data_compiled.probe_shots))
                warning('size not equal')
            end
            
            drift_data_compiled.grad.val=cat(1,drift_data_compiled.grad.val,to_fit_seg.fit_trimmed.slope.val);
            drift_data_compiled.grad.unc=cat(1,drift_data_compiled.grad.unc,to_fit_seg.fit_trimmed.slope.unc);

            main_data_compiled.wp.qwp=cat(1,main_data_compiled.wp.qwp,wp_angles.qwp_ang);
            main_data_compiled.wp.hwp=cat(1,main_data_compiled.wp.hwp,wp_angles.hwp_ang);

            main_data_compiled.lin.to.val=cat(1, main_data_compiled.lin.to.val,to_fit_all.fit_trimmed.to_freq(1).val);
            main_data_compiled.lin.to.unc=cat(1, main_data_compiled.lin.to.val,to_fit_all.fit_trimmed.to_freq(1).unc);
            main_data_compiled.lin.grad.val=cat(1,main_data_compiled.lin.grad.val,to_fit_all.fit_trimmed.slope.val(1));
            main_data_compiled.lin.grad.unc=cat(1,main_data_compiled.lin.grad.unc,to_fit_all.fit_trimmed.slope.unc(1));

            main_data_compiled.quad.to.val=cat(1, main_data_compiled.quad.to.val,to_fit_all.fit_trimmed.to_freq(2).val);
            main_data_compiled.quad.to.unc=cat(1, main_data_compiled.quad.to.unc,to_fit_all.fit_trimmed.to_freq(2).unc);
            main_data_compiled.quad.grad.val=cat(1,main_data_compiled.quad.grad.val,to_fit_all.fit_trimmed.slope.val(2));
            main_data_compiled.quad.grad.unc=cat(1,main_data_compiled.quad.grad.unc,to_fit_all.fit_trimmed.slope.unc(2));

            main_data_compiled.shots_probe = cat(1,main_data_compiled.shots_probe,to_fit_all.num_shots);
            main_data_compiled.shots_all = cat(1,main_data_compiled.shots_probe,numel(to_fit_all.fit_mask));
            main_data_compiled.shots_all = cat(1,main_data_compiled.shots_probe,numel(to_fit_all.fit_mask));

            main_data_compiled.start_time = cat(1,main_data_compiled.start_time,to_fit_all.start_time);
            main_data_compiled.atom_num.mean = cat(1,main_data_compiled.atom_num.mean,to_fit_all.atom_num.mean);
            main_data_compiled.atom_num.std = cat(1, main_data_compiled.atom_num.std,to_fit_all.atom_num.std);


            if ~isequal(size(drift_data_compiled.to.val),size(drift_data_compiled.wp.qwp))
                error('things are not the same size') 
            end

        end
    end



end

data.drift = drift_data_compiled;
data.main = main_data_compiled;

end

%%break_idx = [break_idx;numel(drift_data.to_time)];

% processing steps


%         temp = cell2mat(drift_data.avg_coefs);
%         c_uncs = reshape(temp(:,2),8,numel(drift_data.avg_coefs))';
%         c_vals = reshape(temp(:,1),8,numel(drift_data.avg_coefs))';
% 
%         drift_data_compiled.avg_coef = [drift_data_compiled.avg_coef;c_vals];
%         drift_data_compiled.avg_coef_unc = [drift_data_compiled.avg_coef_unc;c_uncs];
% 
%         try
%         R = drift_data.avg_coefs_cal;
%         R = R(~cellfun('isempty',R));
%         temp = cell2mat(R);
%         c_uncs = reshape(temp(:,2),8,numel(R))';
%         c_vals = reshape(temp(:,1),8,numel(R))';
% 
%         drift_data_compiled.avg_coef_cal = [drift_data_compiled.avg_coef_cal;c_vals];
%         drift_data_compiled.avg_coef_cal_unc = [drift_data_compiled.avg_coef_cal_unc;c_uncs];
%         catch
%         end

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
