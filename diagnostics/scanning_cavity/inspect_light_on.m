function data = inspect_light_on(ld,config)  
% close all

%% Process & display
   

%     all_peaks = light_on_data.vals{:};
%     nz_peaks = all_peaks(all_peaks>0);
%     all_backs = light_on_data.all_vals(:,2);
%     all_pts = light_on_data.raw;
%     alldat = cell_horzcat(all_pts');
%     max_height = max(alldat{1});
%     peak_thresh = 0.01*max_height;
%     all_sort = sort(alldat{1});
% 
%     sort_pk_mask = all_sort>peak_thresh;
%     raw_sort_pks = all_sort(sort_pk_mask);
%     raw_sort_bkg = all_sort(~sort_pk_mask);
%     N_below = length(raw_sort_bkg);
    num_files = numel(ld);
    leb_ratios = cell(num_files,1);
    varf_ratios = cell(num_files,1);
    varf_back = cell(num_files,1);
    varf_axis = cell(num_files,1);
    L_back = cell(num_files,1);
    for nn = 1:num_files
        
        % Lebesgue method
        leb_ratios{nn} = cellfun(@(x) x.lebesgue.P_ratio,ld{nn});
        L_back{nn} = cellfun(@(x) x.lebesgue.P_back,ld{nn});     
        varf_ratios{nn} = cellfun(@(x) x.lebesgue.varf_ratios,ld{nn},'UniformOutput',false);
        varf_axis{nn} = cellfun(@(x) x.lebesgue.varf_axis,ld{nn},'UniformOutput',false);
        varf_back{nn} = cellfun(@(x) x.lebesgue.varf_back,ld{nn},'UniformOutput',false);
    end
    
    data.L_ratios = leb_ratios;
    L_ratios = cell_horzcat(leb_ratios');
    L_ratios = L_ratios{1};
    L_ratios = L_ratios(~isnan(L_ratios) & ~isinf(L_ratios));
    data.all_ratios = L_ratios;
    data.L_mean = nanmean(data.all_ratios);
    data.L_std = nanstd(data.all_ratios);
    data.varf_ratios = varf_ratios;
    data.varf_axis = varf_axis;
    L_back = cell_horzcat(L_back');
    L_back = L_back{1};
    L_back = L_back(~isnan(L_back) & ~isinf(L_back));
    data.L_back = L_back;
    data.varf_back = varf_back;


%% PLOTS






end