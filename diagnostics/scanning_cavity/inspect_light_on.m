function data = inspect_light_on(ld,config)  
% close all

%% Process & display
   
    data = [];

    num_files = numel(ld);
    pv_peaks = cell(num_files,1);
    pv_backs = cell(num_files,1);
    pv_ratios = cell(num_files,1);
    data.pv_stats = cell(num_files,1);
    
    if config.pv_method
        for nn = 1:num_files
            % Get the things
            data.pv_peaks{nn} = cellfun(@(x) x.pv_data.peaks,ld{nn},'UniformOutput',false);
            data.pv_backs{nn} = cellfun(@(x) x.pv_data.backs,ld{nn},'UniformOutput',false);     
            data.pv_ratios{nn} = cellfun(@(x) x.pv_data.ratio,ld{nn},'UniformOutput',false); 
        end
        
        all_ratios = cell_horzcat(data.pv_ratios');
        all_ratios = all_ratios{1};
        all_ratios = cell_vertcat(all_ratios');
        all_ratios = all_ratios{1};
        all_ratios = all_ratios(~isnan(all_ratios) & ~isinf(all_ratios));
        data.pv_ratios{nn} = all_ratios;
%         data.pv_stats{nn} = [mean(all_ratios),std(all_ratios)];
        if config.insp_hists
            sfigure(1337);
            binedges = linspace(-3e-3,10e-3,20);
            histogram(all_ratios,binedges)
            hold on        
        end
        pv_all = cell_horzcat(data.pv_ratios);
        pv_all = cell_vertcat(pv_all{1}');
        data.pv_all = pv_all{1};
        data.pv_stats = [mean(data.pv_all),std(data.pv_all),length(data.pv_all)];
    end
    
    if config.lebesgue_method
        
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
        data.L_num = length(data.all_ratios);
        data.varf_ratios = varf_ratios;
        data.varf_axis = varf_axis;
        L_back = cell_horzcat(L_back');
        L_back = L_back{1};
        L_back = L_back(~isnan(L_back) & ~isinf(L_back));
        data.L_back = L_back;
        data.varf_back = varf_back;
    end
    


%% PLOTS






end