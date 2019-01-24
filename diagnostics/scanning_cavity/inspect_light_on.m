function inspect_light_on(light_on_data,config)  
%% Process & display


    all_peaks = light_on_data.all_vals(:,1);
    nz_peaks = all_peaks(all_peaks>0);
    all_backs = light_on_data.all_vals(:,2);
    all_pts = light_on_data.raw_data;
    alldat = cell_horzcat(all_pts');
    max_height = max(alldat{1});
    peak_thresh = 0.01*max_height;
    all_sort = sort(alldat{1});

    sort_pk_mask = all_sort>peak_thresh;
    raw_sort_pks = all_sort(sort_pk_mask);
    raw_sort_bkg = all_sort(~sort_pk_mask);
    N_below = length(raw_sort_bkg);
    figure();



    % Look at stuff!
    subplot(4,1,1)
    plot(nz_peaks,'.')
    title('Peak powers v time')
    ylim([0,max(nz_peaks)])

    subplot(4,1,2)
    plot(all_backs,'.')
    title('Background v time')
    % ylim([0,max(all_backs)])

    subplot(4,4,9)
    histogram(nz_peaks)
    title('Peak powers')

    subplot(4,4,10)
    histogram(all_backs)
    title('Background powers')
    suptitle(sprintf('ALL files in %s',config.dirname))
    
        % Show the Lebesgue technique

    % Inspect the ratio statistics
    subplot(4,4,11)
    histogram(light_on_data.dir_ratios,50);
    title('Peak/back ratios')
end