%Script that scrapes the analysed data from dirs (currently messy but works)
clear all
%setup directories you wish to loop over
loop_config.dir = {'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20181205_baseline_nuller_on_always\',
        'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20181204_baseline_1\',
        'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20181201_filt_skew_neg111ghz\',
        'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20181203_filt_skew_pos50ghz\',
        'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20181202_filt_skew_pos110ghz\',
        'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20181123_3_filt_align_dep_36.8um\',
        'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20181122_alignment_dep_34_5\',
        'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20181011_to_drift_2\',
        'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20181120_filt_dep_3filt\',
        'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20181202_filt_skew_neg50ghz\',
        'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20181123_3_filt_align_dep_44.9_um\',
        'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20181026_wp_out_stab\',
        'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20181026_wp_out_stab2\',
        'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20181123_3_filt_align_dep_31um\',
        'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20181203_filt_skew_pos50ghz_bad_setpt\',
        'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20190110_baseline_to_1'};
    loop_config.set_pt = [8.0, 3.0, 8.0, 3.0, 3.0, 5.0, 5.0, 2.0, 5.0, 1.0];
selected_dirs = 1:numel(loop_config.dir); %which files to loop over (currently all)
TO_st_pt = 7.257355*1e14;
drift_data_compiled.to_val{:,1}=[];
drift_data_compiled.to_val{:,2}=[];
drift_data_compiled.to_time=[];
drift_data_compiled.atom_num=[];
drift_data_compiled.grad{:,1}=[];
drift_data_compiled.grad{:,2}=[];
drift_data_compiled.avg_coef = [];
drift_data_compiled.avg_coef_unc = [];
main_data_compiled.lin_fit{1} = [];
main_data_compiled.lin_fit{2} = [];
main_data_compiled.quad_fit{1} = [];
main_data_compiled.quad_fit{2} = [];
main_data_compiled.shots = [];
main_data_compiled.time = [];
break_idx = [];
set_pts = [];
sfigure(3712);
clf
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])
subplot(2,1,1)
hold on
set(gcf,'color','w')
for loop_idx=selected_dirs
    current_dir = loop_config.dir{loop_idx};
    out_dirs=dir([current_dir,'out\']);
    offset=0; %to make sure file contains the mat files
    most_recent_dir=out_dirs(end-offset,1);
    %runs through all the out put dirs for a given run and looks for saved data, if none is there
    %skips that data dir
    try
        while ~and(isfile([current_dir,'out\',most_recent_dir.name,'\main_data.mat']),isfile([current_dir,'out\',most_recent_dir.name,'\drift_data.mat']))
            offset = offset + 1;
            most_recent_dir=out_dirs(end-offset,1);
            check = drift_data.avg_coef; %check if it has the avg coefs update
        end
    catch
        fprintf('\n dir: %s didnt work \n',current_dir)
        continue
    end
    load([current_dir,'out\',most_recent_dir.name,'\main_data.mat'])
    load([current_dir,'out\',most_recent_dir.name,'\drift_data.mat'])
    
    %scrapes set points if you want to redo some analysis
    try
        %fprintf('\n dir: %s , st pt = %u \n',current_dir,main_data.set_pt)
        set_pts = [set_pts,main_data.set_pt];
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
    
%     drift_data_compiled.model=to_seg_fits.fit_trimmed.model;

    main_data_compiled.lin_fit{1} = [main_data_compiled.lin_fit{1};main_data.lin_fit{1}];
    main_data_compiled.lin_fit{2} = [main_data_compiled.lin_fit{2};main_data.lin_fit{2}];
    main_data_compiled.quad_fit{1} = [main_data_compiled.quad_fit{1};main_data.quad_fit{1}];
    main_data_compiled.quad_fit{2} = [main_data_compiled.quad_fit{2};main_data.quad_fit{2}];
    main_data_compiled.shots = [main_data_compiled.shots;main_data.shots];
    main_data_compiled.time = [main_data_compiled.time;nanmean(drift_data.to_time)];
    
    break_idx = [break_idx;numel(drift_data.to_time)];
    
    to_time=drift_data.to_time;
    to_val=drift_data.to_val{:,1};
    to_unc=drift_data.to_val{:,2};
    
    plot(to_time,(to_val-TO_st_pt)./1e9,'k')
    hold on
    plot(to_time,(to_val-to_unc-TO_st_pt)./1e9,'b.-')
    plot(to_time,(to_val+to_unc-TO_st_pt)./1e9,'b.-')

end
title('Segment Scan data')
xlabel('time since epoch (h)')
ylabel([' TO value - ',num2str(TO_st_pt./1e9),' (GHz)'])
box on
hold off
%% Plot the analysis of whole runs
subplot(2,1,2)
errorbar(main_data_compiled.time,(main_data_compiled.quad_fit{1}-TO_st_pt)./1e9,main_data_compiled.quad_fit{2}./1e9,'kx')
hold on
errorbar(main_data_compiled.time,(main_data_compiled.lin_fit{1}-TO_st_pt)./1e9,main_data_compiled.lin_fit{2}./1e9,'bo')
legend('Quad','Lin')
xlabel('time since epoch (h)')
ylabel([' TO value - ',num2str(TO_st_pt./1e9),' (GHz)'])
title('Quad fit data for whole runs')
set(gcf,'color','w')
ylim([-3 5])
%% Concatonated Drift data, right now is very messy could do with some cleaning, but it does work
%Big plot over times
sfigure(3221);
clf
%sorts data in time
[~,tsort]=sort(drift_data_compiled.to_time,'ascend');
to_time_concat = [];
break_idx_sort = [];
%finds the indx of the transitions between runs
for kk = 1:numel(break_idx)
    break_idx_sort(kk) = find(tsort==sum(break_idx(1:kk)));
end
break_idx_sort = sort(break_idx_sort);
break_idx_sort = [0,break_idx_sort];
to_time=drift_data_compiled.to_time(tsort);
shift_t = to_time(1);
%removes the time between runs to 0.25 hours
for kk = 2:numel(break_idx)+1
    to_time_concat = [to_time_concat; ones(break_idx_sort(kk)-break_idx_sort(kk-1),1).*shift_t];
    if kk == numel(break_idx)+1
        break
    end
    shift_t = shift_t + to_time(break_idx_sort(kk)+1)-to_time(break_idx_sort(kk))-0.25;
end
%plots thet data and the break lines
to_time=drift_data_compiled.to_time(tsort)-to_time_concat;
hold on
for kk = 2:numel(break_idx_sort)-1
    x_break = (to_time(break_idx_sort(kk))+to_time(break_idx_sort(kk)+1))/2;
    plot([x_break,x_break],[-12 12],'k-','LineWidth',0.5)
end
to_val=drift_data_compiled.to_val{:,1};
to_val=to_val(tsort);
to_unc=drift_data_compiled.to_val{:,2};
to_unc=to_unc(tsort);
set(gcf,'color','w')
for kk = 2:numel(break_idx_sort)
    idx_range = (break_idx_sort(kk-1)+1):break_idx_sort(kk);
    plot(to_time(idx_range),(to_val(idx_range)-TO_st_pt)./1e9,'k')
    plot(to_time(idx_range),(to_val(idx_range)-to_unc(idx_range)-TO_st_pt)./1e9,'b.-')
    plot(to_time(idx_range),(to_val(idx_range)+to_unc(idx_range)-TO_st_pt)./1e9,'b.-')
end
hold off
xlabel('time since epoch (h)')
ylabel([' TO value - ',num2str(TO_st_pt./1e9),' (GHz)'])
ylim([-5 5])
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])
%% Plot tune out value against parameters in the analysis (to see if there is any underling corralations)
%Should turn this into a loop, didn't get around to it
sfigure(979);
clf
subplot(4,4,1)
scatter(drift_data_compiled.grad{:,1},(drift_data_compiled.to_val{:,1}-TO_st_pt)./1e9,'kx')
xlabel('Gradient of Signal')
ylabel([' TO value - ',num2str(TO_st_pt./1e9),' (GHz)'])
set(gcf,'color','w')
box on
subplot(4,4,2)
scatter(drift_data_compiled.avg_coef(:,2),(drift_data_compiled.to_val{:,1}-TO_st_pt)./1e9,'kx')
xlabel('Mean Trap Freq (Hz)')
ylabel([' TO value - ',num2str(TO_st_pt./1e9),' (GHz)'])
set(gcf,'color','w')
box on
subplot(4,4,3)
scatter(drift_data_compiled.avg_coef(:,1),(drift_data_compiled.to_val{:,1}-TO_st_pt)./1e9,'kx')
xlabel('Mean Osc Amp')
ylabel([' TO value - ',num2str(TO_st_pt./1e9),' (GHz)'])
set(gcf,'color','w')
box on
subplot(4,4,4)
scatter(drift_data_compiled.atom_num,(drift_data_compiled.to_val{:,1}-TO_st_pt)./1e9,'kx')
xlabel('Atom Number')
ylabel([' TO value - ',num2str(TO_st_pt./1e9),' (GHz)'])
set(gcf,'color','w')
box on
subplot(4,4,5)
scatter(drift_data_compiled.avg_coef(:,3),(drift_data_compiled.to_val{:,1}-TO_st_pt)./1e9,'kx')
xlabel('Phase')
ylabel([' TO value - ',num2str(TO_st_pt./1e9),' (GHz)'])
set(gcf,'color','w')
box on
subplot(4,4,6)
scatter(drift_data_compiled.avg_coef(:,4),(drift_data_compiled.to_val{:,1}-TO_st_pt)./1e9,'kx')
xlabel('Detector offset')
ylabel([' TO value - ',num2str(TO_st_pt./1e9),' (GHz)'])
set(gcf,'color','w')
box on
subplot(4,4,7)
scatter(drift_data_compiled.avg_coef(:,5),(drift_data_compiled.to_val{:,1}-TO_st_pt)./1e9,'kx')
xlabel('x ramp')
ylabel([' TO value - ',num2str(TO_st_pt./1e9),' (GHz)'])
set(gcf,'color','w')
box on
subplot(4,4,8)
scatter(drift_data_compiled.avg_coef(:,6),(drift_data_compiled.to_val{:,1}-TO_st_pt)./1e9,'kx')
xlabel('y ramp')
ylabel([' TO value - ',num2str(TO_st_pt./1e9),' (GHz)'])
set(gcf,'color','w')
box on
subplot(4,4,9)
scatter(drift_data_compiled.avg_coef(:,7),(drift_data_compiled.to_val{:,1}-TO_st_pt)./1e9,'kx')
xlabel('Damping rate')
ylabel([' TO value - ',num2str(TO_st_pt./1e9),' (GHz)'])
set(gcf,'color','w')
box on
subplot(4,4,10)
scatter(drift_data_compiled.avg_coef(:,8),(drift_data_compiled.to_val{:,1}-TO_st_pt)./1e9,'kx')
xlabel('y ramp')
ylabel([' TO value - ',num2str(TO_st_pt./1e9),' (GHz)'])
set(gcf,'color','w')
box on
subplot(4,4,11)
scatter(drift_data_compiled.avg_coef_unc(:,2),(drift_data_compiled.to_val{:,1}-TO_st_pt)./1e9,'kx')
xlabel('unc in osc freq')
ylabel([' TO value - ',num2str(TO_st_pt./1e9),' (GHz)'])
set(gcf,'color','w')
box on
subplot(4,4,12)
scatter(drift_data_compiled.avg_coef_unc(:,1),(drift_data_compiled.to_val{:,1}-TO_st_pt)./1e9,'kx')
xlabel('unc in osc amplitude')
ylabel([' TO value - ',num2str(TO_st_pt./1e9),' (GHz)'])
set(gcf,'color','w')
box on
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])