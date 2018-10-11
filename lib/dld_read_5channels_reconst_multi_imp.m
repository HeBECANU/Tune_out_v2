function[three_channel_output_combined_sorted]= dld_read_5channels_reconst_multi_imp(filename_call,normalise_time_flag_call,reconst_4_corners_and_mcp_flag_call,reconst_4_corners_nomcp_flag_call,reconst_3_corners_flag_call)

%Reads raw dld data from 4 corners and MCP and converts to t,x,y
%BRYCE CHANGE LOG

%changed from dlmread to the tdc_importer which uses strong typing to speed
%up reads
%runs 7-10% faster than basic

%added 5ch reconstruction tag which speeds things up 45%


%timing notes
%through this code i have refered to a benchmark which takes 0.28s to run
%(with only 4 corrner)
%and often quote the time to run some code relative to that

% tic

format long e

% %variables which will probably go in function call
% filename = 'real_mcp_pulse_testdata1_reconst';
% normalise_time_flag = 0;    %if 1 subtracts hardware trigger time from all event times
% reconst_4_corners_nomcp_flag = 1;   %if 1 reconstructs events with 4 corners but no matching MCP pulse
% reconst_3_corners_flag = 1; %if 1 reconstructs events with only 3 corners (both with and without MCP pulse)
% %%%%%%%%%%%%%%%%

filename = filename_call; %must be a string, without the .txt if the file has this extension
normalise_time_flag = normalise_time_flag_call; %boolean variable: if 1 subtracts hardware trigger time from all event times
reconst_4_corners_nomcp_flag = reconst_4_corners_nomcp_flag_call;   %boolean variable: if 1 reconstructs events with 4 corners but no matching MCP pulse
reconst_3_corners_flag = reconst_3_corners_flag_call; %boolean variable: if 1 reconstructs events with only 2/3 corners (both with and without MCP pulse)

max_group_time = 3400; %maximum time between first and last event in group
dead_time = 400; %time after 1 group to wait before looking for next group
tsum = 3200;
tolerance = 200;    %tolerance in bins

bin_time = 25e-12;                                          %% DLD bin size of 25 ps
v_perp_x = 5.2632e+005;
v_perp_y = 5.2632e+005;

filename_input = [filename,'.txt'];
%BRYCE CHANGE
%saves ~10%
dld_output_raw =tdc_importer(filename_input);
%dld_output_raw = dlmread(filename_input, ',');

%dld_output_raw = dld_output_raw(1:1000,:);   %only used during debugging to reduce file size.  Delete at other times.

%bryce change single step
% number_triggers_matrix = size(dld_output_raw);
% num_trigs = number_triggers_matrix(1);  %number of clicks on TDC
%no time improvement
 num_trigs = size(dld_output_raw,1);


[~, NewRowNumber] = sort(dld_output_raw(:,2));

dld_output_sorted = dld_output_raw(NewRowNumber,:);

t_zero = dld_output_raw(1,2); %time of triggerpulse

if normalise_time_flag == 1
    dld_output_sorted(:,2) = dld_output_sorted(:,2)-t_zero;
end

%dld_output_sorted;

%tsum reconstructing values.  tsumy correcting factor = a_y*x^3+b_y*x^2+c_y*x +d_y 
% Then tsumy = tsum + tsumy correcting factor
a_x = 0.000897;%2.0434e-18;
b_x = 0.021122;%3.6569e-12;
c_x = -2.8726;%-3.7798e-5;
d_x = 54.091;

a_y = 0.0012839;%2.9248e-18;
b_y = 0.014443;%2.5006e-12;
c_y = -2.6612;%-3.60163-5;
d_y = -104.81;








%-------------------------5 hit reconstruction----------------------------
%the 5 hit part takes a huge 140ms which is massive 50 percent of the time
%to run
%timing notes


%BRYCE CHANGE
%added in this toggle as we dont really use the plate pulses
if reconst_4_corners_and_mcp_flag_call
    
    %%%%%%%First, read off hits with MCP pulse and exactly 1 hit on each corner

    event_counter_5trigs = 0;   %number of events with exactly 5 triggers
    five_channel_events_5trigs = NaN(floor(num_trigs/5),5);    %overestimates size, but should still save time

    %dld_raw_print = dld_output_sorted(1:867,:);

    %now sort out events with exactly 5 triggers
    for n=2:num_trigs-5
        check1 = dld_output_sorted(n+4,2)-dld_output_sorted(n,2);   %check if there are 5 events with close trigger times

        if check1 < max_group_time
            check2 = dld_output_sorted(n,2)-dld_output_sorted(n-1,2);   %check this isn't the 2nd event in a multihit group

            if check2 > max_group_time
                check3 = dld_output_sorted(n+5,2)-dld_output_sorted(n,2);   %ensures there's only 5 events in group (implement multihit stuff later if necessary)
                if check3 > check1 + dead_time  %was max_group_time instead of check1, but now this checks based on spread of events in group

                    group_matrix_1 = [dld_output_sorted(n:n+4,1),dld_output_sorted(n:n+4,2)];    %pull out each event as an individual matrix

                    [~, NewRowNum_gm1] = sort(group_matrix_1(:,1));
                    group_matrix_1_sorted = group_matrix_1(NewRowNum_gm1,:); %orders matrix by trigger channel

                    %consistency check, to ensure channels 0-3 and 7 present
                    if group_matrix_1_sorted(1,1) ==0; %Ensure first channel is 0
                        magic_number = group_matrix_1_sorted(2,1)*group_matrix_1_sorted(3,1)*group_matrix_1_sorted(4,1)*group_matrix_1_sorted(5,1);
                        if magic_number==42  %ensure next 4 are 1,2,3,7
                            event_counter_5trigs = event_counter_5trigs + 1;    %increment event counter

                            five_channel_events_5trigs(event_counter_5trigs,1) = group_matrix_1_sorted(5,2);   %mcp pulse (channel 7) is first column in 5_channel matrix
                            five_channel_events_5trigs(event_counter_5trigs,2:5) = group_matrix_1_sorted(1:4,2);   %channels 0-3 are next 4 columns

                            dld_output_sorted(n:n+4,:) = NaN; % nan's the data read to eliminate later
                        else
                            continue
                        end

                    else
                        continue
                    end
                else
                    continue
                end
            else
                continue
            end
        else
            continue
        end
    end

    dld_output_sorted_rowstokeep =~ isnan(dld_output_sorted(:,1));
    dld_output_sorted = dld_output_sorted(dld_output_sorted_rowstokeep,:);   % Removes already read hits from input data matrix

    five_channel_events_5trigs_rowstokeep =~ isnan(five_channel_events_5trigs(:,1));
    five_channel_events_5trigs = five_channel_events_5trigs(five_channel_events_5trigs_rowstokeep,:); %  clear extra pre-allocated zeros in matrix

    % tcalc_x = (five_channel_events_5trigs(:,2)+five_channel_events_5trigs(:,3)-tsum)*.5;
    %
    % tcalc_y = (five_channel_events_5trigs(:,4)+five_channel_events_5trigs(:,5)-tsum)*.5;
    %
    % average_tcalc = (tcalc_x(:)+tcalc_y(:))/2;
    % t_mcp_pulse = five_channel_events_5trigs(:,1);
    % t_diff =  t_mcp_pulse - average_tcalc;
    %calc_mean_t_diff = round(mean(t_diff));   %Seems to be 3400....

    % average_tcalc = (tcalc_x(:)+tcalc_y(:))/2;        %%%%%%%USE LATER
    % t_mcp_pulse = five_channel_events_5trigs(:,1);
    % t_diff_y =  t_mcp_pulse - tcalc_y(:);
    % calc_mean_t_diff_y = round(mean(t_diff_y))
    % t_diff_x =  t_mcp_pulse - tcalc_x(:);
    % calc_mean_t_diff_x = round(mean(t_diff_x))

    %calc_mean_t_diff = round(mean(t_diff));

    mean_t_diff = 3404; %so we don't calculate it each time...
    mean_t_diff_x = 3363;   %calculated from 106868 5 channel hits in real_mcp_pulse_testdata1.txt data file 10/3/10
    x_t_diff_correcting_factor = mean_t_diff_x - mean_t_diff;
    mean_t_diff_y = 3445;   %calculated from 106868 5 channel hits in real_mcp_pulse_testdata1.txt data file 10/3/10
    y_t_diff_correcting_factor = mean_t_diff_y - mean_t_diff;

    five_channel_events_5trigs(:,1) = five_channel_events_5trigs(:,1)-mean_t_diff;  %subtract time off: t channel is now consistent with other data extraction methods

    %five_channel_events_5trigs

    %difference_tcalcxy = (tcalc_x(:)-tcalc_y(:))



    %%%%%%%%%%%% Now convert 5 channels to t,x,y %%%%%%%%%%%%%%%%


    howmanyrows=size(five_channel_events_5trigs);               %% Number of successful counts
    three_channel_output_5trigs=zeros(howmanyrows(1),3);               %% Initialise matrix
    Num_hits_5trigs = howmanyrows(1);


    three_channel_output_5trigs(:,1)=five_channel_events_5trigs(:,1)*bin_time;
    three_channel_output_5trigs(:,2)=(five_channel_events_5trigs(:,2)-five_channel_events_5trigs(:,3))*v_perp_x*bin_time;
    three_channel_output_5trigs(:,3)=(five_channel_events_5trigs(:,4)-five_channel_events_5trigs(:,5))*v_perp_y*bin_time;
    % three_channel_output_5trigs(:,4)=(-2*five_channel_events_5trigs(:,1)+five_channel_events_5trigs(:,2)+five_channel_events_5trigs(:,3));  %tsum_x
    % three_channel_output_5trigs(:,5)=(-2*five_channel_events_5trigs(:,1)+five_channel_events_5trigs(:,4)+five_channel_events_5trigs(:,5));  %tsum_y
    %three_channel_output_5trigs(:,6)=(five_channel_events_5trigs(:,2)+five_channel_events_5trigs(:,3)-2*(five_channel_events_5trigs(:,1)-mean_t_diff_x))*v_perp_x*bin_time;%tof_x calculated from x1 and x2
    %three_channel_output_5trigs(:,7)=(five_channel_events_5trigs(:,4)+five_channel_events_5trigs(:,5)-2*(five_channel_events_5trigs(:,1)-mean_t_diff_y))*v_perp_y*bin_time;%tof_y calculated


    %three_channel_output_5trigs;
    %%%%%%At this point have removed and processed all unique 5 channel hits
    %%%%%%Now go through and find all hits with multi hits on 1 or more
    %%%%%%channels, but only 1 MCP pulse

    number_triggers_left_matrix = size(dld_output_sorted);
    num_trigs_left = number_triggers_left_matrix(1);  %number of clicks on TDC

    max_group_time_multi = 20000;    %was 3800.  Needs to be 1e9 to get all...

    event_counter_5trigs_multi = 1;   %number of events with exactly 5 triggers
    five_channel_events_5trigs_multi = NaN(floor(num_trigs_left/6),5);    %overestimates size, but should still save time

    for m = 2:num_trigs_left-6
        check1 = dld_output_sorted(m+4,2)-dld_output_sorted(m,2);   %check if there are 6 events with close trigger times.  There may be more

        if check1 < max_group_time_multi
            multi_group_hits_counter = 5;
            for q = 5:num_trigs_left

                check2 = q+m;   %this just ensures we don't try to access a matrix cell in dld_output_sorted out of range
                if check2 > num_trigs_left
                    break
                end

                check3 = dld_output_sorted(m+q,2)-dld_output_sorted(m,2);   %check if next cell is within range
                if check3 > max_group_time_multi
                    break
                end

                multi_group_hits_counter = multi_group_hits_counter+1;  %since we've checked that hit is within range, increment counter

            end

            multi_group_hits_raw = dld_output_sorted(m:m+multi_group_hits_counter-1,:);

            multi_group_hits_location_matrix = zeros(multi_group_hits_counter,5);
            multi_group_hits_time_matrix = zeros(multi_group_hits_counter,5);

            for p = 1:multi_group_hits_counter
                if multi_group_hits_raw(p,1) == 7   %time channel is 7, hence is tricky to deal with...
                    multi_group_hits_location_matrix(p,1) = 1;  %index location of hit within group
                    multi_group_hits_time_matrix(p,1) = multi_group_hits_raw(p,2); %and then put actual value in
                else
                    multi_group_hits_location_matrix(p,multi_group_hits_raw(p,1)+2) = 1;    %other channels are easy, as 0 -> x1 (column 2) etc
                    multi_group_hits_time_matrix(p,multi_group_hits_raw(p,1)+2) = multi_group_hits_raw(p,2);
                end

            end

            multi_group_hits_time_matrix;
            multi_group_hits_location_matrix;

            multi_group_no_hits_matrix = sum(multi_group_hits_location_matrix,1);   %1D array with number of hits on each channel in group
            min_hits_single_channel = min(multi_group_no_hits_matrix);

            %         if multi_group_no_hits_matrix(1,1) == 0    %we're only interested in events with mcp pulses for now
            %             continue
            %         end

            if min_hits_single_channel == 0 %Not reconstructing events at this stage
                continue
            end

            multi_group_hits_time_matrix_sorted = sort(multi_group_hits_time_matrix,'descend');    %puts zeros at bottom

            for t_count = 1:multi_group_no_hits_matrix(1)   %loop over number of mcp triggers
                good_hit = 0;
                for x1_count = 1:multi_group_no_hits_matrix(2)  %loop over number of hits on x1
                    for x2_count = 1:multi_group_no_hits_matrix(3)  %loop over hits on x2
                        check1 = multi_group_hits_time_matrix_sorted(x1_count,2)+multi_group_hits_time_matrix_sorted(x2_count,3)-2*(multi_group_hits_time_matrix_sorted(t_count,1)-mean_t_diff)-tsum; %x1+x2-2t = tsum within tolerance
                        if abs(check1) < tolerance
                            five_channel_events_5trigs_multi(event_counter_5trigs_multi,2) = multi_group_hits_time_matrix_sorted(x1_count,2);   %assign x1
                            five_channel_events_5trigs_multi(event_counter_5trigs_multi,3) = multi_group_hits_time_matrix_sorted(x2_count,3);   %assign x2
                            good_hit = good_hit + 1;    %will become >1 if multi events match
                            x1_count_good = multi_group_no_hits_matrix(2) - x1_count+1;   %tag the good count (if multi events matching then it will be overwritten)
                            x2_count_good = multi_group_no_hits_matrix(3)-x2_count+1;    %note that since we re-ordered the matrix in reverse order, have to do the reverse for the indices
                        end

                    end

                end


                %good_hit

                if good_hit ~= 1 %if unable to match x, clear assigned x and move to next t iteration
                    five_channel_events_5trigs_multi(event_counter_5trigs_multi,2) = NaN;
                    five_channel_events_5trigs_multi(event_counter_5trigs_multi,3) = NaN;
                    %crapometer = 1

                    continue
                end

                %now check y

                for y1_count = 1:multi_group_no_hits_matrix(4)  %loop over number of hits on y1
                    for y2_count = 1:multi_group_no_hits_matrix(5)  %loop over hits on y2
                        check2 = multi_group_hits_time_matrix_sorted(y1_count,4)+multi_group_hits_time_matrix_sorted(y2_count,5)-2*(multi_group_hits_time_matrix_sorted(t_count,1)-mean_t_diff)-tsum; %x1+x2-2t = tsum within tolerance. Time at this stage is raw (uncorrected)
                        if abs(check2) < tolerance
                            five_channel_events_5trigs_multi(event_counter_5trigs_multi,4) = multi_group_hits_time_matrix_sorted(y1_count,4);   %assign y1
                            five_channel_events_5trigs_multi(event_counter_5trigs_multi,5) = multi_group_hits_time_matrix_sorted(y2_count,5);   %assign y2
                            good_hit = good_hit + 1;    %will become > 2 if multi events match
                            y1_count_good = multi_group_no_hits_matrix(4) - y1_count + 1;   %tag the good count (if multi events matching then it will be overwritten)
                            y2_count_good = multi_group_no_hits_matrix(5) - y2_count + 1;
                        end

                    end

                end



                if good_hit == 2
                    five_channel_events_5trigs_multi(event_counter_5trigs_multi,1) = multi_group_hits_time_matrix_sorted(t_count,1);    %assign t
                    event_counter_5trigs_multi = event_counter_5trigs_multi + 1;    %increment hit counter
                    t_count_good = t_count;   %tag the good count (if multi events matching then it will be overwritten)

                    t_hit_array_indices = find(multi_group_hits_location_matrix(:,1));
                    raw_data_remove_t = t_hit_array_indices(t_count_good)+m;   %t event to remove from raw data matrix
                    dld_output_sorted(raw_data_remove_t,2) = NaN; %and now remove it

                    x1_hit_array_indices = find(multi_group_hits_location_matrix(:,2));
                    raw_data_remove_x1 = x1_hit_array_indices(x1_count_good)+m;   %x1 event to remove from raw data matrix
                    dld_output_sorted(raw_data_remove_x1,2) = NaN; %and now remove it

                    x2_hit_array_indices = find(multi_group_hits_location_matrix(:,3));
                    raw_data_remove_x2 = x2_hit_array_indices(x2_count_good)+m;   %x2 event to remove from raw data matrix
                    dld_output_sorted(raw_data_remove_x2,2) = NaN; %and now NaN it

                    y1_hit_array_indices = find(multi_group_hits_location_matrix(:,4));
                    raw_data_remove_y1 = y1_hit_array_indices(y1_count_good)+m;   %y1 event to remove from raw data matrix
                    dld_output_sorted(raw_data_remove_y1,2) = NaN; %and now NaN it

                    y2_hit_array_indices = find(multi_group_hits_location_matrix(:,5));
                    raw_data_remove_y2 = y2_hit_array_indices(y2_count_good)+m;   %y2 event to remove from raw data matrix
                    dld_output_sorted(raw_data_remove_y2,2) = NaN; %and now NaN it


                    %m = max([raw_data_remove_t, raw_data_remove_x1, raw_data_remove_x2, raw_data_remove_y1, raw_data_remove_y2])-1; %increment m so we don't go over the NaN'ed data again

                else
                    five_channel_events_5trigs_multi(event_counter_5trigs_multi+1,4) = NaN;   %clear any assigned y values
                    five_channel_events_5trigs_multi(event_counter_5trigs_multi+1,5) = NaN;

                end

            end

            %any events are now recorded

        else
            continue
        end

    end
    





    dld_output_sorted_rowstokeep =~ isnan(dld_output_sorted(:,2));  %clear NaN rows from dld_output_sorted
    dld_output_sorted = dld_output_sorted(dld_output_sorted_rowstokeep,:);

    five_channel_events_5trigs_multi_rowstokeep =~ isnan(five_channel_events_5trigs_multi(:,1));    %remove excess pre-allocated hits
    five_channel_events_5trigs_multi = five_channel_events_5trigs_multi(five_channel_events_5trigs_multi_rowstokeep,:);
    


    %five_channel_events_5trigs_multi

    five_channel_events_5trigs_multi(:,1) = five_channel_events_5trigs_multi(:,1)-mean_t_diff;
    
    


    %%%Now convert to t,x,y

    howmanyrows_multi=size(five_channel_events_5trigs_multi);               %% Number of successful counts
    three_channel_output_5trigs_multi=zeros(howmanyrows_multi(1),3);               %% Initialise matrix
    Num_multi_hits = howmanyrows_multi(1);

    three_channel_output_5trigs_multi(:,1)=five_channel_events_5trigs_multi(:,1)*bin_time;
    three_channel_output_5trigs_multi(:,2)=(five_channel_events_5trigs_multi(:,2)-five_channel_events_5trigs_multi(:,3))*v_perp_x*bin_time;
    three_channel_output_5trigs_multi(:,3)=(five_channel_events_5trigs_multi(:,4)-five_channel_events_5trigs_multi(:,5))*v_perp_y*bin_time;

    %three_channel_output_5trigs_multi
    %three_channel_output_5trigs_multi
    %three_channel_output_5trigs_multi
%BRYCE CHANGE
%added in a toggle for 5ch reconstrucation so had to fake as if it had run
%for the rest of the code
    
else
    five_channel_events_5trigs_multi=zeros(0,3);
    three_channel_output_5trigs_multi=zeros(0,3);
    three_channel_output_5trigs=zeros(0,3);
end

    %%%%%%%%%Now we should have removed every event with exactly 5 unique triggers

    hits_left = size(dld_output_sorted);










%-------------------------4 hit reconstruction----------------------------
%takes 

%%%%%Next stage is to run the old code which finds and sorts events with 4
%%%%%corners only

if reconst_4_corners_nomcp_flag == 1    %only reconstruct 4 corner events if we've flagged that we will

    number_detections_matrix = size(dld_output_sorted);                 %% Will equal [5*n + 1 2] for n detections ideally, the +1 is a master trigger to throw out
    number_detections = number_detections_matrix(1);                    %% Possibly overestimates size since errors will reduce this below what it should be
    number_successes = 0;                                               %% Tally successful hits
    which_row = 0;                                                      %% Index of matrix row to write to
    T_sum = tsum;    %Already defined                                                   %% Is precisely the time taken for signal to travel 8cm, speed 1e6m/s, bins 25e-12
    tolerance_throw = 200;                                              %% Tolerance in what data to throw away
    tolerance_keep = 200;                                               %% Tolerance in what data to keep
    search_no = 36;                                                     %% Seach over 9 (=36/4) complete hits
    T_spread = 0;                                                       %% Spreads in times reconstructed
    T_sum_spread = 0;
    T_sum_x = 0;
    T_sum_y = 0;

    %%%%%%Loop Constants%%%%%
    T_sum_tol_throw = T_sum + tolerance_throw;
    T_sum_tol_keep = T_sum + tolerance_keep;
    two_tolerance_keep = 2*tolerance_keep;

    five_channel_output_4corners = NaN(floor(number_detections/4),5);

    count = int32(0);
    dummy = int8(0);
    count2 = int8(0);
    for count = 2:number_detections                                     %% For each detection event. First hit is a trigger to ignore
        if dld_output_sorted(count,1) == 0                              %% Pick out values of x1 (channel = 0) to match everything to
            x1_val = dld_output_sorted(count,2);                        %% This is x1 to match everything else to
            dummy = 0;                                                  %% Dummy variable to keep search local
            how_many_x = 0;
            %out_of_range_x = [0,0];

            for count2 = 1:search_no                                    %% Now looks at other detections near x1
                dummy = dummy+count2*(-1)^(count2);                     %% Go back 1 row, forward 2, back 3, to keep search efficient

                if (count+dummy) <= 0 ;                                  %% Quit if we go out of range
                    continue
                end
                if (count+dummy) > number_detections
                    continue
                end
                if abs(dld_output_sorted(count+dummy,2) - x1_val) > T_sum;                                  %% Quit if we go out of range
                    %                 out_of_range_x(mod(count2,2)+1) = out_of_range_x(mod(count2,2)+1) + 1;
                    %                 out_of_range_check_x = out_of_range_x(1)*out_of_range_x(2);
                    %                 if out_of_range_check_x ~= 0
                    %                     break
                    %                 end
                    continue
                end

                if dld_output_sorted(count+dummy,1) == 1                %% Pick out values of x2 (channel = 1)
                    x2_val = dld_output_sorted(count+dummy,2);

                    if abs(x1_val - x2_val) < (T_sum_tol_throw) %% See if x2 is close enough to check
                        how_many_x = how_many_x + 1;
                        if how_many_x > 1                               %% Breaks if multihit
                            break
                        else
                            x1 = x1_val;
                            x2 = x2_val;
                            which_row = which_row + 1;
                            tx = 0.5*(x1+x2-T_sum);                     %% Stores the value if everything goes correctly
                            five_channel_output_4corners(which_row,2) = x1 - tx;
                            five_channel_output_4corners(which_row,3) = x2 - tx;

                            %x1_index = count;
                            x2_index = count + dummy;   %if not frozen now, then risk dummy moving on

                            %five_channel_output_4corners(which_row,1) = tx;    %can lead to problems if no y's to match the x. Plus is unnecessary
                        end
                    end
                end
            end                                                         %% At this point, each X is now defined

            if how_many_x == ~1                                         %% If cant find X or is not unique, go to next x1
                continue
            end

            dummy2 = 0;
            how_many_y = 0;

            for count3 = 1:search_no                                    %% For each x1 now match the y1 and y2
                dummy2 = dummy2+count3*(-1)^(count3);
                %out_of_range_y = [0,0];

                if (count+dummy2) <= 0
                    continue
                end
                if (count+dummy2) > number_detections
                    continue
                end
                if abs(dld_output_sorted(count+dummy2,2) - x1_val) > T_sum;                                  %% Quit if we go out of range
                    %                 out_of_range_y(mod(count3,2)+1) = out_of_range_y(mod(count3,2)+1) + 1;
                    %                 out_of_range_check_y = out_of_range_y(1)*out_of_range_y(2);
                    %                 if out_of_range_check_y ~= 0
                    %                     break
                    %                 end
                    continue
                end

                if dld_output_sorted(count+dummy2,1) == 2               %% Pick out values of y1 (channel = 2)
                    y1_val = dld_output_sorted(count+dummy2,2);

                    dummy3 = 0;

                    for count4 = 1:search_no
                        dummy3 = dummy3+count4*(-1)^(count4);

                        if (count+dummy3) <= 0
                            continue
                        end
                        if (count+dummy3) > number_detections
                            continue
                        end
                        if abs(dld_output_sorted(count+dummy3,2) - x1_val) > T_sum;                                  %% Quit if we go out of range
                            continue
                        end

                        if dld_output_sorted(count+dummy3,1) == 3       %% Pick out values of y2 (channel = 3)
                            y2_val = dld_output_sorted(count+dummy3,2);
                            if and(abs(y1_val-y2_val) < (T_sum_tol_keep), abs(x1+x2-y1_val-y2_val) < two_tolerance_keep)
                                how_many_y = how_many_y + 1;
                                if how_many_y > 1
                                    break
                                else
                                    y1 = y1_val;
                                    y2 = y2_val;
                                    ty = 0.5*(y1+y2-T_sum);
                                    t = (tx + ty)/2;
                                    five_channel_output_4corners(which_row,4) = y1 - ty;
                                    five_channel_output_4corners(which_row,5) = y2 - ty;
                                    five_channel_output_4corners(which_row,1) = t;

                                    dld_output_sorted(count,1) = NaN;   %now that we've read the data, NaN the entries in dld_output_sorted we've read for later removal
                                    dld_output_sorted(x2_index,1) = NaN;
                                    dld_output_sorted(count+dummy2,1) = NaN;
                                    dld_output_sorted(count+dummy3,1) = NaN;

                                    number_successes = number_successes + 1;

                                    %                                T_spread = T_spread + abs(tx-ty)^2;     %this is just used for statistics on tsum.  Comment out normally
                                    %                                T_sum_x = T_sum_x +(x1+x2-2*tx);
                                    %                                T_sum_y = T_sum_y +(y1+y2-2*ty);
                                    %                                T_sum_spread = T_sum_spread + abs((x1+x2-2*t) - T_sum)^2 + abs((y1+y2-2*t) - T_sum)^2;
                                end
                            end
                        end
                    end
                end
            end                                                         %% At this point, each Y is now defined
            if how_many_y == 0
                which_row = which_row - 1;
            end
        end
    end

    good_five_channel_output_4corners =~ isnan(five_channel_output_4corners(:,1));
    five_channel_output_4corners = five_channel_output_4corners(good_five_channel_output_4corners,:);

    dld_output_sorted_rowstokeep =~ isnan(dld_output_sorted(:,1));  %clear NaN rows from dld_output_sorted
    dld_output_sorted = dld_output_sorted(dld_output_sorted_rowstokeep,:);

    howmanyrows_4corners=size(five_channel_output_4corners);               %% Number of successful counts
    three_channel_output_4corners=zeros(howmanyrows_4corners(1),3);               %% Initialise matrix
    Num_4corners_hits = howmanyrows_4corners(1);

    three_channel_output_4corners(:,1)=five_channel_output_4corners(:,1)*bin_time;
    three_channel_output_4corners(:,2)=(five_channel_output_4corners(:,2)-five_channel_output_4corners(:,3))*v_perp_x*bin_time;
    three_channel_output_4corners(:,3)=(five_channel_output_4corners(:,4)-five_channel_output_4corners(:,5))*v_perp_y*bin_time;

    %three_channel_output_4corners;


else
    three_channel_output_4corners = [];
end




%-------------------------3 hit reconstruction----------------------------

%
%%%At this point we have removed all hits with 4 corners present.
%%%From now on it will be a re-construction zone

%%%%%%%%%%%% Now try to match up other mcp pulses with signals from 3 corners (or 2)

if reconst_3_corners_flag == 1  %match events with a missing corner only if we've flagged that we will

    size_rem = size(dld_output_sorted);
    num_trigs_remaining = size_rem(1,1);

    dld_output_sorted;


    reconstruct_search_range = 16;
    max_group_time_reconst = 20000; %maximum time after a mcp trigger that we'll look for hits to match up

    event_counter_reconst_mcp = 1;
    five_channel_events_mcp_reconst = nan(ceil(num_trigs_remaining/3),5);
    %three_channel_output_reconst_mcp = zeros(ceil(num_trigs_remaining/3),3);

    for j=1:num_trigs_remaining-3   %loop over remaining hits
        if dld_output_sorted(j,1)==7    %only interested in matching up mcp pulses

            out_of_range_counter = [0,0];
            %BRYCE CHANGE
            %should pre alocate to allow for c conversion
            %unclear if this size is right..
            %hit_group_to_reconst = NaN(10,2);
            num_hits_in_group = 1;
            hit_group_to_reconst(1,:)= dld_output_sorted(j,:);
            k_index_group_array = zeros(17);


            for k=1:reconstruct_search_range
                k_indice = j+((-1)^k)*ceil(k/2);    %check next hit, 1 hit back, 2 forward etc. Need to go back and forth as the MCP pulse isn't always first.

                if k_indice  <2  %skip to next iteration and mark direction of search if too small
                    %out_of_range_counter(mod(k,2)+1) = out_of_range_counter(mod(k,2)+1)+1;
                    continue
                end

                if k_indice > num_trigs_remaining  %skip to next iteration and mark direction of search if range exceeded
                    %out_of_range_counter(mod(k,2)+1) = out_of_range_counter(mod(k,2)+1)+1;
                    continue
                end



                if abs(dld_output_sorted(k_indice ,2)-dld_output_sorted(j,2)) < max_group_time_reconst
                    %blah = 42
                    num_hits_in_group = num_hits_in_group + 1;
                    hit_group_to_reconst(num_hits_in_group,:) = dld_output_sorted(k_indice,:);
                    k_index_group_array(num_hits_in_group) = k_indice;
                else
                    out_of_range_counter(mod(k,2)+1) = out_of_range_counter(mod(k,2)+1)+1;
                    out_of_range_check = out_of_range_counter(1)*out_of_range_counter(2);
                    if  out_of_range_check ~= 0 %if both positive and negative search terms go out of range, then break
                        break
                    end
                end
            end

            if num_hits_in_group <3
                continue   %if less than 3 triggers, we can't reconstruct the event, so move on to next trigger
            end

            %have now removed raw hit matrix from data. Next,need to sort it

            reconst_group_hits_location_matrix = zeros(num_hits_in_group,5);
            reconst_group_hits_time_matrix = zeros(num_hits_in_group,5);
            k_index_group_matrix = zeros(num_hits_in_group,5);

            %num_hits_in_group

            for p = 1:num_hits_in_group
                if hit_group_to_reconst(p,1) == 7   %time channel is 7, hence is tricky to deal with...
                    reconst_group_hits_location_matrix(p,1) = 1;  %index location of hit within group
                    reconst_group_hits_time_matrix(p,1) = hit_group_to_reconst(p,2); %and then put actual value in
                    if p ==1
                        k_index_group_matrix(1,1) = j;  %1st hit on t channel doesn't have a k index
                    else
                        k_index_group_matrix(p,1) = k_index_group_array(p); %record k index for removal later
                    end
                else
                    reconst_group_hits_location_matrix(p,hit_group_to_reconst(p,1)+2) = 1;    %other channels are easy, as 0 -> x1 (column 2) etc
                    reconst_group_hits_time_matrix(p,hit_group_to_reconst(p,1)+2) = hit_group_to_reconst(p,2);
                    k_index_group_matrix(p,hit_group_to_reconst(p,1)+2) = k_index_group_array(p);
                end

            end

            %reconst_group_hits_location_matrix;
            %reconst_group_hits_time_matrix
            %k_index_group_matrix
            %dld_output_sorted(k_index_group_matrix(1,1),2)

            reconst_group_no_hits_matrix = sum(reconst_group_hits_location_matrix,1);   %1D array with number of hits on each channel in group

            reconst_group_hits_time_matrix_sorted = sort(reconst_group_hits_time_matrix,'descend');    %puts zeros at bottom

            if (reconst_group_no_hits_matrix(2)== 0 && reconst_group_no_hits_matrix(3) == 0) || (reconst_group_no_hits_matrix(4)== 0 && reconst_group_no_hits_matrix(5) == 0)

                continue    %if there isn't at least 1 hit on an x and y channel, then we can't reconstruct the group
            end

            %Now we have the sorted time array, plus know number of hits on each channel in the group to reconstruct

            for t_count = 1:reconst_group_no_hits_matrix(1) %loop over number of hits on t channel.  I could probably ignore this, since we check every t anyway in j...
                if reconst_group_no_hits_matrix(2) && reconst_group_no_hits_matrix(3)   %if x1 and x2 present, try to match them up, and then check for y's
                    good_hit_2x = 0;
                    t_mcp_hit = reconst_group_hits_time_matrix_sorted(t_count,1);
                    for x1_count = 1:reconst_group_no_hits_matrix(2)   %loop over x1 hits
                        for x2_count = 1:reconst_group_no_hits_matrix(3)

                            check1 = reconst_group_hits_time_matrix_sorted(x1_count,2)+reconst_group_hits_time_matrix_sorted(x2_count,3)-2*(t_mcp_hit-mean_t_diff)-tsum; %x1+x2-2t = tsum within tolerance
                            if abs(check1) < tolerance
                                x1_hit_reconst = reconst_group_hits_time_matrix_sorted(x1_count,2);   %assign x1
                                x2_hit_reconst = reconst_group_hits_time_matrix_sorted(x2_count,3);   %assign x2
                                good_hit_2x = good_hit_2x + 1;    %will become >1 if multi events match
                                x1_count_good = reconst_group_no_hits_matrix(2) - x1_count+1;   %tag the good count (if multi events matching then it will be overwritten)
                                x2_count_good = reconst_group_no_hits_matrix(3)-x2_count+1;    %note that since we re-ordered the matrix in reverse order, have to do the reverse for the indices
                            end

                        end
                    end

                    if good_hit_2x == 1   %only 1 x1/2 pair found
                        if reconst_group_no_hits_matrix(4) > 0   %if there's at least 1 trigger on y1 try to match it up with our x1/2 pair
                            good_hit_2x_y1 = 0;
                            for y1_count = 1:reconst_group_no_hits_matrix(4)
                                check2 = x1_hit_reconst + x2_hit_reconst - reconst_group_hits_time_matrix_sorted(y1_count,4) - (t_mcp_hit-mean_t_diff); %x1+x2 -y1 -2(t-tdiff) -tsum
                                if abs(check2) < max_group_time
                                    y1_hit_reconst = reconst_group_hits_time_matrix_sorted(y1_count,4);  %store y1 if good
                                    %tsum_y_correction = 0;
                                    x_val_on_vt = (x1_hit_reconst-x2_hit_reconst);
                                    tsum_y_correction = a_y*x_val_on_vt^3+b_y*x_val_on_vt^2+c_y*x_val_on_vt+d_y;    %find correcting factor for tsum
                                    y2_hit_reconst = 2*(t_mcp_hit-mean_t_diff)+tsum + tsum_y_correction - y1_hit_reconst; %reconstruct y2 to store. Need to use tsum array eventually
                                    y1_count_good = reconst_group_no_hits_matrix(4) - y1_count+1;   %tag the good count
                                    good_hit_2x_y1 = good_hit_2x_y1 + 1; %increment counter
                                end

                            end
                            if good_hit_2x_y1 == 1   %if only 1 unique y1 hit matches, then we've got a winner! Record the hit group
                                five_channel_events_mcp_reconst(event_counter_reconst_mcp,1) = t_mcp_hit;    %record time

                                five_channel_events_mcp_reconst(event_counter_reconst_mcp,2) = x1_hit_reconst;%record x1 and x2
                                five_channel_events_mcp_reconst(event_counter_reconst_mcp,3) = x2_hit_reconst;

                                five_channel_events_mcp_reconst(event_counter_reconst_mcp,4) = y1_hit_reconst;%record y1 and y2
                                five_channel_events_mcp_reconst(event_counter_reconst_mcp,5) = y2_hit_reconst;

                                event_counter_reconst_mcp = event_counter_reconst_mcp +1; %increment hit counter

                                t_count_good = t_count;   %tag the good count

                                t_hit_array_indices = find(reconst_group_hits_location_matrix(:,1));  %this should work, just need to convert to k index
                                marker_t = t_hit_array_indices(t_count_good);   %t event to remove from raw data matrix
                                raw_data_remove_t = k_index_group_matrix(marker_t,1);
                                dld_output_sorted(raw_data_remove_t,2) = NaN; %and now remove it

                                x1_hit_array_indices = find(reconst_group_hits_location_matrix(:,2));
                                marker_x1 = x1_hit_array_indices(x1_count_good);   %x1 event to remove from raw data matrix
                                raw_data_remove_x1 = k_index_group_matrix(marker_x1,2);
                                dld_output_sorted(raw_data_remove_x1,2) = NaN; %and now remove it

                                x2_hit_array_indices = find(reconst_group_hits_location_matrix(:,3));
                                marker_x2 = x2_hit_array_indices(x2_count_good);   %x2 event to remove from raw data matrix
                                raw_data_remove_x2 = k_index_group_matrix(marker_x2,3);
                                dld_output_sorted(raw_data_remove_x2,2) = NaN; %and now NaN it

                                y1_hit_array_indices = find(reconst_group_hits_location_matrix(:,4));
                                marker_y1 = y1_hit_array_indices(y1_count_good);   %y1 event to remove from raw data matrix
                                raw_data_remove_y1 = k_index_group_matrix(marker_y1,4);
                                dld_output_sorted(raw_data_remove_y1,2) = NaN; %and now NaN it

                                %                             y2_hit_array_indices = find(reconst_group_hits_location_matrix(:,5));
                                %                             marker_y2 = y2_hit_array_indices(y2_count_good);   %y2 event to remove from raw data matrix
                                %                             raw_data_remove_y2 = k_index_group_matrix(marker_y2);
                                %                             dld_output_sorted(raw_data_remove_y2,2) = NaN; %and now NaN it

                                continue %and since we've matched that mcp trigger, move along to the next one

                            end
                            %I think there still could be a hit matching this x pair, as there may be y2 and y1. Actually, can just remove else below I think...

                        end    %was else. now we know that there has to be at least one trigger on y2
                        good_hit_2x_y2 = 0;
                        for y2_count = 1:reconst_group_no_hits_matrix(5)

                            check3 = x1_hit_reconst + x2_hit_reconst - reconst_group_hits_time_matrix_sorted(y2_count,5) - (t_mcp_hit-mean_t_diff); %x1+x2 -y2 -t
                            if abs(check3) < max_group_time
                                y2_hit_reconst = reconst_group_hits_time_matrix_sorted(y2_count,5);  %store y2 if good
                                %tsum_y_correction = 0;
                                x_val_on_vt = (x1_hit_reconst-x2_hit_reconst)*v_perp_x*bin_time;
                                tsum_y_correction = a_y*x_val_on_vt^3+b_y*x_val_on_vt^2+c_y*x_val_on_vt+d_y;
                                y1_hit_reconst = 2*(t_mcp_hit-mean_t_diff)+tsum + tsum_y_correction - y2_hit_reconst; %reconstruct y1 to store
                                y2_count_good = reconst_group_no_hits_matrix(5) - y2_count+1;   %tag the good count
                                good_hit_2x_y2 = good_hit_2x_y2 + 1; %increment counter
                            end

                        end
                        if good_hit_2x_y2 == 1   %if only 1 unique y2 hit matches, then we've got a winner! Record the hit group
                            five_channel_events_mcp_reconst(event_counter_reconst_mcp,1) = t_mcp_hit;    %record time

                            five_channel_events_mcp_reconst(event_counter_reconst_mcp,2) = x1_hit_reconst;%record x1 and x2
                            five_channel_events_mcp_reconst(event_counter_reconst_mcp,3) = x2_hit_reconst;

                            five_channel_events_mcp_reconst(event_counter_reconst_mcp,4) = y1_hit_reconst;%record y1 and y2
                            five_channel_events_mcp_reconst(event_counter_reconst_mcp,5) = y2_hit_reconst;

                            event_counter_reconst_mcp = event_counter_reconst_mcp +1; %increment hit counter

                            t_count_good = t_count;   %tag the good count

                            t_hit_array_indices = find(reconst_group_hits_location_matrix(:,1));  %this should work, just need to convert to k index
                            marker_t = t_hit_array_indices(t_count_good);   %t event to remove from raw data matrix
                            raw_data_remove_t = k_index_group_matrix(marker_t,1);
                            dld_output_sorted(raw_data_remove_t,2) = NaN; %and now remove it

                            x1_hit_array_indices = find(reconst_group_hits_location_matrix(:,2));
                            marker_x1 = x1_hit_array_indices(x1_count_good);   %x1 event to remove from raw data matrix
                            raw_data_remove_x1 = k_index_group_matrix(marker_x1,2);
                            dld_output_sorted(raw_data_remove_x1,2) = NaN; %and now remove it

                            x2_hit_array_indices = find(reconst_group_hits_location_matrix(:,3));
                            marker_x2 = x2_hit_array_indices(x2_count_good);  %x2 event to remove from raw data matrix
                            raw_data_remove_x2 = k_index_group_matrix(marker_x2,3);
                            dld_output_sorted(raw_data_remove_x2,2) = NaN; %and now NaN it

                            %                             y1_hit_array_indices = find(reconst_group_hits_location_matrix(:,4));  %don't need to NaN this, since we reconstructed it
                            %                             marker_y1 = y1_hit_array_indices(y1_count_good);   %y1 event to remove from raw data matrix
                            %                             raw_data_remove_y1 = k_index_group_matrix(marker_y1,4);
                            %                             dld_output_sorted(raw_data_remove_y1,2) = NaN; %and now NaN it

                            y2_hit_array_indices = find(reconst_group_hits_location_matrix(:,5));
                            marker_y2 = y2_hit_array_indices(y2_count_good);   %y2 event to remove from raw data matrix
                            raw_data_remove_y2 = k_index_group_matrix(marker_y2,5);
                            dld_output_sorted(raw_data_remove_y2,2) = NaN; %and now NaN it

                            continue %and since we've matched that mcp trigger, move along to the next one
                        end
                        %deleted the end on left, this line

                    end




                end

                %if 1 x is missing, see if both y1 and y2 are present

                if reconst_group_no_hits_matrix(4) && reconst_group_no_hits_matrix(5)
                    good_hit_2y = 0;
                    t_mcp_hit = reconst_group_hits_time_matrix_sorted(t_count,1);
                    for y1_count = 1:reconst_group_no_hits_matrix(4)   %loop over x1 hits
                        for y2_count = 1:reconst_group_no_hits_matrix(5)

                            check1 = reconst_group_hits_time_matrix_sorted(y1_count,4)+reconst_group_hits_time_matrix_sorted(y2_count,5)-2*(t_mcp_hit-mean_t_diff)-tsum; %y1+y2-2t = tsum within tolerance
                            if abs(check1) < tolerance
                                y1_hit_reconst = reconst_group_hits_time_matrix_sorted(y1_count,4);   %assign y1
                                y2_hit_reconst = reconst_group_hits_time_matrix_sorted(y2_count,5);   %assign y2
                                good_hit_2y = good_hit_2y + 1;    %will become >1 if multi events match
                                y1_count_good = reconst_group_no_hits_matrix(4) - y1_count+1;   %tag the good count (if multi events matching then it will be overwritten)
                                y2_count_good = reconst_group_no_hits_matrix(5)-y2_count+1;    %note that since we re-ordered the matrix in reverse order, have to do the reverse for the indices
                            end

                        end
                    end

                    if good_hit_2y == 1   %only 1 y1/2 pair found
                        if reconst_group_no_hits_matrix(2) > 0   %if there's at least 1 trigger on x1 try to match it up with our y1/2 pair
                            good_hit_2y_x1 = 0;
                            for x1_count = 1:reconst_group_no_hits_matrix(2)
                                check2 = y1_hit_reconst + y2_hit_reconst - reconst_group_hits_time_matrix_sorted(x1_count,2) - (t_mcp_hit-mean_t_diff); %y1+y2 -x1 -2(t-tdiff) -tsum
                                if abs(check2) < max_group_time
                                    x1_hit_reconst = reconst_group_hits_time_matrix_sorted(x1_count,2);  %store x1 if good
                                    %tsum_x_correction = 0;
                                    y_val = (y1_hit_reconst-y2_hit_reconst)*v_perp_y*bin_time;
                                    tsum_x_correction = a_x*y_val^3+b_x*y_val^2+c_x*y_val+d_x;
                                    x2_hit_reconst = 2*(t_mcp_hit-mean_t_diff)+tsum + tsum_x_correction - x1_hit_reconst; %reconstruct x2 to store. Need to use tsum array eventually
                                    x1_count_good = reconst_group_no_hits_matrix(2) - x1_count+1;   %tag the good count
                                    good_hit_2y_x1 = good_hit_2y_x1 + 1; %increment counter
                                end

                            end
                            if good_hit_2y_x1 == 1   %if only 1 unique y1 hit matches, then we've got a winner! Record the hit group
                                five_channel_events_mcp_reconst(event_counter_reconst_mcp,1) = t_mcp_hit;    %record time

                                five_channel_events_mcp_reconst(event_counter_reconst_mcp,2) = x1_hit_reconst;%record x1 and x2
                                five_channel_events_mcp_reconst(event_counter_reconst_mcp,3) = x2_hit_reconst;

                                five_channel_events_mcp_reconst(event_counter_reconst_mcp,4) = y1_hit_reconst;%record y1 and y2
                                five_channel_events_mcp_reconst(event_counter_reconst_mcp,5) = y2_hit_reconst;

                                event_counter_reconst_mcp = event_counter_reconst_mcp +1; %increment hit counter

                                t_count_good = t_count;   %tag the good count

                                t_hit_array_indices = find(reconst_group_hits_location_matrix(:,1));  %this should work, just need to convert to k index
                                marker_t = t_hit_array_indices(t_count_good);   %t event to remove from raw data matrix
                                raw_data_remove_t = k_index_group_matrix(marker_t,1);
                                dld_output_sorted(raw_data_remove_t,2) = NaN; %and now remove it

                                x1_hit_array_indices = find(reconst_group_hits_location_matrix(:,2));
                                marker_x1 = x1_hit_array_indices(x1_count_good);   %x1 event to remove from raw data matrix
                                raw_data_remove_x1 = k_index_group_matrix(marker_x1,2);
                                dld_output_sorted(raw_data_remove_x1,2) = NaN; %and now remove it

                                %                                 x2_hit_array_indices = find(reconst_group_hits_location_matrix(:,3));
                                %                                 marker_x2 = x2_hit_array_indices(x2_count_good);   %x2 event to remove from raw data matrix
                                %                                 raw_data_remove_x2 = k_index_group_matrix(marker_x2);
                                %                                 dld_output_sorted(raw_data_remove_x2,2) = NaN; %and now NaN it

                                y1_hit_array_indices = find(reconst_group_hits_location_matrix(:,4));
                                marker_y1 = y1_hit_array_indices(y1_count_good);   %y1 event to remove from raw data matrix
                                raw_data_remove_y1 = k_index_group_matrix(marker_y1,4);
                                dld_output_sorted(raw_data_remove_y1,2) = NaN; %and now NaN it

                                y2_hit_array_indices = find(reconst_group_hits_location_matrix(:,5));
                                marker_y2 = y2_hit_array_indices(y2_count_good);   %y2 event to remove from raw data matrix
                                raw_data_remove_y2 = k_index_group_matrix(marker_y2,5);
                                dld_output_sorted(raw_data_remove_y2,2) = NaN; %and now NaN it

                                continue %and since we've matched that mcp trigger, move along to the next one
                            end
                        end    %was else. now we know that there has to be at least a trigger on x2
                        good_hit_2y_x2 = 0;
                        for x2_count = 1:reconst_group_no_hits_matrix(3)

                            check3 = y1_hit_reconst + y2_hit_reconst - reconst_group_hits_time_matrix_sorted(x2_count,3) - (t_mcp_hit-mean_t_diff); %y1+y2 -x2 -t
                            if abs(check3) < max_group_time
                                x2_hit_reconst = reconst_group_hits_time_matrix_sorted(x2_count,3);  %store y2 if good
                                y_val = (y1_hit_reconst-y2_hit_reconst)*v_perp_y*bin_time;
                                tsum_x_correction = a_x*y_val^3+b_x*y_val^2+c_x*y_val+d_x;
                                x1_hit_reconst = 2*(t_mcp_hit-mean_t_diff)+tsum + tsum_x_correction - x2_hit_reconst; %reconstruct y1 to store
                                x2_count_good = reconst_group_no_hits_matrix(3) - x2_count+1;   %tag the good count
                                good_hit_2y_x2 = good_hit_2y_x2 + 1; %increment counter
                            end

                        end
                        if good_hit_2y_x2 == 1   %if only 1 unique y1 hit matches, then we've got a winner! Record the hit group
                            five_channel_events_mcp_reconst(event_counter_reconst_mcp,1) = t_mcp_hit;    %record time

                            five_channel_events_mcp_reconst(event_counter_reconst_mcp,2) = x1_hit_reconst;%record x1 and x2
                            five_channel_events_mcp_reconst(event_counter_reconst_mcp,3) = x2_hit_reconst;

                            five_channel_events_mcp_reconst(event_counter_reconst_mcp,4) = y1_hit_reconst;%record y1 and y2
                            five_channel_events_mcp_reconst(event_counter_reconst_mcp,5) = y2_hit_reconst;

                            event_counter_reconst_mcp = event_counter_reconst_mcp +1; %increment hit counter

                            t_count_good = t_count;   %tag the good count

                            t_hit_array_indices = find(reconst_group_hits_location_matrix(:,1));  %this should work, just need to convert to k index
                            marker_t = t_hit_array_indices(t_count_good);   %t event to remove from raw data matrix
                            raw_data_remove_t = k_index_group_matrix(marker_t,1);
                            dld_output_sorted(raw_data_remove_t,2) = NaN; %and now remove it

                            %                                 x1_hit_array_indices = find(reconst_group_hits_location_matrix(:,2));
                            %                                 marker_x1 = x1_hit_array_indices(x1_count_good);   %x1 event to remove from raw data matrix
                            %                                 raw_data_remove_x1 = k_index_group_matrix(marker_x1);
                            %                                 dld_output_sorted(raw_data_remove_x1,2) = NaN; %and now remove it

                            x2_hit_array_indices = find(reconst_group_hits_location_matrix(:,3));
                            marker_x2 = x2_hit_array_indices(x2_count_good);   %x2 event to remove from raw data matrix
                            raw_data_remove_x2 = k_index_group_matrix(marker_x2,3);
                            dld_output_sorted(raw_data_remove_x2,2) = NaN; %and now NaN it

                            y1_hit_array_indices = find(reconst_group_hits_location_matrix(:,4));
                            marker_y1 = y1_hit_array_indices(y1_count_good);   %y1 event to remove from raw data matrix
                            raw_data_remove_y1 = k_index_group_matrix(marker_y1,4);
                            dld_output_sorted(raw_data_remove_y1,2) = NaN; %and now NaN it

                            y2_hit_array_indices = find(reconst_group_hits_location_matrix(:,5));
                            marker_y2 = y2_hit_array_indices(y2_count_good);   %y2 event to remove from raw data matrix
                            raw_data_remove_y2 = k_index_group_matrix(marker_y2,5);
                            dld_output_sorted(raw_data_remove_y2,2) = NaN; %and now NaN it

                            continue %and since we've matched that mcp trigger, move along to the next one
                        end
                        %deleted end on this line

                    end


                end
                %if both y and x are missing a channel, then have to work harder to reconstruct event.
                %end

                %         good_hit_2channels = good_hit_2x_y2 + good_hit_2x_y1 + good_hit_2y_x1 + good_hit_2y_x2;
                %
                %         if good_hit_2channels > 0
                %             continue
                %         end
                %
                %continue

                good_hit_counter_x_y = 0;  % We only need to worry about x/y, not channels individually, as we only want a unique hit (if there's more than one of x or y channel for this mcp pulse we should have found it already)


                if reconst_group_no_hits_matrix(2) > 0  %if there are hits on x1, try to match them with a y
                    
                    t_mcp_hit = reconst_group_hits_time_matrix_sorted(t_count,1);
                    %dummy = 42

                    for x1_count = 1:reconst_group_no_hits_matrix(2)

                        %good_hit_counter_x1_y2 = 0;

                        check1 = reconst_group_hits_time_matrix_sorted(x1_count,2)-(t_mcp_hit-mean_t_diff);  %x1 - t
                        if abs(check1) < max_group_time    %check if x1 is within group time of t

                            if reconst_group_no_hits_matrix(4) > 0  % if hits on y1, try to match them with x1
                                for y1_count = 1:reconst_group_no_hits_matrix(4)

                                    check2 = reconst_group_hits_time_matrix_sorted(y1_count,4)-(t_mcp_hit-mean_t_diff);  %y1 - t
                                    if abs(check2) < max_group_time    %check if y1 is within group time of t.  If so we may have a hit to reconstruct
                                        x1_hit_reconst = reconst_group_hits_time_matrix_sorted(x1_count,2);  %store x1
                                        x2_hit_reconst = 2*(t_mcp_hit-mean_t_diff)+tsum - x1_hit_reconst; %reconstruct x2 to store
                                        y1_hit_reconst = reconst_group_hits_time_matrix_sorted(y1_count,4);  %store y1 if good
                                        y2_hit_reconst = 2*(t_mcp_hit-mean_t_diff)+tsum - y1_hit_reconst; %reconstruct y2 to store
                                        good_hit_counter_x_y = good_hit_counter_x_y + 1;
                                        x1_count_good = reconst_group_no_hits_matrix(2) - x1_count+1;   %tag the good count
                                        y1_count_good = reconst_group_no_hits_matrix(4) - y1_count+1;   %tag the good count
                                        x_channel = 1;
                                        y_channel = 1;
                                    end

                                end
                            end

                            if reconst_group_no_hits_matrix(5) > 0  % if hits on y2, try to match them with x1

                                for y2_count = 1:reconst_group_no_hits_matrix(5)

                                    check3 = reconst_group_hits_time_matrix_sorted(y1_count,4)-(t_mcp_hit-mean_t_diff);  %y1 - t
                                    if abs(check3) < max_group_time    %check if y1 is within group time of t.  If so we may have a hit to reconstruct
                                        x1_hit_reconst = reconst_group_hits_time_matrix_sorted(x1_count,2);  %store x1
                                        x2_hit_reconst = 2*(t_mcp_hit-mean_t_diff)+tsum - x1_hit_reconst; %reconstruct x2 to store
                                        y2_hit_reconst = reconst_group_hits_time_matrix_sorted(y2_count,5);  %store y2 if good
                                        y1_hit_reconst = 2*(t_mcp_hit-mean_t_diff)+tsum - y2_hit_reconst; %reconstruct y1 to store
                                        good_hit_counter_x_y = good_hit_counter_x_y + 1;
                                        x1_count_good = reconst_group_no_hits_matrix(2) - x1_count+1;   %tag the good count
                                        y2_count_good = reconst_group_no_hits_matrix(5) - y2_count+1;   %tag the good count
                                        x_channel = 1;
                                        y_channel = 2;
                                    end

                                end

                            end

                        end

                    end

                end

                if good_hit_counter_x_y > 1 %if multi hit we're unable to distinguish them, so move along to next mcp trigger
                    continue
                end

                if reconst_group_no_hits_matrix(3) > 0  %if there are hits on x2, try to match them with a y.  The if statement is necessary as there could be a failed hit on both of x1 + yn and x2 + yn within the same group
                    
                    t_mcp_hit = reconst_group_hits_time_matrix_sorted(t_count,1);

                    for x2_count = 1:reconst_group_no_hits_matrix(3)

                        check1 = reconst_group_hits_time_matrix_sorted(x2_count,3)-(t_mcp_hit-mean_t_diff);  %x2 - t
                        if abs(check1) < max_group_time    %check if x1 is within group time of t

                            if reconst_group_no_hits_matrix(4) > 0  % if hits on y1, try to match them with x2
                                for y1_count = 1:reconst_group_no_hits_matrix(4)

                                    check5 = reconst_group_hits_time_matrix_sorted(y1_count,4)-(t_mcp_hit-mean_t_diff);  %y1 - t
                                    if abs(check5) < max_group_time    %check if y1 is within group time of t.  If so we may have a hit to reconstruct
                                        x2_hit_reconst = reconst_group_hits_time_matrix_sorted(x2_count,3);  %store x2
                                        x1_hit_reconst = 2*(t_mcp_hit-mean_t_diff)+tsum - x2_hit_reconst; %reconstruct x1 to store
                                        y1_hit_reconst = reconst_group_hits_time_matrix_sorted(y1_count,4);  %store y1 if good
                                        y2_hit_reconst = 2*(t_mcp_hit-mean_t_diff)+tsum - y1_hit_reconst; %reconstruct y2 to store
                                        good_hit_counter_x_y = good_hit_counter_x_y + 1;
                                        x2_count_good = reconst_group_no_hits_matrix(3) - x2_count+1;   %tag the good count
                                        y1_count_good = reconst_group_no_hits_matrix(4) - y1_count+1;   %tag the good count
                                        x_channel = 2;
                                        y_channel = 1;
                                    end

                                end
                            end

                            if reconst_group_no_hits_matrix(5) > 0  % if hits on y2, try to match them with x2
                                for y2_count = 1:reconst_group_no_hits_matrix(5)

                                    check6 = reconst_group_hits_time_matrix_sorted(y2_count,5)-(t_mcp_hit-mean_t_diff);  %y2 - t
                                    if abs(check6) < max_group_time    %check if y2 is within group time of t.  If so we may have a hit to reconstruct
                                        x2_hit_reconst = reconst_group_hits_time_matrix_sorted(x2_count,3);  %store x2
                                        x1_hit_reconst = 2*(t_mcp_hit-mean_t_diff)+tsum - x2_hit_reconst; %reconstruct x1 to store
                                        y2_hit_reconst = reconst_group_hits_time_matrix_sorted(y2_count,5);  %store y2 if good
                                        y1_hit_reconst = 2*(t_mcp_hit-mean_t_diff)+tsum - y2_hit_reconst; %reconstruct y1 to store
                                        good_hit_counter_x_y = good_hit_counter_x_y + 1;
                                        x2_count_good = reconst_group_no_hits_matrix(3) - x2_count+1;   %tag the good count
                                        y2_count_good = reconst_group_no_hits_matrix(5) - y2_count+1;   %tag the good count
                                        x_channel = 2;
                                        y_channel = 2;
                                    end

                                end

                            end

                        end

                    end

                end

                if good_hit_counter_x_y == 1 %if we have only 1 unique x1/2 and 1 unique y1/2, we can store it.  Otherwise move along to next mcp trigger

                    five_channel_events_mcp_reconst(event_counter_reconst_mcp,1) = t_mcp_hit;    %record time

                    five_channel_events_mcp_reconst(event_counter_reconst_mcp,2) = x1_hit_reconst;%record x1 and x2
                    five_channel_events_mcp_reconst(event_counter_reconst_mcp,3) = x2_hit_reconst;

                    five_channel_events_mcp_reconst(event_counter_reconst_mcp,4) = y1_hit_reconst;%record y1 and y2
                    five_channel_events_mcp_reconst(event_counter_reconst_mcp,5) = y2_hit_reconst;

                    event_counter_reconst_mcp = event_counter_reconst_mcp +1; %increment hit counter

                    t_count_good = t_count;   %tag the good count

                    t_hit_array_indices = find(reconst_group_hits_location_matrix(:,1));  %this should work, just need to convert to k index
                    marker_t = t_hit_array_indices(t_count_good);   %t event to remove from raw data matrix
                    raw_data_remove_t = k_index_group_matrix(marker_t,1);
                    dld_output_sorted(raw_data_remove_t,2) = NaN; %and now remove it

                    if x_channel == 1   %only remove the x channel we didn't reconstruct from the datastream
                        x1_hit_array_indices = find(reconst_group_hits_location_matrix(:,2));
                        marker_x1 = x1_hit_array_indices(x1_count_good);   %x1 event to remove from raw data matrix
                        raw_data_remove_x1 = k_index_group_matrix(marker_x1,2);
                        dld_output_sorted(raw_data_remove_x1,2) = NaN; %and now remove it
                    else
                        x2_hit_array_indices = find(reconst_group_hits_location_matrix(:,3));
                        marker_x2 = x2_hit_array_indices(x2_count_good);   %x2 event to remove from raw data matrix
                        raw_data_remove_x2 = k_index_group_matrix(marker_x2,3);
                        dld_output_sorted(raw_data_remove_x2,2) = NaN; %and now NaN it
                    end

                    if y_channel ==1
                        y1_hit_array_indices = find(reconst_group_hits_location_matrix(:,4));
                        marker_y1 = y1_hit_array_indices(y1_count_good);   %y1 event to remove from raw data matrix
                        raw_data_remove_y1 = k_index_group_matrix(marker_y1,4);
                        dld_output_sorted(raw_data_remove_y1,2) = NaN; %and now NaN it
                    else
                        y2_hit_array_indices = find(reconst_group_hits_location_matrix(:,5));
                        marker_y2 = y2_hit_array_indices(y2_count_good);   %y2 event to remove from raw data matrix
                        raw_data_remove_y2 = k_index_group_matrix(marker_y2,5);
                        dld_output_sorted(raw_data_remove_y2,2) = NaN; %and now NaN it
                    end

                    continue

                else
                    continue    %And we're not going to any more effort to try to find an event, so may as well skip to next mcp trigger
                end

            end

        else
            continue
        end
    end

    good_five_channel_events_mcp_reconst =~ isnan(five_channel_events_mcp_reconst(:,1));    %clear pre-allocated NaN's
    five_channel_events_mcp_reconst = five_channel_events_mcp_reconst(good_five_channel_events_mcp_reconst,:);

    five_channel_events_mcp_reconst;

    dld_output_sorted_rowstokeep =~ isnan(dld_output_sorted(:,2));  %clear NaN rows from dld_output_sorted
    dld_output_sorted = dld_output_sorted(dld_output_sorted_rowstokeep,:);

    size_rem = size(dld_output_sorted);
    num_trigs_remaining = size_rem(1,1);

    howmanyrows_mcp_reconst=size(five_channel_events_mcp_reconst);               %% Number of successful counts
    three_channel_output_mcp_reconst=zeros(howmanyrows_mcp_reconst(1),3);               %% Initialise matrix
    Num_mcp_reconst_hits = howmanyrows_mcp_reconst(1);

    three_channel_output_mcp_reconst(:,1)=five_channel_events_mcp_reconst(:,1)*bin_time;
    three_channel_output_mcp_reconst(:,2)=(five_channel_events_mcp_reconst(:,2)-five_channel_events_mcp_reconst(:,3))*v_perp_x*bin_time;
    three_channel_output_mcp_reconst(:,3)=(five_channel_events_mcp_reconst(:,4)-five_channel_events_mcp_reconst(:,5))*v_perp_y*bin_time;

    three_channel_output_mcp_reconst;

%     toc

    dld_output_sorted;

else
    three_channel_output_mcp_reconst = [];
end


%%%%Now we have matched every mcp trigger where possible, reconstructing if necessary.  
%%%%All that remains is to try to match up events with 2x or 2y plus one other corner.

if reconst_3_corners_flag == 1    %again, check flag to see if we're reconstructing

    event_counter_nomcp_reconst_xpair = 1;
    event_counter_nomcp_reconst_ypair = 1;
    five_channel_output_nomcp_reconst = nan(ceil(num_trigs_remaining/3),5);
    which_row = 0;
    number_successes = 0;

    for count = 2:num_trigs_remaining   %loop over all remaining events
        if dld_output_sorted(count,1) == 0                              %% Pick out values of x1 (channel = 0) to match everything to
            x1_val = dld_output_sorted(count,2);                        %% This is x1 to match everything else to
            dummy = 0;                                                  %% Dummy variable to keep search local
            how_many_x = 0;

            for count2 = 1:search_no                                    %% Now looks at other detections near x1
                dummy = dummy+count2*(-1)^(count2);                     %% Go back 1 row, forward 2, back 3, to keep search efficient

                if (count+dummy) <= 0 ;                                  %% Quit if we go out of range
                    continue
                end
                if (count+dummy) > num_trigs_remaining
                    continue
                end
                if abs(dld_output_sorted(count+dummy,2) - x1_val) > T_sum;                                  %% Quit if we go out of range
                    continue
                end

                if dld_output_sorted(count+dummy,1) == 1                %% Pick out values of x2 (channel = 1)
                    x2_val = dld_output_sorted(count+dummy,2);

                    if abs(x1_val - x2_val) < (T_sum_tol_throw) %% See if x2 is close enough to check
                        how_many_x = how_many_x + 1;
                        if how_many_x > 1                               %% Breaks if multihit
                            break
                        else
                            x1 = x1_val;
                            x2 = x2_val;
                            which_row = which_row + 1;
                            tx = 0.5*(x1+x2-T_sum);                     %% Stores the value if everything goes correctly
                            %                         five_channel_output_nomcp_reconst(which_row,2) = x1 - tx;
                            %                         five_channel_output_nomcp_reconst(which_row,3) = x2 - tx;

                            %x1_index = count;
                            x2_index = count + dummy;   %if not frozen now, then risk dummy moving on

                            %five_channel_output_4corners(which_row,1) = tx;    %can lead to problems if no y's to match the x. Plus is unnecessary
                        end
                    end
                end
            end

            if how_many_x == ~1     %% If cant find X or is not unique, go to next x1
                %blah = 42
                continue
            end

            dummy2 = 0;
            dummy3 = 0;
            how_many_y = 0;

            for count3 = 1:search_no                                    %% For this x pair try to match a y1
                dummy2 = dummy2+count3*(-1)^(count3);

                if (count+dummy2) <= 0
                    continue
                end
                if (count+dummy2) > num_trigs_remaining
                    continue
                end
                if abs(dld_output_sorted(count+dummy2,2) - tx) > max_group_time%x1_val) > T_sum;                                  %% Skip if we go out of range
                    continue
                end

                %got_to_here = 1
                %channnel = dld_output_sorted(count+dummy2,1)

                if dld_output_sorted(count+dummy2,1) == 2               %% Pick out values of y1 (channel = 2)
                    y1_val = dld_output_sorted(count+dummy2,2);         %% Since we've already done a range check, to have a y2 here means it must be good

                    how_many_y = how_many_y +1;
                    which_y_good = 1;

                    y1_index = count + dummy2;  %tag the event for later removal
                end
            end

            %got_to_here = 1

            for count4 = 1:search_no                                    %% For this x pair try to match a y2
                dummy3 = dummy3+count4*(-1)^(count4);
                %out_of_range_y = [0,0];

                if (count+dummy3) <= 0
                    continue
                end
                if (count+dummy3) > num_trigs_remaining
                    continue
                end
                if abs(dld_output_sorted(count+dummy3,2) - tx) > max_group_time %x1_val) > T_sum;                                  %% Skip if we go out of range
                    continue
                end

                %channel2 = dld_output_sorted(count+dummy3,1)

                if dld_output_sorted(count+dummy3,1) == 3               %% Pick out values of y2 (channel = 3)
                    y2_val = dld_output_sorted(count+dummy3,2);         %% Since we've already done a range check, to have a y2 here means it must be good

                    how_many_y = how_many_y +1;
                    which_y_good = 2;

                    y2_index = count + dummy3;  %tag the event for later removal
                end
            end

            if how_many_y == 1  %Only store the hit if unique
                five_channel_output_nomcp_reconst(which_row,1) = tx;    %store reconstructed time

                five_channel_output_nomcp_reconst(which_row,2) = x1 - tx;   %x1 and x2 from before
                five_channel_output_nomcp_reconst(which_row,3) = x2 - tx;

                dld_output_sorted(count,2) = NaN;   %x1_index = count
                dld_output_sorted(x2_index,2) = NaN;    %NaN the read data

                if which_y_good == 1    %if y1, reconstruct y2
                    five_channel_output_nomcp_reconst(which_row,4) = y1_val - tx;   %store y1 as it's good...
                    five_channel_output_nomcp_reconst(which_row,5) = tx-T_sum - y1_val; %and reconstruct y2

                    dld_output_sorted(y1_index,2) = NaN;    %NaN the read data

                    event_counter_nomcp_reconst_xpair = event_counter_nomcp_reconst_xpair + 1;  %increment counter
                else    %y2 must be good, so reconstruct y1
                    five_channel_output_nomcp_reconst(which_row,5) = y2_val - tx;   %store y2 as it's good...
                    five_channel_output_nomcp_reconst(which_row,4) = tx-T_sum - y2_val; %and reconstruct y1

                    dld_output_sorted(y2_index,2) = NaN;    %NaN the read data

                    event_counter_nomcp_reconst_xpair = event_counter_nomcp_reconst_xpair + 1;  %increment counter
                end
            end

        end
    end

    %%Now we've matched teh x pairs, let's check the y pairs

    which_row = 0;
    number_successes = 0;

    for count = 2:num_trigs_remaining   %loop over all remaining events
        if dld_output_sorted(count,1) == 2                              %% Pick out values of y1 (channel = 2) to match everything to
            y1_val = dld_output_sorted(count,2);                        %% This is y1 to match everything else to
            dummy = 0;                                                  %% Dummy variable to keep search local
            how_many_y = 0;

            for count2 = 1:search_no                                    %% Now looks at other detections near y1
                dummy = dummy+count2*(-1)^(count2);                     %% Go back 1 row, forward 2, back 3, to keep search efficient

                if (count+dummy) <= 0 ;                                  %% Quit if we go out of range
                    continue
                end
                if (count+dummy) > num_trigs_remaining
                    continue
                end
                if abs(dld_output_sorted(count+dummy,2) - y1_val) > T_sum;                                  %% Quit if we go out of range
                    continue
                end

                if dld_output_sorted(count+dummy,1) == 3                %% Pick out values of x2 (channel = 1)
                    y2_val = dld_output_sorted(count+dummy,2);

                    if abs(y1_val - y2_val) < (T_sum_tol_throw) %% See if x2 is close enough to check
                        how_many_y = how_many_y + 1;
                        if how_many_y > 1                               %% Breaks if multihit
                            break
                        else
                            y1 = y1_val;
                            y2 = y2_val;
                            which_row = which_row + 1;
                            ty = 0.5*(y1+y2-T_sum);                     %% Stores the value if everything goes correctly
                            %                         five_channel_output_nomcp_reconst(which_row,2) = x1 - tx;
                            %                         five_channel_output_nomcp_reconst(which_row,3) = x2 - tx;

                            %x1_index = count;
                            y2_index = count + dummy;   %if not frozen now, then risk dummy moving on

                            %five_channel_output_4corners(which_row,1) = tx;    %can lead to problems if no y's to match the x. Plus is unnecessary
                        end
                    end
                end
            end

            if how_many_y == ~1     %% If cant find Y or is not unique, go to next y1
                %blah = 42
                continue
            end

            dummy2 = 0;
            dummy3 = 0;
            how_many_x = 0;

            for count3 = 1:search_no                                    %% For this y pair try to match an x1
                dummy2 = dummy2+count3*(-1)^(count3);

                if (count+dummy2) <= 0
                    continue
                end
                if (count+dummy2) > num_trigs_remaining
                    continue
                end
                if abs(dld_output_sorted(count+dummy2,2) - ty) > max_group_time%x1_val) > T_sum;                                  %% Skip if we go out of range
                    continue
                end

                %got_to_here = 1
                %channnel = dld_output_sorted(count+dummy2,1)

                if dld_output_sorted(count+dummy2,1) == 0               %% Pick out values of x1 (channel = 0)
                    x1_val = dld_output_sorted(count+dummy2,2);         %% Since we've already done a range check, to have an x1 here means it must be good

                    how_many_x = how_many_x +1;
                    which_x_good = 1;

                    x1_index = count + dummy2;  %tag the event for later removal
                end
            end

            %got_to_here = 1

            for count4 = 1:search_no                                    %% For the y pair try to match an x2
                dummy3 = dummy3+count4*(-1)^(count4);
                %out_of_range_y = [0,0];

                if (count+dummy3) <= 0
                    continue
                end
                if (count+dummy3) > num_trigs_remaining
                    continue
                end
                if abs(dld_output_sorted(count+dummy3,2) - ty) > max_group_time %x1_val) > T_sum;                                  %% Skip if we go out of range
                    continue
                end

                %channel2 = dld_output_sorted(count+dummy3,1)

                if dld_output_sorted(count+dummy3,1) == 1               %% Pick out values of x2 (channel = 1)
                    x2_val = dld_output_sorted(count+dummy3,2);         %% Since we've already done a range check, to have an x2 here means it must be good

                    how_many_x = how_many_x +1;
                    which_x_good = 2;

                    x2_index = count + dummy2;  %tag the event for later removal
                end
            end

            if how_many_x == 1  %Only store the hit if unique
                five_channel_output_nomcp_reconst(which_row,1) = ty;    %store reconstructed time

                five_channel_output_nomcp_reconst(which_row,4) = y1 - ty;   %y1 and y2 from before
                five_channel_output_nomcp_reconst(which_row,5) = y2 - ty;

                dld_output_sorted(count,2) = NaN;   %y1_index = count
                dld_output_sorted(y2_index,2) = NaN;    %NaN the read data

                if which_x_good == 1    %if x1, reconstruct x2
                    five_channel_output_nomcp_reconst(which_row,2) = x1_val - ty;   %store x1 as it's good...
                    five_channel_output_nomcp_reconst(which_row,3) = ty-T_sum - x1_val; %and reconstruct x2

                    dld_output_sorted(x1_index,2) = NaN;    %NaN the read data

                    event_counter_nomcp_reconst_ypair = event_counter_nomcp_reconst_ypair + 1;  %increment counter
                else    %x2 must be good, so reconstruct x1
                    five_channel_output_nomcp_reconst(which_row,3) = x2_val - ty;   %store x2 as it's good...
                    five_channel_output_nomcp_reconst(which_row,2) = ty-T_sum - x2_val; %and reconstruct yx

                    dld_output_sorted(x2_index,2) = NaN;    %NaN the read data

                    event_counter_nomcp_reconst_ypair = event_counter_nomcp_reconst_ypair + 1;  %increment counter
                end
            end

        end
    end

    event_counter_nomcp_reconst_xpair;
    event_counter_nomcp_reconst_ypair;


    good_five_channel_output_nomcp_reconst =~ isnan(five_channel_output_nomcp_reconst(:,1));    %clear pre-allocated NaN's
    five_channel_output_nomcp_reconst = five_channel_output_nomcp_reconst(good_five_channel_output_nomcp_reconst,:);

    five_channel_output_nomcp_reconst;

    dld_output_sorted_rowstokeep =~ isnan(dld_output_sorted(:,2));  %clear NaN rows from dld_output_sorted
    dld_output_sorted = dld_output_sorted(dld_output_sorted_rowstokeep,:);

    size_rem = size(dld_output_sorted);
    num_trigs_remaining = size_rem(1,1);

    howmanyrows_nomcp_reconst=size(five_channel_output_nomcp_reconst);               %% Number of successful counts
    three_channel_output_nomcp_reconst=zeros(howmanyrows_nomcp_reconst(1),3);               %% Initialise matrix
    Num_nomcp_reconst_hits = howmanyrows_nomcp_reconst(1);

    three_channel_output_nomcp_reconst(:,1)=five_channel_output_nomcp_reconst(:,1)*bin_time;
    three_channel_output_nomcp_reconst(:,2)=(five_channel_output_nomcp_reconst(:,2)-five_channel_output_nomcp_reconst(:,3))*v_perp_x*bin_time;
    three_channel_output_nomcp_reconst(:,3)=(five_channel_output_nomcp_reconst(:,4)-five_channel_output_nomcp_reconst(:,5))*v_perp_y*bin_time;

    %three_channel_output_nomcp_reconst;

else
    three_channel_output_nomcp_reconst = [];
end

three_channel_output_combined = cat(1,three_channel_output_5trigs,three_channel_output_5trigs_multi,three_channel_output_mcp_reconst,three_channel_output_4corners,three_channel_output_nomcp_reconst);


[~, NewRowNumber3chan] = sort(three_channel_output_combined(:,1));

three_channel_output_combined_sorted = three_channel_output_combined(NewRowNumber3chan,:);


%toc

% good_three_channel_output_reconst_mcp = three_channel_output_reconst_mcp(:,1)~=0;
% three_channel_output_reconst_mcp = three_channel_output_reconst_mcp(good_three_channel_output_reconst_mcp,:) %  clear extra pre-allocated zeros in matrix
%





%%%%%%%%Code for calculating tsummap. Ignore normally.

% bin_size_tsummap = 2;   %in mm
% bin_size_tsummap_norm = 80/bin_size_tsummap;
%
% Matrix_tsumx = zeros(bin_size_tsummap_norm,bin_size_tsummap_norm);
% Matrix_tsum_count = zeros(bin_size_tsummap_norm,bin_size_tsummap_norm);   %since both will have same number of hits
% Matrix_tsumy = zeros(bin_size_tsummap_norm,bin_size_tsummap_norm);
%
%
% for j = 1:event_counter_5trigs
%     if abs((three_channel_output_5trigs(j,2))) >.04
%         continue
%     end
%     if abs((three_channel_output_5trigs(j,3))) >.04
%         continue
%     end
%
%     x_index = ceil((three_channel_output_5trigs(j,2)+.04)*1000/bin_size_tsummap);   %convert x pos to matrix index
%     y_index = ceil((three_channel_output_5trigs(j,3)+.04)*1000/bin_size_tsummap);
%     Matrix_tsumx(x_index,y_index) = Matrix_tsumx(x_index,y_index) + three_channel_output_5trigs(j,4);  %add on tsumx to correct location in map matrix
%     Matrix_tsum_count(x_index,y_index) = 1 + Matrix_tsum_count(x_index,y_index);  %increment tsum count index
%     Matrix_tsumy(x_index,y_index) = Matrix_tsumy(x_index,y_index) + three_channel_output_5trigs(j,5);   %ditto for tsumy
%
% end
% Matrix_tsum_count;
%
%
%
% Matrix_tsumx_av = Matrix_tsumx./Matrix_tsum_count;
% Matrix_tsumy_av = Matrix_tsumy./Matrix_tsum_count;
%
%
% %figure1 = surf(1:80/bin_size_tsummap,1:80/bin_size_tsummap,Matrix_tsumx_av)
% figure2 = surf(1:80/bin_size_tsummap,1:80/bin_size_tsummap,Matrix_tsumy_av);
%
% % Matrix_tsumx_av_bad_cols = zeros(1,bin_size_tsummap_norm);
% % Matrix_tsumx_av_bad_rows = zeros(bin_size_tsummap_norm,1);
% % Matrix_tsumy_av_bad_cols = zeros(1,bin_size_tsummap_norm);
% % Matrix_tsumy_av_bad_rows = zeros(bin_size_tsummap_norm,1);
%
% % for r=1:bin_size_tsummap_norm
% %     for c=1:bin_size_tsummap_norm
% %         if Matrix_tsumx_av(r,c) == NaN
% %             Matrix_tsumx_av(r,c)=0;
% %             Matrix_tsumx_av_bad_rows(r,1) =  Matrix_tsumx_av_bad_rows(r,1)+1;
% %             Matrix_tsumx_av_bad_cols(c,1) =  Matrix_tsumx_av_bad_cols(1,c)+1;
% %         end
% %         if Matrix_tsumy_av(r,c) == NaN
% %             Matrix_tsumy_av(r,c)=0;
% %             Matrix_tsumy_av_bad_rows(r,1) =  Matrix_tsumy_av_bad_rows(r,1)+1;
% %             Matrix_tsumy_av_bad_cols(c,1) =  Matrix_tsumy_av_bad_cols(1,c)+1;
% %         end
% %     end
% % end
%
%
% bad_hits = isnan(Matrix_tsumx_av);
% tsumx_array_good = Matrix_tsumx_av;
% tsumy_array_good = Matrix_tsumy_av;
%
% Matrix_tsumx_av_bad_cols = sum(isnan(Matrix_tsumx_av),1);
% Matrix_tsumx_av_bad_rows = sum(isnan(Matrix_tsumx_av),2);
% Matrix_tsumy_av_bad_cols = sum(isnan(Matrix_tsumy_av),1);
% Matrix_tsumy_av_bad_rows = sum(isnan(Matrix_tsumy_av),2);
%
% for r=1:bin_size_tsummap_norm
%     for c=1:bin_size_tsummap_norm
%         if bad_hits(r,c) == 1
%             tsumx_array_good(r,c)=0;
%             tsumy_array_good(r,c)=0;
%
%         end
%
%     end
% end
%
% %Matrix_tsumx_av_bad_cols;
% %Matrix_tsumx_av_bad_rows;
% %blah = (bin_size_tsummap_norm)./(bin_size_tsummap_norm-Matrix_tsumx_av_bad_cols);
% %tsumx_array = zeros(80/bin_size_tsummap,1);
% %tsumx_array = mean(Matrix_tsumx_av,1)
% tsumx_array = mean(tsumx_array_good,2);
% tsumx_array = tsumx_array.*((bin_size_tsummap_norm)./(bin_size_tsummap_norm-Matrix_tsumx_av_bad_rows));
% %tsumx_array = mean(tsumx_array_good,1) %Wrong dimension; doesn't vary
% %tsumx_array = tsumx_array.*((bin_size_tsummap_norm)./(bin_size_tsummap_norm-Matrix_tsumx_av_bad_cols))
%
% tsumy_array = mean(tsumy_array_good,1);
% tsumy_array =tsumy_array.*((bin_size_tsummap_norm)./(bin_size_tsummap_norm-Matrix_tsumy_av_bad_cols));
% tsumy_array_output = tsumy_array.';
%
% event_counter_5trigs;
%
% [blah1, blah2] = min(three_channel_output_5trigs,[],1);
%
% [blah1, blah2] = max(three_channel_output_5trigs,[],1);

%%%%%END
end


%TDC_importer
%a super fast import function for the TDC data
%this does not handle the wrong format gracefully
%fairly sure the format is u(one decimal),u(64)
%its a bit messy to give the output as a 64bit matrix but previouse
%code is set up to theis type, also matlab wont speed up from u8 vs u64
function  data=tdc_importer(filepath)
fileID = fopen(filepath);
data = textscan(fileID,'%f64%f64',...%f32
    'Delimiter', ',', ...
    'MultipleDelimsAsOne',0,...
    'ReturnOnError', true,...
    'NumCharactersToSkip',0,... 
    'CollectOutput', true);
    %'EndOfLine','\n')
fclose(fileID);
data=data{1};
end