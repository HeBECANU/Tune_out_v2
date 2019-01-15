

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%DOCUMENTATION
%bryce:what needs to be in the path for this to run properly
%explain briefly what this program does and what the user will need to edit

%to do
%Add conditional clauses to check for errors, write to log
%Add conditional to enable/disable settings updates - without having to restart LV? I guess one can
%just revert to v5 and reboot
%Move num_instr to instructions.m to put all customization in one file

%INPUT FROM LABVIEW
%mloop,boolean : do you want to interface with mloop or use matlab to scan over variables
%file,string : tells program where to look for setting files
%i ,integer:itteration number

%return
%exact_line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exact_line=''; %legacy to stop crashing

%setpoints = linspace(362848466.40,362871075.28,50);
%setpoints = linspace(362853060,362880320,50); 
power_setpt_list = linspace(8,0,50); %to ~ 362868100(40)


dir_this=fileparts(mfilename('fullpath'));   % get dir of this script
%path_log=strcat(dir_this,'\logs\','log_LabviewMatlab.txt');    % path of log file
path_log='Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\log_LabviewMatlab.txt';

%% Sequential setpoint
%avoid setting pointer to zero, i=iteration number passed from LabView
%add some random noise in the pointer to average out drifts a bit %round((rand(1)-0.5)*10
calibrate_interval=2;

if mod(i,calibrate_interval)==0
    %new_path='c:\remote\settings201819Nov193558.xml'; %probe off Lieutenant Angler
    new_path='c:\remote\settings201826Nov224020.xml'; %Lieutenant Fangtooth
    %write log entry
    f_log=fopen(path_log,'a');  % append to log-file
    nowdt=datetime('now');
    fprintf(f_log,'%.3f,%s,interfacev8,calibrate,itt,%u\n',posixtime(nowdt),datestr(nowdt,'yyyy-mm-ddTHH:MM:SS.FFF'),i);
    fclose(f_log);
else
    new_path='c:\remote\settings201826Nov224252.xml'; %Lieutenant Fangtooth
    pointer=floor(i/calibrate_interval)*(calibrate_interval-1)+ rem(i,calibrate_interval);
    pointer = mod(pointer-1,length(power_setpt_list))+1; 
    power_setpt = power_setpt_list(pointer);

    %write log entry
    f_log=fopen(path_log,'a');  % append to log-file
    nowdt=datetime('now');
    fprintf(f_log,'%.3f,%s,interfacev8,measure_probe,%f,itt,%u\n',posixtime(nowdt),datestr(nowdt,'yyyy-mm-ddTHH:MM:SS.FFF'),power_setpt,i);
    fclose(f_log);

    %List of all paths of control variables

    %path for ramp opt
    paths = {
        {{{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6733 M (Ch 0..7)',{'<Cluster>',33},'Final value, Amplitude (exp)',''},...
        {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6733 M (Ch 0..7)',{'<Cluster>',34},'Initial Value ',''},...
        {{'<Cluster>',1},{'<Cluster>',3},'Waveform Config Array For 6733 M (Ch 0..7)',{'<Cluster>',35},'Initial Value ',''}}};      % 14-param shunting profile

    % read in the xml settings file
    file = new_path;
    file_id = fopen(file,'r');
    j=0;
    while ischar(line) || j==0
        j = j+1;
        line = fgetl(file_id);
        set_file_lines{j} = line;
    end
    fclose(file_id);
    % make adjustments
    %load the mat file with predefine parameter sets
    param=[power_setpt];   %get the parameter set array
    param_limits=[[0,10];[0,10];[0,10]];

    conjoined_indx = [];

    for idx=1:numel(paths)
        %check if param values are within the hard coded safety range
        if param(idx)>param_limits(idx,2)
            param(idx) = param_limits(idx,2);
        elseif param(idx)<param_limits(idx,1)
            param(idx) = param_limits(idx,1);
        end

        %place the variable inputs into the paths to be changed
        if iscell(paths{idx}{end})
            % variables are bunched (i.e. conjoined)
            conjoined_indx = [conjoined_indx, idx];
            for j = 1:numel(paths{idx})
                paths{idx}{j}{end} = num2str(param(idx));
            end
        else
            % solitary variable
            paths{idx}{end} = num2str(param(idx));
        end
    end

    %split up conjoined variables
    if numel(conjoined_indx)>0
        temp = paths;
        paths = {};
        r=1;
        for j = 1:numel(temp)
            if any(j==conjoined_indx)
                for k = 1:numel(temp{j})
                    paths{r} = temp{j}{k};
                    r=r+1;
                end
            else
                paths{r} = temp{j};
                r=r+1;
            end
        end
    end

    read_paths={};

    %if we have multiple paths we can steup a looping structure as follows
    %paths = {all the different paths}
    %for j = 1:length(paths)
    %then replace path with paths{j}
    if true%i == 1
        exact_line =zeros(1,numel(paths));
        comment_line = zeros(1,numel(read_paths));
    end
    for k = 1:numel(paths)
        %bryce: what happens when the path string cant be found? is there an
        %error message?

        %select the path we want to change
        path = paths{k};

        %we only wish to find the exact path once, as it is computationally
        %expensive
        %if i==1
        if true %hack to see if finding each time helps
            comment_line = 1;
            for j = 1:length(path)-1
                count = 0;
                toggle = 1; 
                test = 0;
                if iscell(path{j})
                    name = path{j}{1};
                    value = path{j}{2};
                    end_name = strrep(name,'<','</');
                elseif ischar(path{j})
                    name = path{j};
                    value = 1;
                    end_name = '';
                end
                while count<value
                    comment_line = comment_line + 1;
                    if ~isempty(strfind(set_file_lines{1,comment_line},name))
                        if test == 0 && count>0
                            toggle = 1 - toggle;
                        end
                        if toggle
                            count = count + 1;
                            toggle = 1 - toggle;
                        end
                        test = test + 1;
                    elseif ~isempty(strfind(set_file_lines{1,comment_line},end_name))
                        test = test - 1;
                    end
                end
            end
            exact_line(1,k) = comment_line+1;
        end

        %set the value to what we want it to be
        start_ps=strfind(set_file_lines{exact_line(1,k)},'>');
        end_ps= strfind(set_file_lines{exact_line(1,k)},'<');
        strt = set_file_lines{exact_line(1,k)}(1:start_ps(1));
        fin = set_file_lines{exact_line(1,k)}(end_ps(2):end);
        set_file_lines{exact_line(1,k)} = strcat(strt,path{end},fin);
    end

    if true %just find each time
        for k = 1:numel(read_paths)
            %bryce: what happens when the path string cant be found? is there an
            %error message?
            %select the path we want to change
            path = paths{k};
            %we only wish to find the exact path once, as it is computationally
            %expensive

            comment_line = 1;
            for j = 1:length(path)-1
                count = 0;
                toggle = 1; 
                test = 0;
                if iscell(path{j})
                    name = path{j}{1};
                    value = path{j}{2};
                    end_name = strrep(name,'<','</');
                elseif ischar(path{j})
                    name = path{j};
                    value = 1;
                    end_name = '';
                end
                while count<value
                    comment_line = comment_line + 1;
                    if ~isempty(strfind(set_file_lines{1,comment_line},name))
                        if test == 0 && count>0
                            toggle = 1 - toggle;
                        end
                        if toggle
                            count = count + 1;
                            toggle = 1 - toggle;
                        end
                        test = test + 1;
                    elseif ~isempty(strfind(set_file_lines{1,comment_line},end_name))
                        test = test - 1;
                    end
                end
            end
            comment_line(1,k) = comment_line+1;

            %or use the exact line is already known
            % exact_line = {};

            %print out what the comment says
            j=0;
            end_of_comment = 1;
            while 0
                %display this line somehow
                %A{comment_line(1,k)+j}
                if isempty(strfind(set_file_lines{comment_line(1,k)+j},'</Val>'))
                    j = j+1;
                else
                    end_of_comment = 0; 
                end
            end
        end
    end
    %% write out new file
    new_path='c:\remote\temp.xml';
    file_id = fopen(new_path,'w');
    for j = 1:numel(set_file_lines)-2
            fprintf(file_id,'%s\n', set_file_lines{j});
    end
	fprintf(file_id,'%s', set_file_lines{end-1});
    fclose(file_id);

    %%clean up
    clear('A','set_file_lines');

    % f_log=fopen(path_log,'a');  % append to log-file
    % fprintf(f_log,[datestr(datetime,'yyyymmdd_HHMMSS'),' interfacev5    : finished. \n']);
    % fclose(f_log);

end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PARAMETER SCANNER HINTS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%Parameter scan
    % %for param scan only create an array of the values for each iteration

    % %Combining variables
    % %option one 
    % %product
    %a=[1 2]
    %b=[5 6]
    %transpose(combvec(a,b))
    %     1     5
    %     2     5
    %     1     6
    %     2     6
    % %which maps through all posible combinations

    % %Option two list
    % %steps though both in sequenc
    % %requires both lists to be same length !
    %transpose([a; b])
    %     1     5
    %     2     6

    %val1=linspace(0.6,3.4,30);
    %val2=linspace(0.6,3.4,30);
    %param_values=transpose([val1; val2]);
    %
    %param_values=param_values(randperm(size(param_values,1)),:);

    % END PARAMETER SCANNER HINTS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %END user settings
