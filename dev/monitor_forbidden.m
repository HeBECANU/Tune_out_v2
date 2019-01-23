%monitor the forbidden transtion run
% monitor an experiment that uses measrument of the trap freq
% It would be usefull to get a decent idea of what the trap frequency is doing during a run without
% the need for a full processing of the data as in main
%  - processing each shot as it is made
%  - plot a history of the trap freq

% Known BUGS/ Possible Improvements
%
%
% Author: Bryce Henson
% email: Bryce.Henson@live.com
% Last revision:2018-10-01
% BEGIN USER VAR-------------------------------------------------
anal_opts.tdc_import.dir='\\amplpc29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\';
anal_opts.tdc_import.file_name='d';
anal_opts.tdc_import.force_load_save=false;   %takes precidence over force_reimport
anal_opts.tdc_import.force_reimport=true;
anal_opts.tdc_import.force_forc=false;
anal_opts.tdc_import.dld_xy_rot=0.61;

tmp_xlim=[-30e-3, 30e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
tmp_ylim=[-30e-3, 30e-3];
tlim=[0,4];
anal_opts.tdc_import.txylim=[tlim;tmp_xlim;tmp_ylim];


anal_opts.atom_laser.pulsedt=8.000e-3;
anal_opts.atom_laser.t0=0.41784; %center i ntime of the first pulse
anal_opts.atom_laser.start_pulse=1; %atom laser pulse to start with
anal_opts.atom_laser.pulses=325;
anal_opts.atom_laser.appr_osc_freq_guess=[52,48,40];
anal_opts.atom_laser.pulse_twindow=anal_opts.atom_laser.pulsedt*0.9;
anal_opts.atom_laser.xylim=anal_opts.tdc_import.txylim(2:3,:); %set same lims for pulses as import

anal_opts.global.fall_time=0.417;
anal_opts.global.qe=0.09;

anal_opts.trig_dld=20.3;
anal_opts.dld_aquire=4;
anal_opts.trig_ai_in=20;


anal_opts.history.shots=inf;

% END USER VAR-----------------------------------------------------------
fclose('all')
%add all subfolders
folder = fileparts(which(mfilename));
folder=strsplit(folder,filesep); %go up a directory
folder=strjoin(folder(1:end-1),filesep);
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

hebec_constants
anal_opts.tdc_import.mat_save=false;
anal_opts.global.velocity=const.g0*anal_opts.global.fall_time;

if anal_opts.tdc_import.dir(end) ~= '\', dirpath = [dirpath '\']; end
if (exist([anal_opts.tdc_import.dir,'out'], 'dir') == 0), mkdir([anal_opts.tdc_import.dir,'out']); end
 
anal_out.dir=sprintf('%sout\\monitor\\',...
    anal_opts.tdc_import.dir);
if (exist(anal_out.dir, 'dir') == 0), mkdir(anal_out.dir); end
anal_opts.global.out_dir=anal_out.dir;



%%
atom_loss_history=[];
atom_loss_history.frac_rem_log_val=[];
atom_loss_history.frac_rem_log_unc=[];
atom_loss_history.shot_num=[];
sfigure(1);
set(gcf,'color','w')
clf;
sfigure(2)
set(gcf,'color','w')
%%

loop_num=0;
while true
    pause(0.1)
    batch_data=[];
    anal_opts.tdc_import.shot_num=find_data_files(anal_opts.tdc_import);

    max_shot_num=max(anal_opts.tdc_import.shot_num);
    anal_opts.tdc_import.shot_num=anal_opts.tdc_import.shot_num(...
        anal_opts.tdc_import.shot_num>(max_shot_num-anal_opts.history.shots));
    %remove processed ones

    anal_opts.tdc_import.shot_num=anal_opts.tdc_import.shot_num(...
        ~ismember(anal_opts.tdc_import.shot_num, atom_loss_history.shot_num ) );

    if numel(anal_opts.tdc_import.shot_num)==0
            if mod(loop_num,4)==0
                pause(.2)
                fprintf('\b\b\b')
                loop_num=1;
            else
                pause(.1) %little wait animation
                fprintf('.')
                loop_num=loop_num+1;
            end
    else
        batch_data.mcp_tdc=import_mcp_tdc_data(anal_opts.tdc_import);
        %data.mcp_tdc=mcp_tdc_data;
        %just to give me a logical vector
        batch_data.mcp_tdc.all_ok=batch_data.mcp_tdc.num_counts>1e4;
        batch_data.mcp_tdc.all_ok(batch_data.mcp_tdc.all_ok)=...
            cellfun(@(x) x(end,1),batch_data.mcp_tdc.counts_txy(batch_data.mcp_tdc.all_ok))>anal_opts.dld_aquire*0.8;
        if sum(batch_data.mcp_tdc.all_ok)==0
            fprintf('waiting for file to be writen\n')
            pause(0.1)
        else
            batch_data.mcp_tdc.al_pulses=bin_al_pulses(anal_opts.atom_laser,batch_data);
            %%
            mask=batch_data.mcp_tdc.all_ok;
            
            %% FITTING THE ATOM NUMBER PART 1
            %use the inital few atom laser pulses in order to determine the atom number
            %not of that much benifit TBH
            anal_opts.atom_num_fit=[];
            anal_opts.atom_num_fit.pulses=[1,15]; %min,max index of pulses
            sfigure(2)
            subplot(2,1,1)
            %only show the fit if there is less than 5 shots
            anal_opts.atom_num_fit.plot.each_shot=sum(batch_data.mcp_tdc.all_ok)<5;

            anal_opts.atom_num_fit.plot.history=false;
            anal_opts.atom_num_fit.qe=anal_opts.global.qe;

            data.num_fit.pre_probe=fit_atom_number(anal_opts.atom_num_fit,batch_data);

            %% FITTING THE ATOM NUMBER PART 2
            sfigure(2)
            subplot(2,1,2)
            anal_opts.atom_num_fit.pulses=[144,210]; %min,max index of pulses
            data.num_fit.post_probe=fit_atom_number(anal_opts.atom_num_fit,batch_data);

            %% Compare

            atom_num_pre=cellfun(@(x) x(16),data.num_fit.pre_probe.fit_predict(mask));
            atom_num_post=cellfun(@(x) x(144),data.num_fit.post_probe.fit_predict(mask));
            atom_num_unc_pre=cellfun(@(x) x(16),data.num_fit.pre_probe.fit_predict_unc(mask));
            atom_num_unc_post=cellfun(@(x) x(144),data.num_fit.post_probe.fit_predict_unc(mask));
            
            frac_rem_vec_lin=atom_num_post./atom_num_pre;
            frac_rem_vec_log_val=log(atom_num_post./atom_num_pre);
            frac_rem_vec_log_unc=abs(sqrt((atom_num_unc_pre./atom_num_pre).^2+(atom_num_unc_post./atom_num_post).^2));
           
            
            atom_loss_history.shot_num=[atom_loss_history.shot_num,anal_opts.tdc_import.shot_num(mask)];
            atom_loss_history.frac_rem_log_val=[atom_loss_history.frac_rem_log_val,frac_rem_vec_log_val'];
            atom_loss_history.frac_rem_log_unc=[atom_loss_history.frac_rem_log_unc,frac_rem_vec_log_unc'];

            %trim the history vectors
            if numel(atom_loss_history.shot_num)>anal_opts.history.shots
                %bit sloppy but will assume they are the same length
                atom_loss_history.shot_num=atom_loss_history.shot_num(end-anal_opts.history.shots:end);
                atom_loss_history.frac_rem_log_val=atom_loss_history.frac_rem_log_val(end-anal_opts.history.shots:end);
                atom_loss_history.frac_rem_log_unc=atom_loss_history.frac_rem_log_unc(end-anal_opts.history.shots:end);
            end

            sfigure(1);
            errorbar(atom_loss_history.shot_num,...
                atom_loss_history.frac_rem_log_val,atom_loss_history.frac_rem_log_unc,...
                'kx-','MarkerSize',7,'CapSize',0,'LineWidth',1)
            xlabel('Shot Number')
            ylabel('Log Number remaining')
            pause(1e-6)

        end
    end

end
