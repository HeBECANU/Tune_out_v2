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
anal_opts.atom_laser.pulses=100;
anal_opts.atom_laser.appr_osc_freq_guess=[52,48,40];
anal_opts.atom_laser.pulse_twindow=anal_opts.atom_laser.pulsedt*0.9;
anal_opts.atom_laser.xylim=anal_opts.tdc_import.txylim(2:3,:); %set same lims for pulses as import

anal_opts.global.fall_time=0.417;
anal_opts.global.qe=0.09;

anal_opts.trig_dld=20.3;
anal_opts.dld_aquire=4;
anal_opts.trig_ai_in=20;


anal_opts.osc_fit.binsx=1000;
anal_opts.osc_fit.blur=1;
anal_opts.osc_fit.xlim=[-20,20]*1e-3;
anal_opts.osc_fit.tlim=[0.86,1.08];
anal_opts.osc_fit.dimesion=2; %Select coordinate to bin. 1=X, 2=Y.

anal_opts.history.shots=50;

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
trap_freq_history=[];
trap_freq_history.trap_freq_val=[];
trap_freq_history.trap_freq_unc=[];
trap_freq_history.shot_num=[];
sfigure(1);
set(gcf,'color','w')
clf;

%%

loop_num=0;
while true
    pause(0.1)
    batch_data=[];
    batch_data.shot_num=[];
    anal_opts.tdc_import.shot_num=find_data_files(anal_opts.tdc_import);

    max_shot_num=max(anal_opts.tdc_import.shot_num);
    anal_opts.tdc_import.shot_num=anal_opts.tdc_import.shot_num(...
        anal_opts.tdc_import.shot_num>(max_shot_num-anal_opts.history.shots));
    %remove processed ones

    anal_opts.tdc_import.shot_num=anal_opts.tdc_import.shot_num(...
        ~ismember(anal_opts.tdc_import.shot_num, trap_freq_history.shot_num ) );

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
        batch_data.mcp_tdc.all_ok=batch_data.mcp_tdc.num_counts>5e3;
        batch_data.mcp_tdc.all_ok(batch_data.mcp_tdc.all_ok)=...
            cellfun(@(x) x(end,1),batch_data.mcp_tdc.counts_txy(batch_data.mcp_tdc.all_ok))>anal_opts.dld_aquire*0.8;
        if sum(batch_data.mcp_tdc.all_ok)==0
            fprintf('waiting for file to be writen\n')
            pause(0.1)
        else
            batch_data.mcp_tdc.al_pulses=bin_al_pulses(anal_opts.atom_laser,batch_data);
            %%
            anal_opts.osc_fit.adaptive_freq=true; %estimate the starting trap freq 
            anal_opts.osc_fit.appr_osc_freq_guess=[52,44,40];
            anal_opts.osc_fit.freq_fit_tolerance=5;
            if sum(batch_data.mcp_tdc.all_ok)>2
                anal_opts.osc_fit.plot_fits=false;
            else
                anal_opts.osc_fit.plot_fits=true;
            end

            anal_opts.osc_fit.plot_err_history=false;
            anal_opts.osc_fit.plot_fit_corr=false;

            anal_opts.osc_fit.global=anal_opts.global;
            batch_data.osc_fit=fit_trap_freq(anal_opts.osc_fit,batch_data);

            batch_data.osc_fit.trap_freq_recons=nan*batch_data.osc_fit.ok.did_fits;
            mask=batch_data.osc_fit.ok.did_fits;
            batch_data.osc_fit.trap_freq_recons(mask)=3*(1/anal_opts.atom_laser.pulsedt)+batch_data.osc_fit.model_coefs(mask,2,1);

            trap_freq_history.shot_num=[trap_freq_history.shot_num,anal_opts.tdc_import.shot_num(mask)];
            trap_freq_history.trap_freq_val=[trap_freq_history.trap_freq_val,batch_data.osc_fit.trap_freq_recons(mask)];
            trap_freq_history.trap_freq_unc=[trap_freq_history.trap_freq_unc,batch_data.osc_fit.model_coefs(mask,2,2)'];

            %trim the history vectors
            if numel(trap_freq_history.shot_num)>anal_opts.history.shots
                %bit sloppy but will assume they are the same length
                trap_freq_history.shot_num=trap_freq_history.shot_num(end-anal_opts.history.shots:end);
                trap_freq_history.trap_freq_val=trap_freq_history.trap_freq_val(end-anal_opts.history.shots:end);
                trap_freq_history.trap_freq_unc=trap_freq_history.trap_freq_unc(end-anal_opts.history.shots:end);
            end

            sfigure(1);
            errorbar(trap_freq_history.shot_num,...
                trap_freq_history.trap_freq_val,trap_freq_history.trap_freq_unc,...
                'kx-','MarkerSize',7,'CapSize',0,'LineWidth',1.5)
            grid on
            h=gca;
            grid on    % turn on major grid lines
            grid minor % turn on minor grid lines
            % Set limits and grid spacing separately for the two directions:
            % Must set major grid line properties for both directions simultaneously:
            h.GridLineStyle='-'; % the default is some dotted pattern, I prefer solid
            h.GridAlpha=1;  % the default is partially transparent
            h.GridColor=[0,0,0]; % here's the color for the major grid lines
            % Idem for minor grid line properties:
            h.MinorGridLineStyle='-';
            h.MinorGridAlpha=0.1;
            h.MinorGridColor=[0,0,0]; % here's the color for the minor grid lines
            xlabel('Shot Number')
            ylabel('Fit Trap Freq')
            pause(1e-6)
            saveas(gcf,fullfile(anal_out.dir,'freq_history.png'))

        end
    end

end
