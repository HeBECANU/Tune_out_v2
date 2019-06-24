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


anal_opts.global.fall_time=0.417;
anal_opts.global.qe=0.09;

anal_opts.trig_dld=20.3;
anal_opts.dld_aquire=4;
anal_opts.trig_ai_in=20;

anal_opts.t_window=[0.45,2];
anal_opts.T_bin_width=0.0005;
anal_opts.mod_freq=420;


anal_opts.history.shots=200;

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



T_bin_num=round((anal_opts.t_window(2)-anal_opts.t_window(1))/anal_opts.T_bin_width);
anal_opts.T_bin_num=2*floor( T_bin_num/2)+1; %round to an odd number
anal_opts.T_bin_width=(anal_opts.t_window(2)-anal_opts.t_window(1))/T_bin_num; %change the value to the rounded one
anal_opts.T_bins=linspace(anal_opts.t_window(1),anal_opts.t_window(1)+anal_opts.T_bin_width*anal_opts.T_bin_num,anal_opts.T_bin_num);
%%
trap_mod_history=[];
trap_mod_history.trap_mod_val=[];
trap_mod_history.shot_num=[];
sfigure(1);
set(gcf,'color','w')
clf;

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
        ~ismember(anal_opts.tdc_import.shot_num, trap_mod_history.shot_num ) );

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
            if sum(batch_data.mcp_tdc.all_ok)>2
                anal_opts.anal_opts.plot_each=false;
            else
                anal_opts.anal_opts.plot_each=true;
            end

            iimax=numel(batch_data.mcp_tdc.shot_num);
            for ii=1:iimax
                txy_shot=batch_data.mcp_tdc.counts_txy{ii};
                if numel(txy_shot)>0
                    [T1d_counts,edges]=histcounts(txy_shot(:,1),anal_opts.T_bins);
                    t_centers=mean([edges(1:end-1);edges(2:end)]);
                    T1d_counts=1e-3*T1d_counts/(T_bin_width);
                    fft_out = fft_tx(t_centers,T1d_counts,'padding',1,'window','hamming');
                    batch_data.fft(ii,:,:)=fft_out;
                    [~,idx]=min(abs(anal_opts.mod_freq-fft_out(1,:)));
                    batch_data.mod_amp(ii)=fft_out(2,idx);
                else
                    batch_data.mod_amp(ii)=nan;
                end
            end
                     
            trap_mod_history.shot_num=[trap_mod_history.shot_num,anal_opts.tdc_import.shot_num];
            trap_mod_history.trap_mod_val=[trap_mod_history.trap_mod_val,batch_data.mod_amp];

            %plot(squeeze(batch_data.fft(1,1,:)),squeeze(mean(abs(batch_data.fft(:,2,:)),1)))
            %trim the history vectors
            if numel(trap_mod_history.shot_num)>anal_opts.history.shots
                %bit sloppy but will assume they are the same length
                trap_mod_history.shot_num=trap_mod_history.shot_num(end-anal_opts.history.shots:end);
                trap_mod_history.trap_mod_val=trap_mod_history.trap_mod_val(end-anal_opts.history.shots:end);
            end
            %%
            sfigure(1);
            subplot(3,1,1)
             plot(trap_mod_history.shot_num,...
                smooth(abs(trap_mod_history.trap_mod_val),5),...
                'r-','MarkerSize',7,'LineWidth',1)
            hold on
            plot(trap_mod_history.shot_num,...
                abs(trap_mod_history.trap_mod_val),...
                'kx-','MarkerSize',7,'LineWidth',1)
            hold off
            xlabel('Shot Number')
            ylabel('Mod amp')
            add_sig_lines(abs(trap_mod_history.trap_mod_val));
            yl=ylim;
            ylim([0,yl(2)])
            subplot(3,1,2)
            plot(trap_mod_history.shot_num,...
                smooth(real(trap_mod_history.trap_mod_val),5),...
                'r-','MarkerSize',7,'LineWidth',1)
            hold on
            plot(trap_mod_history.shot_num,...
                real(trap_mod_history.trap_mod_val),...
                'kx-','MarkerSize',7,'LineWidth',1)
            hold off
            xlabel('Shot Number')
            ylabel('Mod IN phase')
            add_sig_lines(real(trap_mod_history.trap_mod_val));
            subplot(3,1,3)
            plot(trap_mod_history.shot_num,...
                smooth(imag(trap_mod_history.trap_mod_val),5),...
                'r-','MarkerSize',7,'LineWidth',1)
            hold on
            plot(trap_mod_history.shot_num,...
                imag(trap_mod_history.trap_mod_val),...
                'kx-','MarkerSize',7,'LineWidth',1)
            hold off
            xlabel('Shot Number')
            ylabel('Mod OUT PHASE')
            add_sig_lines(imag(trap_mod_history.trap_mod_val));
            pause(1e-6)
            
            

        end
    end

end


function add_sig_lines(signal)
    xl=xlim;
    mean_std_sig=[nanmean(signal),nanstd(signal)];
    line(xl,[1,1]*...
        mean_std_sig(1)+mean_std_sig(2),'Color' ,'b')
    line(xl,[1,1]*...
        mean_std_sig(1)-mean_std_sig(2),'Color' ,'b')
    line(xl,[1,1]*...
        mean_std_sig(1)+mean_std_sig(2)*2,'Color','g')
    line(xl,[1,1]*...
        mean_std_sig(1)-mean_std_sig(2)*2,'Color','g')
end
