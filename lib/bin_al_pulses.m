function al_pulses=bin_al_pulses(anal_opt_al,data)

%tic
iimax=size(data.mcp_tdc.counts_txy,2);
al_pulses=[];
al_pulses.pulsedt=anal_opt_al.pulsedt;
al_pulses.window=nan(anal_opt_al.pulses,3,2); %initalize
al_pulses.num_counts=nan(iimax,anal_opt_al.pulses);
  
plots=false;
fprintf('binning pulses in files %04u:%04u',size(data.mcp_tdc.counts_txy,2),0)
first_good_shot=true;
for shot=1:iimax
        if data.mcp_tdc.all_ok(shot)
            for pulse=1:anal_opt_al.pulses
                %set up time window centered arround t0
                trange=anal_opt_al.t0+anal_opt_al.pulsedt...
                    *(anal_opt_al.start_pulse+pulse-2)+...
                    anal_opt_al.pulse_twindow*[-0.5,0.5];
                pulse_win_txy=[trange;anal_opt_al.xylim]; 
                counts_pulse=masktxy(data.mcp_tdc.counts_txy{shot},pulse_win_txy);
                if plots
                    sfigure(79);
                    set(gcf,'Color',[1 1 1]);
                    subplot(3,1,1)
                    hist(counts_pulse(:,1),100)
                    xlabel('t')
                    title('full')
                    subplot(3,1,2)
                    hist(counts_pulse(:,2),100)
                    xlabel('x')
                    title('full')
                    subplot(3,1,3)
                    hist(counts_pulse(:,3),100)
                    xlabel('y')
                    title('full')
                    pause(0.01)
                end
                if first_good_shot
                    %only need to store this on first shot becasue the same for
                    %all shots
                    al_pulses.window(pulse,:,:)=pulse_win_txy; 
                    al_pulses.time(pulse,:)=(pulse+anal_opt_al.start_pulse-2)*anal_opt_al.pulsedt;
                end
                al_pulses.num_counts(shot,pulse)=size(counts_pulse(:,3),1);
                al_pulses.pos_stat(shot,pulse,:)=[...
                                           mean(counts_pulse(:,1)),...
                                           mean(counts_pulse(:,2)),...
                                           mean(counts_pulse(:,3)),...
                                           std(counts_pulse(:,1)),...
                                           std(counts_pulse(:,2)),...
                                           std(counts_pulse(:,3))]; 
            end%pulse
        if first_good_shot,first_good_shot=false; end
        end%is data.mcp_tdc.all_ok
        if mod(shot,10)==0,fprintf('\b\b\b\b%04u',shot),end    
%to set the pulse t0 right it can be handy to uncomment the next line
%fprintf('\nmean time %3.5f            \n ',mean(al_pulses.pos_stat(shot,:,1)-al_pulses.time(:)'))
end%shots
%toc
fprintf('...Done\n') 


end