%function out=fit_temperature(anal_opts,data)
%anal_opts.xylim=anal_opts.tdc_import.txylim(2:3,:);


pulse_num=1;
time_cen=data.mcp_tdc.al_pulses.time_cen(pulse_num);
time_width=data.mcp_tdc.al_pulses.pulsedt*0.2;

% subtract the com from each 

% for ii=1:num_shots
%     txy_shot=data.mcp_tdc.counts_txy{ii};
%     if ~isempty(txy_shot)
%         txy_shot=masktxy_square(txy_shot,anal_opts.hotspot_mask.square_mask);
%         txy_shot=masktxy_2d_circle(txy_shot,anal_opts.hotspot_mask.circ_mask);
%         data.mcp_tdc.masked.num_counts(ii)=numel(txy_shot);
%         data.mcp_tdc.masked.counts_txy{ii}=txy_shot;
%     else
%         warning('empty shot')
%     end
%     if mod(ii,10)==0,fprintf('\b\b\b\b%04u',ii),end 
% end


all_counts=cat(1,data.mcp_tdc.counts_txy{:});
tlim=time_cen+[-1,1]*0.5*time_width;
tmp_xlim=[-50e-3, 50e-3];    
tmp_ylim=[-50e-3, 50e-3];
anal_opts=[];
txy_window=[tlim;tmp_xlim;tmp_ylim];
all_counts_masked=masktxy_square(all_counts,txy_window);

%
s1=smooth_hist(all_counts_masked(:,2),'sigma',1e-4)

stfig('asdf')
mask=s1.count_rate.smooth<1e5;
plot(s1.bin.centers(mask),s1.count_rate.smooth(mask))

%%



%end