




%%

function osc_fit=fit_trap_freq_mod(anal_opts_osc_fit,data)
    %using the binned data we fit the trap freq to each shot
    %loop over every shot in data.mcp_tdc but output nothing if
    %data.mcp_tdc.all_ok(ii)=false
    %for compactness could use max((1:numel(data.mcp_tdc.all_ok)).*data.mcp_tdc.all_ok');
    %find the last good shot, but this will fuck up the mask && mask operations
    %later
    num_pulses = size(data.mcp_tdc.al_pulses.pos_stat,2);
    
%         num_segments = anal_opts_osc_fit.num_seg;
%         segment_size = (num_pulses-mod(num_pulses,num_segments))/num_segments;
    seg_data = anal_opts_osc_fit.segment_borders;
    num_segments = size(seg_data,1);
    
    %% MODIFIED OUTPUT
    % Now returns osc_fit.mod_fit_coefs with an extra dimension,
    % corresponding to the segment that is fitted. To look for anharmonic
    % effects, we should collect these by amplitude and see if the fitted
    % frequency varies. 
    iimax=size(data.mcp_tdc.counts_txy,2); 
     %try and guess the trap freq based on the peak of the fft, needed when the
     %kick amp is changing
    %ignore some of the fit errors
    warning('off','stats:nlinfit:ModelConstantWRTParam');
    warning('off','MATLAB:rankDeficientMatrix');
    osc_fit=[]; %clear the output struct
    %prealocate so that i can do the logics later
    osc_fit.dld_shot_idx=nan(iimax,1);
    osc_fit.model_coefs=nan(iimax,8,2);
    osc_fit.fit_rmse=nan(iimax,1);
    osc_fit.model=cell(iimax,num_segments);
    fprintf('Fitting oscillations in shots %04i:%04i',iimax,0)

     if anal_opts_osc_fit.plot_fits
                sfigure(51);
                clf
                set(gcf,'color','w')
     end



%     segment_freqs = zeros(num_segments);
%     iimax = 50
    for ii=1:iimax
        %position that data appears in data.mcp_tdc, not ness shot number
        %specified because we will remove any elements of osc_fit that did not
        %fit because of all_ok condition
        osc_fit.dld_shot_idx(ii)=ii;
        %shot number eg d123.txt as recorded by the tdc computer, not ness lv
        %number
        dld_shot_num=data.mcp_tdc.shot_num(ii);
        if data.mcp_tdc.all_ok(ii)
            osc_fit.dld_shot_num(ii)=dld_shot_num;
            %construct a more convinent temp variable txyz_tmp wich is the position in mm for use in the fit
            x_tmp=1e3*squeeze(data.mcp_tdc.al_pulses.pos_stat(ii,:,2));
            x_tmp=x_tmp-nanmean(x_tmp);
            y_tmp=1e3*squeeze(data.mcp_tdc.al_pulses.pos_stat(ii,:,3));
            y_tmp=y_tmp-nanmean(y_tmp);
            z_tmp=data.mcp_tdc.al_pulses.time(:)-squeeze(data.mcp_tdc.al_pulses.pos_stat(ii,:,1))';
            z_tmp=z_tmp'*anal_opts_osc_fit.global.velocity*1e3;
            z_tmp=z_tmp-nanmean(z_tmp);
            txyz_tmp_master=[data.mcp_tdc.al_pulses.time';x_tmp;y_tmp;z_tmp];
            sqrtn=sqrt(data.mcp_tdc.al_pulses.num_counts(ii,:)); %find the statistical uncert in a single shot
            xerr=1e3*squeeze(data.mcp_tdc.al_pulses.pos_stat(ii,:,5))./sqrtn;
            yerr=1e3*squeeze(data.mcp_tdc.al_pulses.pos_stat(ii,:,6))./sqrtn;
            zerr=1e3*squeeze(data.mcp_tdc.al_pulses.pos_stat(ii,:,4))*anal_opts_osc_fit.global.velocity./sqrtn;
            xyzerr_tmp_master=[xerr;yerr;zerr];
            xyzerr_tmp_master(:,sqrtn<2)=nan;
            %remove any data pts with nan position
            mask_master=(sum(isnan(txyz_tmp_master),1)==0);
            for segment = 1:num_segments
                
                %segment mask
%                 mask = mask_master&(1:num_pulses>segment_size*(segment-1))...
%                     &(1:num_pulses<=segment_size*(segment));
                mask = mask_master&(1:num_pulses>seg_data(segment,1))...
                    &(1:num_pulses<=seg_data(segment,2));                
                
                xyzerr_tmp=xyzerr_tmp_master(:,mask);
                txyz_tmp=txyz_tmp_master(:,mask);
                txyz_tmp(1,:)=txyz_tmp(1,:)-min(txyz_tmp(1,:));
                %try to find the peak osc freq to start the fit there
                if anal_opts_osc_fit.adaptive_freq
                    out=fft_tx(txyz_tmp(1,:),txyz_tmp(anal_opts_osc_fit.dimesion+1,:),10);
                    [~,nearest_idx]=max(abs(out(2,:)));
                    fit_freq=out(1,nearest_idx);
                    %fft_phase=angle(out(2,nearest_idx))+0.535;
                else
                    fit_freq=anal_opts_osc_fit.appr_osc_freq_guess(anal_opts_osc_fit.dimesion);
                end
                modelfun = @(b,x) exp(-x(:,1).*max(0,b(7))).*b(1).*sin(b(2)*x(:,1)*pi*2+b(3)*pi*2)+b(4)+b(8)*x(:,1)+0*b(5)*x(:,2)+0*b(6)*x(:,3);
                beta0=[std(txyz_tmp(anal_opts_osc_fit.dimesion+1,:))*8, fit_freq, 0, 1,0,0,2,0.01];
                cof_names={'amp','freq','phase','offset','ycpl','zcpl','damp','grad'};
                opt = statset('TolFun',1e-10,'TolX',1e-10,'MaxIter',1e4,...
                    'UseParallel',1);
                %select the aproapriate values to go in the response variable
                idx=1:4;
                idx(anal_opts_osc_fit.dimesion+1)=[];
                predictor=txyz_tmp(idx,:)';
                weights=1./(xyzerr_tmp(anal_opts_osc_fit.dimesion,:).^2);
                weights(isnan(weights))=1e-20; %if nan then set to smallest value you can
                weights=weights/sum(weights);
                %predictor=[tvalues,xvalues,zvalues];
                
                %% Fit trials
                num_variations = 2;
                fitobject_temp = cell(num_variations+1,1);
                beta_stash =beta0.*([ones(size(beta0));1+0.1*randn(num_variations,size(beta0,2))]);
                for i=1:num_variations+1
                    beta_fit = beta_stash(i,:);
                    fitobject_temp{i}=fitnlm(predictor,txyz_tmp(anal_opts_osc_fit.dimesion+1,:)',modelfun,beta_fit,...
                        'Weights',weights,'options',opt,...
                        'CoefficientNames',cof_names);
                end
                
                errs = cellfun(@(x) x.MSE,fitobject_temp);
                bestidx = find(min(errs) == errs);
                fitobject = fitobject_temp{bestidx};
                osc_fit.model{ii,segment}=fitobject;
                %%
%                 fitobject=fitnlm(predictor,txyz_tmp(anal_opts_osc_fit.dimesion+1,:)',modelfun,beta0,...
%                     'Weights',weights,'options',opt,...
%                     'CoefficientNames',cof_names);

                fitparam=fitobject.Coefficients;
                
                osc_fit.model_coefs(ii,:,:,segment)=[fitparam.Estimate,fitparam.SE];
                osc_fit.fit_rmse(ii,segment)=fitobject.RMSE;
    
            
    %             osc_fit.model{ii} = model_store;
                %limiting frequnecy prediction from http://adsabs.harvard.edu/full/1999DSSN...13...28M
                %meanwidth=sqrt(mean(squeeze(data.mcp_tdc.al_pulses.pos_stat(ii,:,5)).^2))*1e3;
%                 frequnclim=sqrt(6/sum(data.mcp_tdc.al_pulses.num_counts(ii,:)))*...
%                     (1/(pi*range(data.mcp_tdc.al_pulses.time)))*...
%                     (meanwidth/fitparam{2,1});
                %fprintf('sampling limit %2.3g Hz, fint unc %2.3g Hz, ratio %2.3g \n',[frequnclim,fitparam{2,2},fitparam{2,2}/frequnclim])
                %osc_fit.fit_sample_limit{ii,segment}=[frequnclim,fitparam{2,2},fitparam{2,2}/frequnclim];
                if anal_opts_osc_fit.plot_fits
                    tplotvalues=linspace(min(data.mcp_tdc.al_pulses.time),...
                        max(data.mcp_tdc.al_pulses.time),1e5)';
                    predictorplot=[tplotvalues,...
                               interp1(predictor(:,1),predictor(:,2),tplotvalues),...
                               interp1(predictor(:,1),predictor(:,3),tplotvalues)];
                    [prediction,ci]=predict(fitobject,predictorplot);
                    sfigure(51);
                    subplot(2,1,1)
                    plot(txyz_tmp(1,:),txyz_tmp(2,:),'kx-')
                    hold on
                    plot(txyz_tmp(1,:),txyz_tmp(3,:),'rx-')
                    plot(txyz_tmp(1,:),txyz_tmp(4,:),'bx-')
                    title(sprintf('%u',segment))
                    hold off
                    ylabel('X Pos (mm)')
                    xlabel('Time (s)')
                    set(gca,'Ydir','normal')
                    set(gcf,'Color',[1 1 1]);
                    legend('x','y','z')

                    subplot(2,1,2)
                    plot(predictorplot(:,1),prediction,'-','LineWidth',1.5,'Color',[0.5 0.5 0.5])
                    ax = gca;
                    set(ax, {'XColor', 'YColor'}, {'k', 'k'});
                    hold on
                    plot(predictorplot(:,1),ci(:,1),'-','LineWidth',1.5,'Color','k')
                    plot(predictorplot(:,1),ci(:,2),'-','LineWidth',1.5,'Color','k')
                    errorbar(predictor(:,1),txyz_tmp(anal_opts_osc_fit.dimesion+1,:)',xyzerr_tmp(anal_opts_osc_fit.dimesion,:),'k.','MarkerSize',10,'CapSize',0,'LineWidth',1,'Color','r') 
                    set(gcf,'Color',[1 1 1]);
                    ylabel('X(mm)')
                    xlabel('Time (s)')
                    hold off
                    ax = gca;
                    set(ax, {'XColor', 'YColor'}, {'k', 'k'});
                    set(gca,'linewidth',1.0)
%                     saveas(gca,sprintf('%sfit_dld_shot_num%04u.png',anal_opts_osc_fit.global.out_dir,dld_shot_num))
                    pause(1e-5)
                end%PLOTS
            end%segments
        end
        fprintf('\b\b\b\b%04u',ii)
    end
    fprintf('...Done\n')
    %if the above did a fit then set the element in the logic vector to true

%     for segment = 1:num_segments
%         osc_fit.ok.did_fits(segment)=~cellfun(@(x) isequal(x,[]),osc_fit.model{:}{segment});
% 
%         %% look for failed fits
%         %look at the distribution of fit errors
%         mean_fit_rmse=nanmean(osc_fit.fit_rmse(osc_fit.ok.did_fits));
%         median_fit_rmse=nanmedian(osc_fit.fit_rmse(osc_fit.ok.did_fits));
%         std_fit_rmse=nanstd(osc_fit.fit_rmse(osc_fit.ok.did_fits));
%         fprintf('fit error: median %f mean %f std %f\n',...
%            median_fit_rmse,mean_fit_rmse,std_fit_rmse)
% 
%         %label anything std_cut_fac*std away from the median as a bad fit
%         std_cut_fac=4;
%         mask=osc_fit.ok.did_fits;
%         osc_fit.ok.rmse(segment)=mask & osc_fit.fit_rmse...
%             < median_fit_rmse+std_cut_fac*std_fit_rmse;
%     end
    %label the fits with trap freq more than 1hz awaay from the mean as bad
    %osc_fit.ok.rmse(osc_fit.ok.rmse)=abs(osc_fit.model_coefs(osc_fit.ok.rmse,2,1)'...
    %    -mean(osc_fit.model_coefs(osc_fit.ok.rmse,2,1)))<1;
    
%     if anal_opts_osc_fit.plot_err_history
%         figure(52)
%         clf
%         set(gcf,'color','w')
%         subplot(2,1,1)
%         hist(osc_fit.fit_rmse(osc_fit.ok.did_fits))
%         xlabel('RMSE')
%         ylabel('counts')
%         yl=ylim;
%         hold on
%         line([1,1]*mean_fit_rmse,yl,'Color','k')
%         line([1,1]*median_fit_rmse,yl,'Color','m')
%         line([1,1]*(mean_fit_rmse+std_cut_fac*std_fit_rmse),yl,'Color','r')
%         line([1,1]*(mean_fit_rmse-std_cut_fac*std_fit_rmse),yl,'Color','r')
%         hold off
%     
%         subplot(2,1,2)
%         plot(osc_fit.fit_rmse(osc_fit.ok.did_fits))
%         xlabel('shot idx')
%         ylabel('RMSE')
%         saveas(gca,sprintf('%sfits_rmse_history.png',anal_opts_osc_fit.global.out_dir))
%         
%     end
%     
%     
%     if anal_opts_osc_fit.plot_fit_corr
%         sfigure(53);
%         clf
%         set(gcf,'color','w')
%         %see if the fit error depends on the probe freq
%         subplot(3,3,1)
%         tmp_probe_freq=data.mcp_tdc.probe.freq.act.mean(osc_fit.ok.rmse);
%         tmp_probe_freq=(tmp_probe_freq-nanmean(tmp_probe_freq))*1e-3;
%         plot(tmp_probe_freq,...
%             osc_fit.fit_rmse(osc_fit.ok.rmse),'xk')
%         xlabel('probe beam freq (GHz)')
%         ylabel('fit error')
%         title('Fit Error')
%         %see if the damping changes
%         subplot(3,3,2)
%         plot(tmp_probe_freq,...
%            osc_fit.model_coefs((osc_fit.ok.rmse),1,1),'xk')
%         xlabel('probe beam freq (GHz)')
%         ylabel('Osc amp (mm)')
%         title('Amp')
%         subplot(3,3,3)
%         plot(tmp_probe_freq,...
%            osc_fit.model_coefs((osc_fit.ok.rmse),3,1),'xk')
%         xlabel('probe beam freq (GHz)')
%         ylabel('Phase (rad)')
%         title('Phase')
%         subplot(3,3,4)
%         plot(tmp_probe_freq,...
%            osc_fit.model_coefs((osc_fit.ok.rmse),4,1),'xk')
%         xlabel('probe beam freq (GHz)')
%         ylabel('offset')
%         title('offset')
%         subplot(3,3,5)
%         plot(tmp_probe_freq,...
%            osc_fit.model_coefs((osc_fit.ok.rmse),5,1),'xk')
%         xlabel('probe beam freq (GHz)')
%         ylabel('ycpl')
%         title('ycpl')
%         subplot(3,3,6)
%         plot(tmp_probe_freq,...
%            osc_fit.model_coefs((osc_fit.ok.rmse),6,1),'xk')
%         xlabel('probe beam freq (GHz)')
%         ylabel('zcpl')
%         title('zcpl')
%         %see if the amp changes
%         subplot(3,3,7)
%         plot(tmp_probe_freq,...
%            1./osc_fit.model_coefs((osc_fit.ok.rmse),7,1),'xk')
%         xlabel('probe beam freq (GHz)')
%         ylabel('Damping Time (s)')
%         title('Damping')
%         subplot(3,3,8)
%         plot(tmp_probe_freq,...
%            osc_fit.model_coefs((osc_fit.ok.rmse),8,1),'xk')
%         xlabel('probe beam freq (GHz)')
%         ylabel('Grad mm/s')
%         title('Grad')
%         %look for anharminicity with a osc amp/freq correlation
%         subplot(3,3,9)
%         plot(abs(osc_fit.model_coefs((osc_fit.ok.rmse),1,1)),...
%             osc_fit.model_coefs((osc_fit.ok.rmse),2,1),'xk')
%         xlabel('osc (mm)')
%         ylabel('fit freq (Hz)')
%         title('anharmonicity')
%     
%         set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])
%         plot_name='fit_correlations';
%         saveas(gcf,[anal_opts_osc_fit.global.out_dir,plot_name,'.png'])
%         saveas(gcf,[anal_opts_osc_fit.global.out_dir,plot_name,'.fig'])
%     end



end