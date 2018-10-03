%trying some more things to see if i can decrease the dependence of the
%probe beam on the fit error
%there are a few effects that could be coming into play here
%   -number dependent damping
%   -probe beam based damping
%   -amplitude dependent trap frequency
%i cant see either as causing the dependence of the probe beam freq on the
%fit error which implies
%   -that the fit error depends on frequency rather steeply
%   -probe beam based damping
%the latter seems to be the most likely

%LOG
%2018-10-02 01:30 
    %tried adding an extra term to the damping as a 1st order
    %expansion of a non constant damping time
    % @(b,x) exp(-x(:,1).*(max(0,b(7))+x(:,1).*b(8))).*b(1).*sin(b(2)*x(:,1)*pi*2+b(3)*pi*2)+b(4)+b(8)*x(:,1)+b(5)*x(:,2)+b(6)*x(:,3);
    %prev mean fit err 0.051912 now 0.051725


%%
iimax=size(data.mcp_tdc.counts_txy,2); 
do_plots_atom_num_fit=false;
pulse_num_min=1;
pulse_num_max=10;
figure(1)
clf
set(gcf,'color','w')
fit2_res=nan(iimax,2,2);
fprintf('Fitting atom num in shots %04i:%04i',iimax,0)
for ii=1:iimax
    %position that data appears in data.mcp_tdc, not ness shot number
    %specified because we will remove any elements of osc_fit that did not
    %fit because of all_ok condition
    datatest.num_fit.dld_shot_idx(ii)=ii;
    %shot number eg d123.txt as recorded by the tdc computer, not ness lv
    %number
    dld_shot_num=data.mcp_tdc.shot_num(ii);
    if data.mcp_tdc.all_ok(ii)
        counts=data.mcp_tdc.al_pulses.num_counts(ii,:);
        pulse_number=1:numel(counts);
        pulse_number=pulse_number(pulse_num_min:pulse_num_max);
        counts=counts(pulse_num_min:pulse_num_max);
        
        fo = statset('TolFun',1e-30,...
            'TolX',1e-50,...
            'MaxIter',1e7,...
            'UseParallel',1); %'DerivStep',1e-5
        modelfun = @(b,x) b(1)*b(2)*((1-b(2)).^(x(:,1)-1));
        beta0=[1e6,1e-2];
        fit_mdl_pow=fitnlm(pulse_number,counts,modelfun,...
           beta0,'CoefficientNames',{'N0','eta'},'Options',fo);  
        
        fit2_res(ii,:,:)=[[fit_mdl_pow.Coefficients.Estimate(2),fit_mdl_pow.Coefficients.SE(2)];...
        [fit_mdl_pow.Coefficients.Estimate(1)/anal_opts.qe,...
        fit_mdl_pow.Coefficients.SE(1)/anal_opts.qe]];
    
        if do_plots_atom_num_fit
            xplot=linspace(min(pulse_number),max(pulse_number),1e3)';
            [ymodel,ci]=predict(fit_mdl_pow,xplot,'Prediction','observation','Alpha',0.3174);
            plot(pulse_number,counts,'kx')
            hold on
            plot(xplot,ymodel,'r-')
            plot(xplot,ci(:,2),'m-')
            plot(xplot,ci(:,1),'m-')
            hold off
            set(gca, 'YScale', 'log')
            xlabel('pulse number')
            ylabel('detected atoms')
            pause(1e-2)
        end
        
    end
    fprintf('\b\b\b\b%04u',ii)
end
fprintf('...Done\n')



%%
tmp_atom_num_fit=fit2_res(data.mcp_tdc.all_ok,2,1);
tmp_out_frac_fit=fit2_res(data.mcp_tdc.all_ok,1,1);
fprintf('out frac mean %.2g,med %.2g sd %.2g norm sd %.2f\n',...
    mean(tmp_out_frac_fit),median(tmp_out_frac_fit),std(tmp_out_frac_fit),...
    std(tmp_out_frac_fit)/mean(tmp_out_frac_fit))
fprintf('atom number mean %.2g,med %.2g sd %.2g norm sd %.2f\n',...
    mean(tmp_atom_num_fit),median(tmp_atom_num_fit),std(tmp_atom_num_fit),...
    std(tmp_atom_num_fit)/mean(tmp_atom_num_fit))
fprintf('predicted mean sd from fits %.2g\n',mean(fit2_res(data.mcp_tdc.all_ok,2,2)))
%find the ratio to just counting
tmp_atom_num_count=data.mcp_tdc.num_counts(data.mcp_tdc.all_ok)/anal_opts.qe;
fprintf('mean ratio of fit method vs counting %.3f\n',mean((tmp_atom_num_count)./tmp_atom_num_fit'))


figure(2)
subplot(2,1,1)
set(gcf,'color','w')
plot(tmp_atom_num_fit)
title('historical count trend')
xlabel('shot idx')
ylabel('BEC number')
subplot(2,1,2)
plot((data.mcp_tdc.num_counts(data.mcp_tdc.all_ok)/anal_opts.qe),tmp_atom_num_fit','x')
title('comparing atom num methods')
xlabel('atom num fit')
ylabel('atom num total counts')

%% Fitting the Trap freq
%using the binned data we fit the trap freq to each shot
%loop over every shot in data.mcp_tdc but output nothing if
%data.mcp_tdc.all_ok(ii)=false
%for compactness could use max((1:numel(data.mcp_tdc.all_ok)).*data.mcp_tdc.all_ok');
%find the last good shot, but this fuck up the mask && mask operations
%later
iimax=size(data.mcp_tdc.counts_txy,2); 
plots=false;
anal_opts.atom_laser.appr_osc_freq_guess=[52,46.7,40];
 %try and guess the trap freq based on the peak of the fft, needed when the
 %kick amp is changing
adaptive_fit_freq=true;
%ignore some of the fit errors
warning('off','stats:nlinfit:ModelConstantWRTParam');
warning('off','MATLAB:rankDeficientMatrix');
data.osc_fit=[]; %clear the output struct
%prealocate so that i can do the logics later
data.osc_fit.dld_shot_idx=nan(1,iimax);
data.osc_fit.model_coefs=nan(iimax,9,2);
data.osc_fit.fit_rmse=nan(1,iimax);
data.osc_fit.model=cell(1,iimax);
fprintf('Fitting oscillations in shots %04i:%04i',iimax,0)
for ii=1:iimax
    %position that data appears in data.mcp_tdc, not ness shot number
    %specified because we will remove any elements of osc_fit that did not
    %fit because of all_ok condition
    data.osc_fit.dld_shot_idx(ii)=ii;
    %shot number eg d123.txt as recorded by the tdc computer, not ness lv
    %number
    dld_shot_num=data.mcp_tdc.shot_num(ii);
    if data.mcp_tdc.all_ok(ii)
        data.osc_fit.dld_shot_num(ii)=dld_shot_num;
        %construct a more convinent temp variable txyz_tmp wich is the position in mm for use in the fit
        x_tmp=1e3*squeeze(data.mcp_tdc.al_pulses.pos_stat(ii,:,2));
        x_tmp=x_tmp-nanmean(x_tmp);
        y_tmp=1e3*squeeze(data.mcp_tdc.al_pulses.pos_stat(ii,:,3));
        y_tmp=y_tmp-nanmean(y_tmp);
        z_tmp=data.mcp_tdc.al_pulses.time-squeeze(data.mcp_tdc.al_pulses.pos_stat(ii,:,1))';
        z_tmp=z_tmp'*anal_opts.velocity*1e3;
        z_tmp=z_tmp-nanmean(z_tmp);
        txyz_tmp=[data.mcp_tdc.al_pulses.time';x_tmp;y_tmp;z_tmp];
        sqrtn=sqrt(data.mcp_tdc.al_pulses.num_counts(ii,:)); %find the statistical uncert in a single shot
        xerr=1e3*squeeze(data.mcp_tdc.al_pulses.pos_stat(ii,:,5))./sqrtn;
        yerr=1e3*squeeze(data.mcp_tdc.al_pulses.pos_stat(ii,:,6))./sqrtn;
        zerr=1e3*squeeze(data.mcp_tdc.al_pulses.pos_stat(ii,:,4))*anal_opts.velocity./sqrtn;
        xyzerr_tmp=[xerr;yerr;zerr];
        xyzerr_tmp(:,sqrtn<2)=nan;

        %remove any data pts with nan position
        mask=sum(isnan(txyz_tmp),1)==0;
        xyzerr_tmp=xyzerr_tmp(:,mask);
        txyz_tmp=txyz_tmp(:,mask);

        %try to find the peak osc freq to start the fit there
        if adaptive_fit_freq
            out=fft_tx(txyz_tmp(1,:),txyz_tmp(histplot.dimesion+1,:),10);
            [~,nearest_idx]=max(abs(out(2,:)));
            fit_freq=out(1,nearest_idx);
            %fft_phase=angle(out(2,nearest_idx))+0.535;
        else
            fit_freq=anal_opts.atom_laser.appr_osc_freq_guess(histplot.dimesion);
        end
        modelfun = @(b,x) exp(-x(:,1).*max([-10,b(7)])-(x(:,1).^2).*max([-10,b(9)])).*...
            b(1).*sin(b(2)*x(:,1)*pi*2+b(3)*pi*2)...
            +b(4)+b(8)*x(:,1)+b(5)*x(:,2)+b(6)*x(:,3);
        beta0=[std(txyz_tmp(histplot.dimesion+1,:))*8, fit_freq, 0, 1,0,0,2,0.01,1e-5];
        cof_names={'amp','freq','phase','offset','ycpl','zcpl','damp','grad','d2'};
        opt = statset('TolFun',1e-10,'TolX',1e-10,'MaxIter',1e4,...
            'UseParallel',1);
        %select the aproapriate values to go in the response variable
        idx=1:4;
        idx(histplot.dimesion+1)=[];
        predictor=txyz_tmp(idx,:)';
        response=txyz_tmp(histplot.dimesion+1,:)';
        weights=1./(xyzerr_tmp(histplot.dimesion,:).^2);
        weights(isnan(weights))=1e-20; %if nan then set to smallest value you can
        weights=weights/sum(weights);
        %predictor=[tvalues,xvalues,zvalues];
        fitobject=fitnlm(predictor,response,modelfun,beta0,...
            'Weights',weights,'options',opt,...
            'CoefficientNames',cof_names);
        data.osc_fit.model{ii}=fitobject;
        fitparam=fitobject.Coefficients;
        data.osc_fit.model_coefs(ii,:,:)=[fitparam.Estimate,fitparam.SE];
        data.osc_fit.fit_rmse(ii)=fitobject.RMSE;
        %limiting frequnecy prediction from http://adsabs.harvard.edu/full/1999DSSN...13...28M
        meanwidth=sqrt(mean(squeeze(data.mcp_tdc.al_pulses.pos_stat(ii,:,5)).^2))*1e3;
        frequnclim=sqrt(6/sum(data.mcp_tdc.al_pulses.num_counts(ii,:)))*...
            (1/(pi*range(data.mcp_tdc.al_pulses.time)))*...
            (meanwidth/fitparam{2,1});
        %fprintf('sampling limit %2.3g Hz, fint unc %2.3g Hz, ratio %2.3g \n',[frequnclim,fitparam{2,2},fitparam{2,2}/frequnclim])
        data.osc_fit.fit_sample_limit{ii}=[frequnclim,fitparam{2,2},fitparam{2,2}/frequnclim];
        if plots
            [prediction_smpl,~]=predict(fitobject,predictor);
            tplotvalues=linspace(min(data.mcp_tdc.al_pulses.time),max(data.mcp_tdc.al_pulses.time),1e5)';
            predictorplot_fine=[tplotvalues,...
                       interp1(predictor(:,1),predictor(:,2),tplotvalues),...
                       interp1(predictor(:,1),predictor(:,3),tplotvalues)];
            [prediction_fine,ci_fine]=predict(fitobject,predictorplot_fine);
            residulas=response-prediction_smpl;
        
            sfigure(2);
            clf
            set(gcf,'color','w')
            subplot(3,1,1)
            plot(txyz_tmp(1,:),txyz_tmp(2,:),'kx-')
            hold on
            plot(txyz_tmp(1,:),txyz_tmp(3,:),'rx-')
            plot(txyz_tmp(1,:),txyz_tmp(4,:),'bx-')
            hold off
            ylabel('X Pos (mm)')
            xlabel('Time (s)')
            set(gca,'Ydir','normal')
            set(gcf,'Color',[1 1 1]);
            legend('x','y','z')
            pause(0.05)

            subplot(3,1,2)
            plot(predictorplot_fine(:,1),prediction_fine,'-','LineWidth',1.5,'Color',[0.5 0.5 0.5])
            ax = gca;
            set(ax, {'XColor', 'YColor'}, {'k', 'k'});
            hold on
            plot(predictorplot_fine(:,1),ci_fine(:,1),'-','LineWidth',1.5,'Color','k')
            plot(predictorplot_fine(:,1),ci_fine(:,2),'-','LineWidth',1.5,'Color','k')
            errorbar(txyz_tmp(1,:),txyz_tmp(histplot.dimesion+1,:)',xyzerr_tmp(histplot.dimesion,:),'k.','MarkerSize',10,'CapSize',0,'LineWidth',1,'Color','r') 
            set(gcf,'Color',[1 1 1]);
            ylabel('X(mm)')
            xlabel('Time (s)')
            hold off
            ax = gca;
            set(ax, {'XColor', 'YColor'}, {'k', 'k'});
            set(gca,'linewidth',1.0)
            
            subplot(3,1,3)
            
            errorbar(predictor(:,1),residulas,xyzerr_tmp(histplot.dimesion,:),'k.','MarkerSize',10,'CapSize',0,'LineWidth',1,'Color','r') 
            ax = gca;
            set(ax, {'XColor', 'YColor'}, {'k', 'k'});
            hold on
            plot(predictorplot(:,1),ci_fine(:,1)-prediction_fine,'-','LineWidth',1.5,'Color','k')
            plot(predictorplot(:,1),ci_fine(:,2)-prediction_fine,'-','LineWidth',1.5,'Color','k')
            set(gcf,'Color',[1 1 1]);
            ylabel('X(mm)')
            xlabel('Time (s)')
            hold off
            ax = gca;
            set(ax, {'XColor', 'YColor'}, {'k', 'k'});
            set(gca,'linewidth',1.0)
            
            saveas(gca,sprintf('.\\out\\fit_dld_shot_num%04u.png',dld_shot_num))
        end% PLOTS
    end
    fprintf('\b\b\b\b%04u',ii)
end
fprintf('...Done\n')

data.osc_fit.ok.did_fits=~cellfun(@(x) isequal(x,[]),data.osc_fit.model);
        

%%
%look for failed fits
figure(1)
clf
subplot(2,1,1)
hist(data.osc_fit.fit_rmse(data.osc_fit.ok.did_fits))
xlabel('RMSE')
ylabel('counts')
subplot(2,1,2)
plot(data.osc_fit.fit_rmse(data.osc_fit.ok.did_fits))
xlabel('shot idx')
mask=data.osc_fit.ok.did_fits;
data.osc_fit.ok.rmse=mask & data.osc_fit.fit_rmse...
    < nanmean(data.osc_fit.fit_rmse(mask))+2*nanstd(data.osc_fit.fit_rmse(mask));
data.osc_fit.ok.rmse(data.osc_fit.ok.rmse)=abs(data.osc_fit.model_coefs(data.osc_fit.ok.rmse,2,1)'...
    -mean(data.osc_fit.model_coefs(data.osc_fit.ok.rmse,2,1)))<1;
%% CHECK IF FIT ERROR DEPENDENCE WITH PROBE BEAM
figure
clf

plot(data.mcp_tdc.probe.freq.act.mean(data.osc_fit.ok.rmse),...
    data.osc_fit.fit_rmse(data.osc_fit.ok.rmse),'x')
xlabel('probe beam freq')
ylabel('fit error')