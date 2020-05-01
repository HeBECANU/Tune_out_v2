function out=fit_atom_number(anal_opts_afit,data)
%this function take in the data structure and returns a estimate of the total atom number by fitting a decay function 
%to the number in the pulsed atom laser 

atom_fit_mdl = @(b,x) b(1)*b(2)*((1-b(2)).^(x(:,1)-1));
exp_qe=anal_opts_afit.qe;




iimax=size(data.mcp_tdc.counts_txy,2); 
fit_out_frac_ntot=nan(iimax,2,2); %value and uncert
fit_predict=cell(iimax,1);
fit_predict_unc=fit_predict;
fprintf('Fitting atom number in shots %04i:%04i',iimax,0)
pulse_num=anal_opts_afit.pulses(1):anal_opts_afit.pulses(2);
for ii=1:iimax
    if data.mcp_tdc.all_ok(ii)
        pulse_counts=data.mcp_tdc.al_pulses.num_counts(ii,pulse_num);

        coefs_guess=[data.mcp_tdc.num_counts(ii),1e-2];
        opts = statset('nlinfit');
        opts.RobustWgtFun = 'bisquare';
        opts.MaxIter=2e3;
        pulse_fit=fitnlm(pulse_num,pulse_counts,atom_fit_mdl,coefs_guess,'Options',opts);
        fit_out_frac_ntot(ii,:,:)=[[pulse_fit.Coefficients.Estimate(2),pulse_fit.Coefficients.SE(2)];...
                                    [pulse_fit.Coefficients.Estimate(1)/exp_qe,...
                                     pulse_fit.Coefficients.SE(1)/exp_qe]];
        %pass a function that can be used to get the atom number at any pulse number
        fit_predict{ii}=@(x) predict(pulse_fit,x)...
            ./pulse_fit.Coefficients.Estimate(2);
        fit_predict_unc{ii}=@(x) make_unc_fun(x,pulse_fit);
        if anal_opts_afit.plot.each_shot
            stfig('atom number fit','add_stack',1)
            fprintf('fit number in shot %u\n',ii)
            xplot=linspace(min(pulse_num),max(pulse_num),1e3)';
            [ymodel,ci]=predict(pulse_fit,xplot,'Prediction','observation','Alpha',1-erf(1/sqrt(2)));
            plot(pulse_num,pulse_counts,'kx')
            hold on
            plot(xplot,ymodel,'r-')
            plot(xplot,ci(:,2),'m-')
            plot(xplot,ci(:,1),'m-')
            hold off
            xlabel('Pulse Index')
            ylabel('Pulse Counts')
            set(gca, 'YScale', 'log')
            fprintf('frac predicted %s\n',...
                string_value_with_unc(fit_out_frac_ntot(ii,1,1),fit_out_frac_ntot(ii,1,2),'type','b'))
            fprintf('Total no. fit %s measured in file %.2e (qe %.1e corrected)\n',...
                string_value_with_unc(fit_out_frac_ntot(ii,2,1),fit_out_frac_ntot(ii,2,2),'type','b')...
                ,data.mcp_tdc.num_counts(ii)/exp_qe,exp_qe)
            pause(1e-9)
        end
    end
    if mod(ii,10)==0, fprintf('\b\b\b\b%04u',ii), end
end
fprintf('...Done\n')

out.fit_out_frac_ntot=fit_out_frac_ntot;
out.fit_predict=fit_predict;
out.fit_predict_unc=fit_predict_unc;
%%
if anal_opts_afit.plot.history
    stfig('atom number history','add_stack',1)
    clf
    set(gcf,'color','w')
    subplot(3,1,1)
    errorbar(data.mcp_tdc.shot_num,fit_out_frac_ntot(:,1,1),fit_out_frac_ntot(:,1,2),'k.','capsize',0)
    ylabel('Fit Out. Frac.')
    xlabel('Shot number')
    subplot(3,1,2)
    errorbar(data.mcp_tdc.shot_num,fit_out_frac_ntot(:,2,1),fit_out_frac_ntot(:,2,2),'k.','capsize',0)
    ylabel('Fit Total Num')
    xlabel('Shot number')
    subplot(3,1,3)
    raw_counts_atom_num=data.mcp_tdc.num_counts(:)/exp_qe;
    errorbar(data.mcp_tdc.shot_num,fit_out_frac_ntot(:,2,1)./raw_counts_atom_num,...
        fit_out_frac_ntot(:,2,2)./raw_counts_atom_num,'k.','capsize',0)
    xl=xlim;
    line(xl,[1,1],'color','k')
    ylabel('Num count/Fit')
    xlabel('Shot number')
end
%%
% we can quantify how much worse this is than shot noise by taking the shot noise in the number of counts in those pulses
% scaling it by the QE and comparing it to the derived uncert in the inital atom number
%data.num_fit.fit_out_frac_ntot(:,2,2)./(sqrt(nansum(data.mcp_tdc.al_pulses.num_counts(:,pulse_num),2))/exp_qe)
% this comes out to ~36 which is not too shabby


end


function fit_pred_unc=make_unc_fun(x,fit_obj)
    ci_one_sd=1-erf(1/sqrt(2));%one sd CI
    [~,ci]=predict(fit_obj,x,'Prediction' ,'curve','Alpha',ci_one_sd);
    fit_pred_unc=(range(ci)/2)./fit_obj.Coefficients.Estimate(2);
end