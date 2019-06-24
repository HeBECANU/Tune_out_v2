function fit_scan(data)

% sfigure(200);
% clf;


FFT=fft_tx(data.T,data.Y);




beta0 = [0.6e-3, 0.99,.0154,-1,0];

%Make plots

sfigure(100);
clf;
subplot(2,1,1)
plot(data.T,data.Y)
hold on
% plot(data.T,p_Airy(beta0,data.T),'k-.')
% fit_data = fitnlm(data.T,data.Y,@p_Airy,beta0);
% plot(data.T,p_Airy(fit_data.Coefficients.Estimate,data.T),'r-.')
title('Single scan fit attempt')
xlabel('PZ voltage')
ylabel('Voltage')
legend('Raw data','Fit seed','Fit')

subplot(2,1,2)
plot(data.T,data.Y)
hold on
% plot(data.T,p_Airy(beta0,data.T),'k-.')
% plot(data.T,p_Airy(fit_data.Coefficients.Estimate,data.T),'r-.')
set(gca,'Yscale','log')
title('Single scan fit attempt')
xlabel('PZ voltage')
ylabel('log(Voltage)')
legend('Raw data','Fit seed','Fit')

sfigure(300);
clf;
subplot(2,1,1)
plot(FFT(1,:),abs(FFT(2,:)))
title('Scan FFT')
xlabel('Freq [/V?]')
ylabel('log Power')
subplot(2,1,2)
plot(FFT(1,:),abs(FFT(2,:)))
set(gca,'Yscale','log')
title('Scan log FFT')
xlabel('Freq [/V?]')
ylabel('Power')
end

function Y_air = p_Airy(p,x)
    Y_air = p(1)*Airy(p(2:end-1),x)+p(end);
end