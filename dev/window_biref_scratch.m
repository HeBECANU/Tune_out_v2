%% For centred window
P0 = [3.03e-6,1.68e-3];
R0 = P0(1)/sum(P0); 
T1  = [20, 32, 40, 50, 60];
P1  = [2.7e-6,1.32e-3;
        2.1e-6,1.47e-3;
        2.5e-6,1.56e-3;
        1.9e-6,1.63e-3;
        2.2e-6,1.56e-3];
R1 = P1(:,1)./sum(P1,2);
<<<<<<< HEAD
plot(T1,R1,'ko')
hold on
plot([20,60],[R0,R0],'k')
=======

>>>>>>> 3ef538a7a011b2dc14aa45d1d07b10a64f06c902
% hold off

%% For offset window
P02 = [2.75e-6,1.427e-3];
R02 = P02(1)/sum(P02); 
T2  = [20, 45, 40, 35, 31];
P2  = [30.9e-6,1.32e-3;
        2.3e-6,1.53e-3;
        7.7e-6,1.54e-3;
        13.4e-6,1.57e-3;
        20.3e-6,1.5e-3];
R2 = P2(:,1)./sum(P2,2);
<<<<<<< HEAD
plot(T2,R2,'ro')
plot([20,60],[R0,R0],'r')
hold off


xlabel('Temperature')
set(gca,'Yscale','log')
ylabel('Output power ratio')
legend('Centred window','No window', 'Offset window','No window')
title('Temperature dependence of window birefringence')
=======

subplot(1,2,1)
plot(T1,R1,'ko')
hold on
plot([20,60],[R0,R0],'k')
plot(T2,R2,'ro')
plot([20,60],[R0,R0],'r')
hold off
xlabel('Temperature [C]')
ylabel('Output power ratio')
legend('Centred window','No window', 'Offset window','No window')
title('Lin scale')

subplot(1,2,2)
plot(T1,R1,'ko')
hold on
plot([20,60],[R0,R0],'k')
plot(T2,R2,'ro')
plot([20,60],[R0,R0],'r')
hold off
set(gca,'Yscale','log')
xlabel('Temperature [C]')
ylabel('Output power ratio')
title('Log scale')

suptitle('Temperature dependence of window birefringence')
>>>>>>> 3ef538a7a011b2dc14aa45d1d07b10a64f06c902

