% sample code using Mueller-Stokes calculus
%
%  File information:
%     version 1.0 (jan 2014)
%     (c) Martin Vogel
%     email: matlab@martin-vogel.info
%
%  Revision history:
%     1.0 (jan 2014) initial release version

%% -------------------------------------------------------------------
% example 1: send light with horizontal linear polarization through 
% a rotating, perfect halfwave plate and subsequent polarizer: 
% final intensity should vary as cos(2*angle)^2.

angles = 0:360;
wps = mueller_rotate(mueller_waveplate(0.5, 'wav'), angles, 'deg');
lightin = stokes_lphorizontal();
lightout = mueller_stokes(mueller_linpolarizer(),wps,lightin);
ilightout = stokes_intensity(lightout);
figure();
plot(angles, ilightout);
title('transmitted intensity [should look like cos(2*a)^2');
xlabel('angle of halfwave plate axis');
ylabel('intensity [a.u.]');
legend('transmitted intensity');

%% -----------------------------------------------------------------
% example 2: send light with horizontal linear polarization through 
% a rotating, non-perfect halfwave plate and subsequent polarizer: 
% final intensity should deviate from the perfect cos(2*angle)^2
% curve, never reaching zero transmission

angles = 0:360;
wps = mueller_rotate(mueller_waveplate(0.5, 'wav'), angles, 'deg');
wps2 = mueller_rotate(mueller_waveplate(0.4, 'wav'), angles, 'deg');
lightin = stokes_lphorizontal();
lightout = mueller_stokes(mueller_linpolarizer(),wps,lightin);
ilightout = stokes_intensity(lightout);
lightout2 = mueller_stokes(mueller_linpolarizer(),wps2,lightin);
ilightout2 = stokes_intensity(lightout2);
figure();
plot(angles, ilightout, angles, ilightout2);
title('transmitted intensity with perfect and non-perfect halfwave plate');
xlabel('angle of halfwave plate axis');
ylabel('intensity [a.u.]');
legend('perfect (0.5-)plate', 'non-perfect (0.4-)plate');

%% -----------------------------------------------------------------
% example 3: send light with horizontal linear polarization through 
% rotating waveplates with increasing delay and subsequent polarizer 

angle = 0:360;
delay = 0:0.05:1;
% angles are in rows, delays in columns
angle_all = repmat(angle, [length(delay), 1]);
delay_all = repmat(delay', [1, length(angle)]);
wps3 = mueller_waveplate(delay_all, 'wav');
wps3 = mueller_rotate(wps3, angle_all, 'deg');
lightin = stokes_lphorizontal();
lightout3 = mueller_stokes(mueller_linpolarizer(),wps3,lightin);
ilightout3 = stokes_intensity(lightout3);
figure();
plot(angle, ilightout3);
title('transmitted intensity with plates of increasing delay');
xlabel('angle of plate axis');
ylabel('intensity [a.u.]');
legend(cellfun(@(x)sprintf('delay=%.2f',x),num2cell(delay),'UniformOutput',false));
