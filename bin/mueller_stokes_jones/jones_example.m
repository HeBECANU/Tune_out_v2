% sample code using Jones calculus
%
%  File information:
%     version 1.0 (jan 2014)
%     (c) Martin Vogel
%     email: matlab@martin-vogel.info
%
%  Revision history:
%     1.0 (jan 2014) initial release version

%% -----------------------------------------------------------------
% example 1: send light with horizontal linear polarization through 
% a rotating, perfect halfwave plate and subsequent polarizer: 
% final intensity should vary as cos(2*angle)^2.

angles = 0:360;
wps = jones_rotate(jones_waveplate(0.5, 'wav'), angles, 'deg');
lightin = jones_lphorizontal();
lightout = jones(jones_linpolarizer(),wps,lightin);
ilightout = jones_intensity(lightout);
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
wps = jones_rotate(jones_waveplate(0.5, 'wav'), angles, 'deg');
wps2 = jones_rotate(jones_waveplate(0.45, 'wav'), angles, 'deg');
lightin = jones_lphorizontal();
lightout = jones(jones_linpolarizer(),wps,lightin);
ilightout = jones_intensity(lightout);
lightout2 = jones(jones_linpolarizer(),wps2,lightin);
ilightout2 = jones_intensity(lightout2);
figure();
plot(angles, ilightout, angles, ilightout2);
title('transmitted intensity with perfect and non-perfect halfwave plate');
xlabel('angle of halfwave plate axis');
ylabel('intensity [a.u.]');
legend('perfect (0.5-)plate', 'non-perfect (0.45-)plate');

%% -----------------------------------------------------------------
% example 3: send light with horizontal linear polarization through 
% rotating waveplates with increasing delay  and subsequent polarizer

angle = 0:360;
delay = 0:0.05:1;
% angles are in rows, delays in columns
angle_all = repmat(angle, [length(delay), 1]);
delay_all = repmat(delay', [1, length(angle)]);
wps3 = jones_waveplate(delay_all, 'wav');
wps3 = jones_rotate(wps3, angle_all, 'deg');
lightin = jones_lphorizontal();
lightout3 = jones(jones_linpolarizer(),wps3,lightin);
ilightout3 = jones_intensity(lightout3);
figure();
plot(angle, ilightout3);
title('transmitted intensity with plates of increasing delay');
xlabel('angle of plate axis');
ylabel('intensity [a.u.]');
legend(cellfun(@(x)sprintf('delay=%.2f',x),num2cell(delay),'UniformOutput',false));

