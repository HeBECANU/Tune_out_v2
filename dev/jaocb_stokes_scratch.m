L_0 = [1, 1, 0, 0]';
% Let's partially depolarize the beam;
eta = .0;
% First we pass the beam through a half-wave plate at angle theta1
% For now we'll fix this at zero
theta_hwp =pi*0;
%Then we pass through a QWP at angle phi
theta_qwp = pi*0;
% And then through a linear polarizer at angle alpha
alpha = pi*0.;

DP = mueller_depolarizer(1-eta);
WP1 = mueller_rotate(mueller_waveplate(0.5),theta_hwp);
WP2 = mueller_rotate(mueller_waveplate(0.25),theta_qwp);
LP_t = mueller_rotate(mueller_linpolarizer(),alpha);
LP_r = mueller_rotate(mueller_linpolarizer(),alpha+pi/2);
% Which gives us two output beams;
L_trans = LP_t*WP2*WP1*DP*L_0;
L_rflct = LP_r*WP2*WP1*DP*L_0;

% Let's look at the transmitted intensities

sfigure(999);
clf;
plot(-1,L_trans(1),'ko')
hold on
plot(1,L_rflct(1),'bx')
hold on
xlim([-1.5,1.5])
ylim([0,2])

%% With input light
% Partially depolarized & circular
fwtext('Running')
th_circ = 90;
% L_0 = [1, cos(pi*th_circ/180), 0, sin(pi*th_circ/180)]';
L_0 = [1 0 1 0]';
eta = .00;
L_0 = mueller_depolarizer(1-eta)*L_0
% Our four measurements are:
I_00 = polz_tomo(0,0)*L_0;
I_045 = polz_tomo(0,45)*L_0;
I_090 = polz_tomo(0,90)*L_0;
I_4590 = polz_tomo(90,45)*L_0;
I = [I_00,I_045,I_090,I_4590]
S = [I_00 + I_090, I_00-I_090, 2*I_045-I_00-I_090, 2*I_4590-(I_00+I_090)];
S_out = S(1,:)'

% polz_tomo(0,45)


function M_tr = polz_tomo(theta_qwp,theta_plz)
    M_tr = mueller_rotate(mueller_linpolarizer(),deg2rad(theta_plz))*mueller_rotate(mueller_waveplate(0.25),deg2rad(theta_qwp));
end