%% the idea of this function is to create a bunch of differenet models of the pol state for hwp/ qwp data angles
pol_opts.location = 'post';
if strcmp(pol_opts.location,'pre')

%Pre (at) window pol data for QWP data
%QWP HWP QWP_anal HWP_Anal :155 Phi_min Phi_max (Transmit,Reflect)...
data=  [130,333,nan,155,343,140,332,131,112,108,347,196;
        134,333,nan,155,300,162,325,160,242,130,352,286;
        138,333,nan,155,273,185,275,180,068,150,314,109;
        142,333,nan,155,252,202,256,202,169,178,285,209;
        146,333,nan,155,245,216,270,185,270,156,307,124;
        150,333,nan,155,242,222,297,165,297,138,335,130;
        154,333,nan,155,235,256,325,138,079,116,380,117;
        162,333,nan,155,195,278,380,66.,039,084,409,084;
        177,333,nan,155,087,282,527,10.2,122,12.,421,165;
        187,333,nan,155,022,406,500,03.6,204,0.8,437,250;
        202,333,nan,155,21.7,412,420,38.8,287,27.8,418,242;
        220,333,nan,155,150,306,350,129,275,125,367,232;
        226,333,nan,155,199,262,309,157,179,159,352,126;
        234,333,nan,155,240.,210,260,208,172,220,282,208;
        246,333,nan,155,242.,211,353,158,226,147,346,272;
        254,333,nan,155,267.,222,411,111,037,100,395,351;
        280,333,nan,155,021.,477,535,03.7,204,01.3,508,249;
        283,333,nan,155,13.7,464,540,05.2,204,04.6,463,158;
        270,340,nan,155,126.,360,412,105,121,076,368,77;
        274,350,nan,155,203.,280,317,207,024,184,291,86;
        268,354,nan,155,299.,153,297,155,066,160,312,21;
        280,010,nan,155,472,045,465,41.7,066,34,437,22;
        290,006,nan,155,326,158,372,148,059,140,366,14];
   
%% complete polarisation data for post window
elseif strcmp(pol_opts.location,'post')
load('C:\Users\BEC Machine\Documents\K\Tune_out_v2_trap_freq\dev\post_pol_1.mat') %load in the data
pol_data_val = cell2mat(pol_data(:,1));
%options measure (what we measured) or interpolated (create model based on all measured values)
pol_opts.predict = 'interp';
    
%if we want to use the observation method
V = pol_data_val(:,7).*2.*sqrt(pol_data_val(:,3).*pol_data_val(:,5))...
    ./(pol_data_val(:,3)+pol_data_val(:,5));%the V parameter for each run
d_p = (pol_data_val(:,3)-pol_data_val(:,5))./(pol_data_val(:,3)+pol_data_val(:,5));%contrast
theta = pol_data_val(:,6).*pi/180; %using min pow angle

%this next bit is tricky, have to spit it up into the qwp and hwp sections
if strcmp(pol_opts.predict,'interp') %the interpolation method
    qwp_ang = pol_data_val(38:63);
    sin_mdl = @(b,x) b(1).*sin(b(2).*x.*pi/180+b(3))+b(4);
    sin_mdl_abs = @(b,x) abs(b(1).*sin(b(2).*x.*pi/180+b(3))+b(4));
    beta0 = [1.0,2.0,0.5,0.0];
    fit_V_qwp = fitnlm(qwp_ang,V(38:63),sin_mdl,beta0);
    fit_d_p_qwp = fitnlm(qwp_ang,d_p(38:63),sin_mdl_abs,beta0);
    
    V(38,63) = sin_mdl(fit_V_qwp.Coefficients.Estimate,qwp_ang);
    d_p(38,63) = real(1-sin_mdl(fit_V_qwp.Coefficients.Estimate,qwp_ang).^2);
    
%     ang_mdl = @(b,x) 0.5.*atan((cos(b(1)-2.*x)+cot(2.*x).*tan(b(2)))...
%         ./(cos(b(1)-2.*x).*cot(2.*x)-tan(b(2))))+b(3);
%     theta = mod(theta(:)+pi/4,pi);
%     beta0 = [5.2,0,0.5];
%     fit_theta = fitnlm(qwp_ang.*pi/180,theta(38:63),ang_mdl,beta0);
%     sfigure(329839282);
%     clf
%     scatter(qwp_ang,theta(38:63)) %interpolate angles bit trickery
%     hold on
%     xsamp = linspace(0,2*pi);
%     plot(xsamp.*180/pi,ang_mdl(fit_theta.Coefficients.Estimate,xsamp'))
%     xlabel('qwp')
%     ylabel('min angle (rads)')
end
end