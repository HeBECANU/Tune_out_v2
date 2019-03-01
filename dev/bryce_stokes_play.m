fclose('all')
%add all subfolders
folder = fileparts(which(mfilename));
folder=strsplit(folder,filesep); %go up a directory
folder=strjoin(folder(1:end-1),filesep);
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

%https://doi.org/10.1364/AO.55.000B14

% define the stokes vector for some eliptical light
stokes_elip=@(phi,alpha) [1,cos(phi).*cos(alpha),cos(phi).*sin(alpha),sin(phi)]'; %stokes parameters

%% start with some eliptically polarized light
%lightin = stokes_lpvertical();
%lightin= mueller_stokes(mueller_rotator(rotate_angle),lightin)
rotate_angle=pi/8; %free parameter
light_in=stokes_elip(0.02*pi/2,rotate_angle); %first argument 0 is pure lin pi/2 is pure circ, second is roatation
% rotate it by some amount
%%
lin_pol_angles=linspace(0,2,1e3)*pi;
lin_pol=mueller_rotate(mueller_linpolarizer(),lin_pol_angles);
lightout=mueller_stokes(lin_pol,light_in);
ilightout = stokes_intensity(lightout);
%predict the min/max power from the light befor the polarizer
I_p_psi_chi=stokes_to_sph(light_in);
p_min=sin(I_p_psi_chi(4)).^2;
p_max=cos(I_p_psi_chi(4)).^2;
a_max=I_p_psi_chi(3);
a_min=I_p_psi_chi(3)+pi/2;
%%
figure(45)
set(gcf,'color','w')
clf
plot(lin_pol_angles/pi, ilightout);
title('transmitted intensity [should look like cos(2*a)^2');
xlabel('angle of analyzer BC');
ylabel('intensity [a.u.]');

xl=xlim;
yl=ylim;
hold on
line(xl,[1,1].*p_min,'color','r')
line(xl,[1,1].*p_max,'color','g')
line([1,1].*a_min/pi,yl,'color','r')
line([1,1].*a_max/pi,yl,'color','g')
hold off
legend('transmitted intensity','pred min','pred max');



%% ok we have now validated our tools lets compare the min/max power through the analyser for two elipticalites
rotate_angle=pi/8; %free parameter
impurity_angle=2*0.07;
lightin_pcirc=stokes_elip(impurity_angle,rotate_angle); %first argument 0 is pure lin pi/2 is pure circ, second is roatation
lightin_ncirc=stokes_elip(-impurity_angle,rotate_angle); %first argument 0 is pure lin pi/2 is pure circ, second is roatation
light_in={lightin_pcirc,lightin_ncirc};
I_p_psi_chi=stokes_to_sph(light_in);
p_min=sin(I_p_psi_chi(:,4)).^2;
p_max=cos(I_p_psi_chi(:,4)).^2;
a_max=I_p_psi_chi(:,3);
a_min=I_p_psi_chi(:,3)+pi/2;

%as expexted there is no measurabl difference

%%

rotate_angle_common=pi/2; %free parameter
rotate_angle_diff=-pi*7/8; %free parameter
impurity_angle=pi*0.001;
qwp_fast_angle_offset=0*(pi/180);
lightin_pcirc=stokes_elip(impurity_angle,rotate_angle_common); %first argument 0 is pure lin pi/2 is pure circ, second is roatation
lightin_ncirc=stokes_elip(-impurity_angle,rotate_angle_common+rotate_angle_diff); %first argument 0 is pure lin pi/2 is pure circ, second is roatation
light_lr_circ={lightin_pcirc,lightin_ncirc};
lr_compare=[];
for ii=1:2
light_in=light_lr_circ{ii};
I_p_psi_chi=stokes_to_sph(light_in);
a_max=I_p_psi_chi(:,3);
a_min_start=I_p_psi_chi(:,3)+pi/2;
%center the qwp rotation on the min or max transmission
qwp_angles=a_min_start+linspace(-1,1,1e4)*pi; %a_min_start

%qwp_angles=a_max+linspace(-1,1,1e3)*pi;
light_qwp=mueller_rotate(mueller_linretarder(2*pi/4),qwp_angles+qwp_fast_angle_offset);
lightout=mueller_stokes(light_qwp,light_in);

I_p_psi_chi=stokes_to_sph(lightout);
p_min=sin(I_p_psi_chi(:,4)).^2;
p_max=cos(I_p_psi_chi(:,4)).^2;
a_max=I_p_psi_chi(:,3);
a_min=I_p_psi_chi(:,3)+pi/2;

%asign to output struct
lr_compare(ii).p_min=p_min;
lr_compare(ii).p_max=p_max;
lr_compare(ii).a_max=a_max;
lr_compare(ii).a_min=a_min;
lr_compare(ii).qwp_angles=qwp_angles-a_min_start;
end
sfigure(46)
clf
set(gcf,'color','w')

angle_conv=360/(2*pi);
angle_conv_label='°';
subplot(2,2,1)
plot(lr_compare(1).qwp_angles*angle_conv,unwrap(lr_compare(1).a_max*4)*angle_conv/4,'r')
hold on
plot(lr_compare(2).qwp_angles*angle_conv,unwrap(lr_compare(2).a_max*4)*angle_conv/4,'b')
hold off
ylabel(sprintf('angle of PBS to give max power %s',angle_conv_label))
xlabel(sprintf('qwp angle %s',angle_conv_label))
subplot(2,2,2)
plot(lr_compare(1).qwp_angles*angle_conv,unwrap(lr_compare(1).a_min*4)*angle_conv/4,'r')
hold on
plot(lr_compare(2).qwp_angles*angle_conv,unwrap(lr_compare(2).a_min*4)*angle_conv/4,'b')
hold off
ylabel(sprintf('angle of PBS to give min power %s',angle_conv_label))
xlabel(sprintf('qwp angle %s',angle_conv_label))
subplot(2,2,3)
plot(lr_compare(1).qwp_angles*angle_conv,lr_compare(1).p_max,'r')
hold on
plot(lr_compare(2).qwp_angles*angle_conv,lr_compare(2).p_max,'b')
hold off
ylabel('max power of PBS')
xlabel(sprintf('qwp angle %s',angle_conv_label))
subplot(2,2,4)
plot(lr_compare(1).qwp_angles*angle_conv,lr_compare(1).p_min,'r')
hold on
plot(lr_compare(2).qwp_angles*angle_conv,lr_compare(2).p_min,'b')
hold off
ylabel('min power of PBS')
xlabel(sprintf('qwp angle %s',angle_conv_label))



%%
%lightin_pcirc=stokes_elip(impurity_angle,rotate_angle);
% qwp_angle=0.2
wp_angles = linspace(0,2*pi,100);
fdlts = zeros(size(wp_angles));
for ii=1:length(wp_angles)
    fdlts(ii) = trace(mueller_rotate(mueller_linretarder(2*pi/4),0)*mueller_rotate(mueller_linretarder(2*pi/4),wp_angles(ii)))/4;
end
sfigure(889);
plot(wp_angles,fdlts)