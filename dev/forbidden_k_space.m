% sim rf knife method for forbidden transition

num_scattered=1e4;
recoil_vel1=1;
recoil_vel2=0.6;
recoil_vel3=0.4;

%% make a spont. scatt. halo of radius recoil_vel1 centered at (1,0,0)*recoil_vel1
% this is for a single photon decay
% k_scattered=randn(num_scattered,3);
% k_scattered=k_scattered./repmat(vecnorm(k_scattered,2,2),1,3); 
% k_scattered=k_scattered.*recoil_vel1;
% k_scattered(:,1)=k_scattered(:,1)+recoil_vel1;
% scatter3(k_scattered(:,1),k_scattered(:,2),k_scattered(:,3),'k.')
% xlabel('x')
% ylabel('y')
% zlabel('z')

%% try to simulate a 2 photon decay
%absorbtion photon
k_photon1=zeros(num_scattered,3);
k_photon1(:,1)=recoil_vel1; % displace in the x dirn
k_photon2= randn(num_scattered,3); % random direction
k_photon2=k_photon2./repmat(vecnorm(k_photon2,2,2),1,3); 
k_photon2=k_photon2.*recoil_vel2;
k_photon3= randn(num_scattered,3); % random direction
k_photon3=k_photon3./repmat(vecnorm(k_photon3,2,2),1,3); 
k_photon3=k_photon3.*recoil_vel3;

% sum up the k space
k_scattered=sum(cat(3,k_photon1,k_photon2,k_photon3),3);

%% plot the k_space
stfig('k space dist. inital');
subplot(2,2,1)
scatter3(k_scattered(:,1),k_scattered(:,2),k_scattered(:,3),'k.')
xlabel('x')
ylabel('y')
zlabel('z')

%% find the distribution of the k space radi
dist_from_origin=vecnorm(k_scattered,2,2);
in_sh.xdat=dist_from_origin;
in_sh.min=0;
in_sh.max=sum([recoil_vel1,recoil_vel2,recoil_vel3]);
in_sh.sigma=1e-4;
out=smooth_hist(in_sh);
stfig('k space dist. inital');
subplot(2,2,2)
plot(out.bin.centers,out.count_rate.smooth)
xlabel('k norm')
ylabel('radial count density')
% should normalize by the shel size at each radial bin to give the count density
%% cdf
cum_density=cumsum(out.counts.smooth)./num_scattered;
stfig('k space dist. inital');
subplot(2,2,3)

plot(out.bin.centers,cum_density)
xlabel('radial distance from origin')
ylabel('cumulative count fraction')



%% try to find the radial distribution about the center
meank=mean(k_scattered);
dist_from_mean=vecnorm(k_scattered-repmat(meank,size(k_scattered,1),1),2,2);
in_sh.xdat=dist_from_mean;
in_sh.min=0;
in_sh.max=1;
in_sh.sigma=1e-5;
out=smooth_hist(in_sh);
stfig('k space dist. inital');
subplot(2,2,4)
plot(out.bin.centers,out.count_rate.smooth)
xlabel('distance from center')
ylabel('radial count density')
% should normalize by the shel size at each radial bin to give the count density

%% use the rf kife to lense in momentum space
knife_vel=1;
k_scatt_norm=vecnorm(k_scattered,2,2);
greater_than_knife_mask=k_scatt_norm>knife_vel;
num_outcoupled=sum(greater_than_knife_mask);
fprintf('outcoupled fraction %f\n',num_outcoupled/num_scattered)

k_outcoupled=k_scattered(greater_than_knife_mask,:);
stfig('k space dist. after knife');
scatter3(k_outcoupled(:,1),k_outcoupled(:,2),k_outcoupled(:,3),'k.')
title('sudo k space')

%% remove the velocity that is given up to the trap
k_outcoupled_norm=vecnorm(k_outcoupled,2,2);
k_out_unit_vec=k_outcoupled./repmat(k_outcoupled_norm,1,3);
k_outcoupled=k_outcoupled-k_out_unit_vec*knife_vel;
stfig('k space dist. after knife');
scatter3(k_outcoupled(:,1),k_outcoupled(:,2),k_outcoupled(:,3),'k.')
title('out k space')

%% find the distirbution on the detector
%z(t)=x(0)+z'(0)*t-1/2 g t^2
%0=fall_dist+z'(0)*t-1/2 g t^2
% use the quadratic formula
%t=
% do it properly later
detector_distance=0.700; % put in actual number later
fall_time=k_outcoupled(:,3)*0+1;

xy_det=k_outcoupled(:,1:2).*repmat(fall_time,1,2);
stfig('detector dist. after knife');
subplot(2,2,1)
scatter(xy_det(:,1),xy_det(:,2),'.')

%%
radial_det_distance=vecnorm(xy_det,2,2);
det_radius=0.5;
num_det_counts=sum(radial_det_distance>det_radius);
fprintf('detected count fraction %f\n',num_det_counts/num_scattered)

%% find the distribution of the detector radi
in_sh.xdat=radial_det_distance;
in_sh.min=0;
in_sh.max=max(radial_det_distance);
in_sh.sigma=1e-4;
out=smooth_hist(in_sh);
stfig('detector dist. after knife');
subplot(2,2,2)
plot(out.bin.centers,out.count_rate.smooth)
xlabel('k norm')
ylabel('radial count density')
% should normalize by the shel size at each radial bin to give the count density
%% cdf
cum_density=cumsum(out.counts.smooth)./num_scattered;
stfig('detector dist. after knife');
subplot(2,2,3)

plot(out.bin.centers,cum_density)
xlabel('radial distance cen of det')
ylabel('cumulative count fraction')



