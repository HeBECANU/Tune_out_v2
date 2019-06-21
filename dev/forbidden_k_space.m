% sim rf knife method for forbidden transition

% find this .m file's path, this must be in the project root dir
this_folder = fileparts(fileparts(which(mfilename)));
% Add that folder plus all subfolders to the path.
addpath(genpath(this_folder));%add all subfolders to the path to find genpath_exclude
path_to_genpath=fileparts(which('genpath_exclude'));
path(pathdef) %clean up the path back to the default state to remove all the .git that were added
addpath(this_folder)
addpath(path_to_genpath)
addpath(genpath_exclude(fullfile(this_folder,'lib'),'\.')) %dont add hidden folders
addpath(genpath_exclude(fullfile(this_folder,'dev'),'\.'))
addpath(genpath_exclude(fullfile(this_folder,'bin'),'\.'))



%%
num_scattered=3e6;


f23s1_33s1=f2wl(427.7e-9);
%1803 p2 transition
%https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.119.263002
f23s1_23p2=276732186526.2e3;
f23p2_33s1=f23s1_33s1-f23s1_23p2;
recoil_vel1=freq_to_recoil_vel_he(f23s1_33s1); %absorb 427nm photon
recoil_vel2=freq_to_recoil_vel_he(f23p2_33s1); %emit 706nm photon
recoil_vel3=freq_to_recoil_vel_he(f23s1_23p2);%emit 1083nm photon
knife_vel=recoil_vel1*1;
det_radius=70e-3/2;

%in the limit of a infinitesimal size detector then you want the recoil velocity at the peak of the radial count density
%vs k norm which is at about 1.2*recoil_vel1
% for our detector the optimal is ~1 recoil

%%


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
hold on
yl=ylim;
line([1,1]*knife_vel,yl,'color','k')
line([1,1]*recoil_vel1*2,yl,'color','r')

hold off
% should normalize by the shel size at each radial bin to give the count density
%% cdf
cum_density=cumsum(out.counts.raw)./num_scattered;
stfig('k space dist. inital');
subplot(2,2,3)
plot(out.bin.centers,cum_density)
xlabel('radial distance from origin')
ylabel('cumulative count fraction')
hold on
line([1,1]*knife_vel,[0,1.1],'color','k')
line([1,1]*recoil_vel1*2,[0,1.1],'color','r')
ylim([0,1.1])
hold off


%% try to find the radial distribution about the center
meank=mean(k_scattered);
dist_from_mean=vecnorm(k_scattered-repmat(meank,size(k_scattered,1),1),2,2);
in_sh.xdat=dist_from_mean;
in_sh.min=min(dist_from_mean);
in_sh.max=max(dist_from_mean);
in_sh.sigma=1e-5;
out=smooth_hist(in_sh);
stfig('k space dist. inital');
subplot(2,2,4)
plot(out.bin.centers,out.count_rate.smooth)
xlabel('distance from center')
ylabel('radial count density')
hold on
yl=ylim;
line([1,1]*(recoil_vel2-recoil_vel3),yl,'color','k')
line([1,1]*recoil_vel1,yl,'color','r')
hold off

% should normalize by the shel size at each radial bin to give the count density

%% use the rf kife to lense in momentum space
k_scatt_norm=vecnorm(k_scattered,2,2);
greater_than_knife_mask=k_scatt_norm>knife_vel;
num_outcoupled=sum(greater_than_knife_mask);
if num_outcoupled==0
    warning('no atoms outcoupled')
end
fprintf('outcoupled fraction %f\n',num_outcoupled/num_scattered)

k_outcoupled=k_scattered(greater_than_knife_mask,:);
stfig('k space dist. after knife');
subplot(2,2,1)
scatter3(k_outcoupled(:,1),k_outcoupled(:,2),k_outcoupled(:,3),'k.')
title('sudo k space')

%% remove the velocity that is given up to the trap
k_outcoupled_norm=vecnorm(k_outcoupled,2,2);
k_out_unit_vec=k_outcoupled./repmat(k_outcoupled_norm,1,3);
k_outcoupled=k_outcoupled-k_out_unit_vec*knife_vel;
stfig('k space dist. after knife');
subplot(2,2,2)
scatter3(k_outcoupled(:,1),k_outcoupled(:,2),k_outcoupled(:,3),'k.')
title('out k space')

%% find the distirbution on the detector
%z(t)=x(0)+z'(0)*t-1/2 g t^2
%0=fall_dist+z'(0)*t-1/2 g t^2
% use the quadratic formula
%t=(z'(0)+sqrt((z'(0))^2+2*g*fall_dist))/g
% do it properly later
hebec_constants
global const
detector_fall_distance=0.8517; % check this value
fall_time=(k_outcoupled(:,3)+sqrt((k_outcoupled(:,3)).^2+2*const.g0*detector_fall_distance))/const.g0;
xy_det=k_outcoupled(:,1:2).*repmat(fall_time,1,2);

%%
radial_det_distance=vecnorm(xy_det,2,2);
detected_mask=radial_det_distance<det_radius;
num_det_counts=sum(detected_mask);
fprintf('detected count fraction %f\n',num_det_counts/num_scattered)
%% plot the dectected and undetected counts
stfig('detector dist. after knife');
subplot(2,2,1)
scatter(xy_det(detected_mask,1),xy_det(detected_mask,2),'.k')
hold on
scatter(xy_det(~detected_mask,1),xy_det(~detected_mask,2),'.b')
hold off

%%
dyn_range_pow=0.5;
spatial_blur=1;
% x_edges=col_vec(linspace(min(min(xy_det(:,1))),max(xy_det(:,1)),1e2));
% y_edges=col_vec(linspace(min(min(xy_det(:,2))),max(xy_det(:,2)),1e2));
x_edges=col_vec(linspace(-det_radius,det_radius,1e2));
y_edges=col_vec(linspace(-det_radius,det_radius,1e2));
bin_area=diff(x_edges(1:2))*diff(y_edges(1:2));
[counts,centers]=hist3(xy_det,'edges',{x_edges,y_edges});
counts=counts/bin_area;

%imagesc seems to plot the wrong way round so we transpose here

if dyn_range_pow~=1
    counts=counts.^dyn_range_pow;
end
if  ~spatial_blur==0
    counts=imgaussfilt(counts,spatial_blur);
end

stfig('detector dist. after knife');
subplot(2,2,2)
imagesc(10^3*centers{1},10^3*centers{2},transpose(counts))
colormap(viridis)
set(gca,'Ydir','normal')
set(gcf,'Color',[1 1 1]);
title('Spatial Dist. TOP')
xlabel('X(mm)')
ylabel('Y(mm)')
c = colorbar;
c.Label.String=sprintf('count density^{%.2f}',dyn_range_pow);

        
%% find the distribution of the detector radi
in_sh.xdat=radial_det_distance;
in_sh.min=0;
in_sh.max=max(radial_det_distance);
in_sh.sigma=1e-4;
out=smooth_hist(in_sh);
stfig('detector dist. after knife');
subplot(2,2,3)
plot(out.bin.centers*1e3,out.count_rate.smooth)
xlabel('radial distance cen of det (mm)')
ylabel('radial count density')
hold on
yl=ylim;
line([1,1]*det_radius*1e3,yl,'color','k')
hold off
% should normalize by the shel size at each radial bin to give the count density
%% cdf
cum_density=cumsum(out.counts.smooth)./num_scattered;
stfig('detector dist. after knife');
subplot(2,2,4)
plot(out.bin.centers*1e3,cum_density)
xlabel('radial distance cen of det (mm)')
ylabel('cumulative count fraction')
hold on
yl=ylim;
line([1,1]*det_radius*1e3,yl,'color','k')
hold off


%% thermal leakage
%one consideration is how much thermal will leak out from having the rf knife on
% one d Maxwell–Boltzmann
% sqrt(m/(2*pi*k*T))*exp(-m v^2 / (2*k*T))
% which if we integrate fomr vmax to inf and -vmax to -inf gives
% see mathematica document thermal_leak_rate
% Erfc[(m vmax)/(Sqrt[2] Sqrt[k m T])]
% erfc(sqrt(const.mhe)*knife_vel/(sqrt(2*const.kb*1e-6)))

leak_rate=@(therm_temp,v_knife) erfc(sqrt(const.mhe)*v_knife/(sqrt(2*const.kb*therm_temp)));
leak_per_atom=leak_rate(1e-6,knife_vel)
total_leak=leak_per_atom*1e6



