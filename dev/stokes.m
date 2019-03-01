%phi = linspace(-pi/5,pi/5,2)';
phi = linspace(-pi/2*0.09,pi/2*0.09,2)'; %approx. the worst impurity that we have
alpha = 0;
S = [cos(phi).*cos(alpha),cos(phi).*sin(alpha),sin(phi)]'; %stokes parameters
p_min = [];
p_max = [];
theta_vec = linspace(-pi/2,pi/2,1e4);
S_path = zeros(3,length(phi),length(theta_vec)); %bryce:change: theta->theta_vec
S_lin = zeros(4,length(phi),length(theta_vec));
jj = 1;
theta_anal = pi/2;
rot_mat = [1, 0, 0, 0; 
    0, cos(2*theta_anal), sin(2*theta_anal), 0;
    0, -sin(2*theta_anal),cos(2*theta_anal),0;
    0, 0, 0, 1];
lin_pol = 0.5.*[1,1,0,0;1,1,0,0;0,0,0,0;0,0,0,0];
chi_p = [];
for theta = theta_vec
    QWP = [cos(2*theta)^2, sin(2*theta)*cos(2*theta), -sin(2*theta);
        sin(2*theta)*cos(2*theta), sin(2*theta)^2, cos(2*theta);
        sin(2*theta), -cos(2*theta), 0];
    HWP = [cos(2*theta)^2-sin(2*theta)^2, 2*cos(2*theta)*sin(2*theta),0;
        2*cos(2*theta)*sin(2*theta),-cos(2*theta)^2+sin(2*theta)^2,0;
        0,0,-1];
    %QWP-inv(QWP)
    Sp = HWP*S;
    S_path(:,:,jj) = Sp;
    S_lin(:,:,jj) = lin_pol*rot_mat*[1,1;Sp];
    chi = 1/2*atan(Sp(3,:)./sqrt(Sp(1,:).^2+Sp(2,:).^2));
    chi_p = [chi_p;chi];
    p_min = [p_min;sin(chi).^2];
    p_max = [p_max;cos(chi).^2];
    jj = jj + 1;
end
figure(1109);
clf
hold on
for ii = 1:length(phi)
    plot(theta_vec,p_min(:,ii))
end
xlabel('Angle of QWP from the vertical')
ylabel('Measured minimum normalised intensity')
grid on
figure(1108);
clf
hold on
for ii = 1:length(phi)
    plot(theta_vec,p_max(:,ii))
end
xlabel('Angle of QWP from the vertical')
ylabel('Measured max normalised intensity')
grid on
figure(1107);
clf
hold on
for ii = 1:length(phi)
    plot(theta_vec,p_min(:,ii)./p_max(:,ii))
end
xlabel('Angle of QWP from the vertical')
ylabel('Measured intensity ratio')
grid on
figure(1111);
clf
plot(theta_vec,squeeze(S_lin(1,1,:)))
hold on
plot(theta_vec,squeeze(S_lin(1,2,:)))
xlabel('Angle of QWP from the vertical')
ylabel('Transmitted power through Pol Analyser')

%% bryce playing with ploting
figure(1110);
[x y z] = sphere(256);
h = surfl(x, y, z);
colormap(viridis)
set(h, 'FaceAlpha', 0.5)
shading interp
hold on
temp = squeeze(S_path(:,2,:));
scatter3(temp(1,:),temp(2,:),temp(3,:))
scatter3(S(1,1),S(2,1),S(3,1),'bx')
temp = squeeze(S_path(:,1,:));

scatter3(temp(1,:),temp(2,:),temp(3,:))
axis equal
title('Poincare Sphere Representation')
xlabel('S_1')
ylabel('S_2')
zlabel('S_3')
hold off
sfigure(2211);
scatter(theta_vec,sin(2.*chi_p(:,1)))
xlabel('ang of QWP')
ylabel('A')
sfigure(2021);
scatter(theta_vec,1/2*atan(S_path(2,1,:)./S_path(1,1,:)))
hold on
scatter(theta_vec,1/2*atan(S_path(2,2,:)./S_path(1,2,:)))
hold off
xlabel('ang of QWP')
ylabel('\psi')