%phi = linspace(-pi/2.1,pi/2.1,2)';
phi = linspace(-pi/2*0.0444,pi/2*0.05,2)'; %approx. the worst impurity that we have
S = [cos(phi),zeros(length(phi),1),sin(phi)]'; %stokes parameters
p_min = [];
p_max = [];
theta_vec = linspace(-pi/2,pi/2,1e4);
S_path = zeros(3,length(phi),length(theta_vec)); %bryce:change: theta->theta_vec
jj = 1;
for theta = theta_vec
    QWP = [cos(2*theta)^2, sin(2*theta)*cos(2*theta), -sin(2*theta);
        sin(2*theta)*cos(2*theta), sin(2*theta)^2, cos(2*theta);
        sin(2*theta), -cos(2*theta), 0];
    %QWP-inv(QWP)
    Sp = QWP*S;
    S_path(:,:,jj) = Sp;
    chi = 1/2*atan(Sp(3,:)./abs(Sp(1,:)));
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