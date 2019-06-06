% sim rf knife method for forbidden transition

num_scattered=1e5;
recoil_vel1=1;
k_scattered=randn(num_scattered,3);

%%
k_scattered=k_scattered./repmat(vecnorm(k_scattered,2,2),1,3); 
k_scattered=k_scattered.*recoil_vel1;
k_scattered(:,1)=k_scattered(:,1)+recoil_vel1;
scatter3(k_scattered(:,1),k_scattered(:,2),k_scattered(:,3),'k.')
xlabel('x')
ylabel('y')
zlabel('z')

%% use the rf kife to lense in momentum space
knife_vel=1;
k_scatt_norm=vecnorm(k_scattered,2,2);
greater_than_knife_mask=k_scatt_norm>knife_vel;
k_outcoupled=k_scattered(greater_than_knife_mask,:);
scatter3(k_outcoupled(:,1),k_outcoupled(:,2),k_outcoupled(:,3),'k.')


%% remove the velocity that is given up to the trap
k_outcoupled_norm=vecnorm(k_outcoupled,2,2);
k_out_unit_vec=k_outcoupled./repmat(k_outcoupled_norm,1,3);
k_outcoupled=k_outcoupled-k_out_unit_vec*knife_vel;
scatter3(k_outcoupled(:,1),k_outcoupled(:,2),k_outcoupled(:,3),'k.')


