function [to_val_out,detailed_out]=two_stage_two_dim_to_fit(full_data,fit_offset,use_fit_weights,fit_order)
%return the value in MHZ away from fit_offset as to_val_out
% return more details in detailed_out

full_data=cell2mat(full_data);

polz_cont=full_data(:,1);
polz_v=full_data(:,2);
polz_theta=full_data(:,3);
to_val_for_polz=full_data(:,4);
to_unc=full_data(:,5);


predictor = [polz_cont,polz_v,polz_theta];
to_val_fit = to_val_for_polz;

w_fit=1./(to_unc.^2);
if ~use_fit_weights
    w_fit=w_fit*0+1;
end
w_fit = w_fit./sum(w_fit);

% so here we will use a single peice of information from the atomic theory, namely the gradient of the reduced tensor
% term, from Li-Yan Tang's document "Response to Bryce's questions" \lambda_to (theta_p=0)>\lambda_to (theta_p=90)
% & as the polarizability decreased with increasing wavelength d \alpha^Scalar(\lambda)/ d \lambda <0 
% ? a increase in angle theta_p decreases the polarizability
% & as the prefactor in front of the tensor term decreases with increasing theta_p 
% ? \alpha^Tensor(\omega?\omega_TO)>0
% as \alpha^Scalar(\omega)/ d \omega >0
%  ? the reduced tensor term is positive
tune_out_vals_scaled=(to_val_fit-fit_offset)*1e-6;



gf_opt=[];
switch fit_order
    case 'scalar_sec'
        gf_opt.domain=[[-1,1]*1e5;...   %tune_out_scalar
            [-1,1]*1e12;...  %reduced_vector*cos_theta_k
            [-1,1]*1e12;...  %reduce_tensor
            [-1,1]*pi/2;...  %angle between polz measurment basis and B cross k
            pi+[-1,1]*pi/4;... % theta k
            [-1,1]*1e12;...
            ];
        gf_opt.start=[0, 1e5, 1e3, (rand(1)-0.5)*pi/4,pi,1e3];%pi+(rand(1)-0.5)*pi/4
        param_names={'tune_out_scalar','reduced_scalar','reduced_vec','phase','thetak','reduce_tensor'};
    case 'all_first'
         gf_opt.domain=[[-1,1]*1e5;...   %tune_out_scalar
            [-1,1]*1e6;...  %reduced_vector*cos_theta_k
            [-1,1]*1e6;...  %reduce_tensor
            [-1,1]*pi/2;...  %angle between polz measurment basis and B cross k
            pi+[-1,1]*pi/4;... % theta k
            [-1,1]*1e0;...
            [-1,1]*1e0
            ];
        gf_opt.start=[0, 1e3, 1e4, (rand(1)-0.5)*pi/4,pi,0,0];%pi+(rand(1)-0.5)*pi/4
        param_names={'tune_out_scalar','reduced_vector','reduce_tensor','phase','thetak','vec_deriv','tens_deriv'};
    case 2
        gf_opt.domain=[[-1,1]*1e5;...   %tune_out_scalar
            [-1,1]*1e6;...  %reduced_vector*cos_theta_k
            [0,1]*1e6;...  %reduce_tensor
            [-1,1]*pi/2;...  %angle between polz measurment basis and B cross k
            pi+[-1,1]*pi/4;... % theta k
            [-1,1]*1e5;...
            [-1,1]*1e5;...
            [-1,1]*1e5
            ];
        gf_opt.start=[0, 1e3, 1e4, (rand(1)-0.5)*pi/4,pi,1,0,0];%pi+(rand(1)-0.5)*pi/4
        param_names={'tune_out_scalar','reduced_vector','reduce_tensor','phase','thetak','scalar_deriv','vec_deriv','tens_deriv'};
    case 3
        gf_opt.domain=[[-1,1]*1e5;...   %tune_out_scalar
            [-1,1]*1e6;...  %reduced_vector*cos_theta_k
            [0,1]*1e6;...  %reduce_tensor
            [-1,1]*pi/2;...  %angle between polz measurment basis and B cross k
            pi+[-1,1]*pi/4;... % theta k
            [-1,1]*1e5;...
            [-1,1]*1e5;...
            [-1,1]*1e5;...
            [0.000000000001,1]*1e5
            ];
        gf_opt.start=[0, 1e3, 1e4, (rand(1)-0.5)*pi/4,pi,1,0,0,0.01];%pi+(rand(1)-0.5)*pi/4
        param_names={'tune_out_scalar','reduced_vector','reduce_tensor','phase','thetak','scalar_deriv','vec_deriv','tens_deriv','scalar_second_deriv'};
    otherwise
        gf_opt.domain=[[-1,1]*1e5;...   %tune_out_scalar
            [-1,1]*1e6;...  %reduced_vector*cos_theta_k
            [0,1]*1e6;...  %reduce_tensor
            [-1,1]*pi/2;...  %angle between polz measurment basis and B cross k
            pi+[-1,1]*pi/4;... % theta k
            ];
        gf_opt.start=[0, 1e3, 1e4, (rand(1)-0.5)*pi/4,pi];%pi+(rand(1)-0.5)*pi/4
        param_names={'tune_out_scalar','reduced_vector','reduce_tensor','phase','thetak'};
end
gf_opt.rmse_thresh=2e3;
gf_opt.plot=false;
gf_opt.level=2;
gf_out=global_fit(predictor,tune_out_vals_scaled,@full_tune_out_polz_model,gf_opt);

%
opts = statset('MaxIter',1e4);
beta0 = gf_out.params;%0.2 or 0.8
%beta0=[0, 1e3, 1*1e3, (rand(1)-0.5)*2*pi]
fit_mdl_full = fitnlm(predictor,tune_out_vals_scaled,@full_tune_out_polz_model,beta0,...
   'Weights',w_fit,...
    'CoefficientNames' ,param_names,...
    'Options',opts);

fit_vals_full = fit_mdl_full.Coefficients.Estimate;
fit_uncs_full = fit_mdl_full.Coefficients.SE;
%


% if fit_vals_full(3)<0
%     error('fit reduced tensor term is less than zero')
% end


%
[tsmht_details.val_predict,tsmht_details.unc_predict]=...
    predict(fit_mdl_full,[-1,0,-fit_vals_full(4)],'Alpha',1-erf(1/sqrt(2)),'Prediction','Curve');
tsmht_details.val_predict=tsmht_details.val_predict*1e6+fit_offset;
tsmht_details.unc_predict=1/2*range(tsmht_details.unc_predict)*1e6;

tsmht_details.val_no_offset_mhz=  predict(fit_mdl_full,[-1,0,-fit_vals_full(4)],'Alpha',1-erf(1/sqrt(2)),'Prediction','Curve');%fit_vals_full(1)+(1/2)*fit_vals_full(3)./(fit_vals_full(6)+(1/2)*fit_vals_full(8));
% fit_vals_full(1)+(1/2)*fit_vals_full(3);
tsmht_details.val=tsmht_details.val_no_offset_mhz*1e6+fit_offset;
tsmht_details.unc=sqrt(fit_uncs_full(1)^2 + (fit_uncs_full(3)/2)^2);
tsmht_details.unc=tsmht_details.unc*1e6;

detailed_out.tsmht=tsmht_details;
detailed_out.fit=fit_mdl_full;

if fit_vals_full(3)<0
    warning('fit has returned a reduced tensor term is less than zero')
    to_val_out=nan;
else
    to_val_out=tsmht_details.val_no_offset_mhz; %return this value in mhz for the bootstrap
end

end
