% first we need to calulate the derivative of the tune out with a change in
% each of the constituent transtion frequencies

% find tune-out using crossing search
fminopt = optimset('TolX',1e-20,'TolFun',1e-6);   %,'Display','iter'
tune_out_freq=fzero(@(x) aprox_he_polz(x),f2wl(413e-9),fminopt);  
%%

to_polz_si_deriv=derivest(@(x) outselect(2,@aprox_he_polz,x),tune_out_freq,'DerivativeOrder',1,'FixedStep',0.1)

%%

fminopt = optimset('TolX',1e-20,'TolFun',1e-6);   %,'Display','iter'
trans_dat=metastable_helium_transition_data
num_trans=numel(trans_dat);
deriv_trans=nan(num_trans,1)
for ii=1:num_trans
    fun_find_to_with_shift=@(x) fzero(shift_transition_he_polz(x,ii),tune_out_freq,fminopt);
    [tune_out_trans_sens_val,tune_out_trans_sens_unc]=derivest(fun_find_to_with_shift,1e3,'vectorized','no','FixedStep',1e3);
    deriv_trans(ii)=tune_out_trans_sens_val;
end

%%
sprintf('the sensitivity of the tune out to a change in the transition frequencies is\m')
labels=cellfun(@(x) x.name,trans_dat,'UniformOutput',false);
cat(2,num2cell(deriv_trans),labels)


%% check the derivative by hand
mf_trans_shift=1e3;
to_mf=fzero(shift_transition_he_polz(mf_trans_shift,1),tune_out_freq,fminopt)
mf_to_shift=to_mf-tune_out_freq;
rough_deriv=mf_to_shift/mf_trans_shift;



%% Calculate the wost case shift

% we will usea a 10mm/s oscillation because thats about what the tune out used
tf_param=[];
tf_param.trap_freq=2*pi*[55,411,415]; 
tf_param.mass=const.mhe;
tf_param.a_scat_len=const.ahe_scat;
tf_param.num=7e5;
tf_details=bec_properties(tf_param.trap_freq/(2*pi),tf_param.num,tf_param.mass,tf_param.a_scat_len);
tf_details.density_mean

(8*pi*const.hb/const.mhe)*tf_details.density_mean*(-const.ahe_scat+300e-9)%5000*const.a0
%%


%%
function fun_handle_out=shift_transition_he_polz(delt,idx)
transitions=metastable_helium_transition_data;
if idx>numel(transitions)
    error('idx not valid')
end
transitions{idx}.frequency_hz.val=transitions{idx}.frequency_hz.val+delt;
fun_handle_out=@(x) outselect(2,@aprox_he_polz,x,transitions);
end