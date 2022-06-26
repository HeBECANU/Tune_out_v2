% tune out sensitivity calculations

%% first lets get the derivative in the polarizability

% find tune-out using crossing search
fminopt = optimset('TolX',1e-20,'TolFun',1e-6);   %,'Display','iter'
tune_out_freq=fzero(@(x) aprox_he_polz(x),f2wl(413e-9),fminopt);  
au_polz_at_to=aprox_he_polz(tune_out_freq);
f2wl(tune_out_freq)
hold on
plot(tune_out_freq*1e-12,au_polz_at_to,'xr')
hold off
fprintf('found tune out at %.6f THz, %.6f nm, polz %g au \n',tune_out_freq*1e-12,f2wl(tune_out_freq)*1e9,aprox_he_polz(tune_out_freq))

% find tune-out using minimization of abs value
% fminopt = optimset('TolX',1,'TolFun',1e-12);  
% tune_out_freq=fminsearch(@(x) abs(aprox_he_polz(x)),600e12,fminopt);  
% hold on
% plot(tune_out_freq*1e-12,aprox_he_polz(tune_out_freq),'xr')
% hold off
% fprintf('found tune out at %.6f THz, %.6f nm, polz %g au \n',tune_out_freq*1e-12,f2wl(tune_out_freq)*1e9,aprox_he_polz(tune_out_freq))
% 

%% derivatives in SI /hz
to_polz_si_deriv=arrayfun(@(y) derivest(@(x) outselect(2,@aprox_he_polz,x),tune_out_freq,'DerivativeOrder',y,'FixedStep',0.1),1:4);
to_polz_si_deriv(1)
to_polz_si_deriv(2)


%%
beam_waist=20e-6;
beam_power=0.1;%0.120;
peak_intensity= 2*beam_power/(pi*(beam_waist^2));
detuning_hz=1e9;
dipole_pot=- (1/(2*const.epsilon0*const.c))*to_polz_si_deriv(1)*peak_intensity*detuning_hz;

dipole_pot/const.kb
dipole_pot

%% calculate the shift from the hyperpolarizability
atom_u=[];
% atom_u.energy=(const.hb)^2/(const.me*(const.a0^2)); 
% atom_u.time=const.hb/atom_u.energy;
atom_u.energy=4.3597447222071e-18 ;%J
atom_u.time=2.4188843265857e-17 ;%s
atom_u.polz=1.64877727436e-41;% C2⋅m2⋅J−1 
atom_u.second_hyp_polz=6.2353799905e-65;% C^4⋅m^4⋅J^-3
atom_u.a0=5.29177210903e-11 ;% m


hyp_pol_au=-1.1e7;% from li-yan personal coresp
hyp_pol_si=hyp_pol_au*atom_u.second_hyp_polz;


freq_shif_w_power=-(2/24)*hyp_pol_si*(2/(const.c*const.epsilon0))*(1/to_polz_si_deriv(1))

freq_shif_w_power*peak_intensity



%% calculate the intensity from the gradient of the probe beam signal
probe_freq_sq_grad=80e-9;
% find the derivative of the probe beam trap frequency as a function of optical freqency, leave off the factor of probe
% intensity
%power_dep=( 1/(4*(pi^2)) )*(1/(2*const.epsilon0*const.c*const.mhe))...
%                    *to_polz_si_deriv(1)...
%                    *( (4)/(beam_waist^2) ) 

%to_polz_si_deriv(1)*
sens_param_A=4*beam_power/(const.mhe*const.epsilon0*const.c*pi*beam_waist^4);
                
%power_dep=(1/(4*(pi^2))) * (1/(2*const.mhe*const.c*const.epsilon0))  *to_polz_si_deriv(1)   *  ( 4/(beam_waist^2) );
fprintf('probe trap freq squared at 4 GHz detuning =%g',sens_param_A*4e9*to_polz_si_deriv(1)/(4*pi^2))
probe_freq_sq_grad/( sens_param_A)

(sens_param_A*to_polz_si_deriv(1)*1e9)

%*2*pi
single_shot_sens=(0.020*2*pi)*2*(420*2*pi)*(1/(sens_param_A*to_polz_si_deriv(1)));
fprintf('signle shot sensitivity %.3f GHz \n',1e-9*single_shot_sens)


flat_pot_detuning=1e-12*((420*2*pi)^2)/(sens_param_A*to_polz_si_deriv(1))
fprintf('flat potential limit %.3f THz \n',flat_pot_detuning)



%%

% %% trap freq for abbas dipole trap
% beam_waist=75e-6;
% beam_power=0.1;
% peak_intensity= 2*beam_power/(pi*(beam_waist^2));
% 
% trap_freq=( 1/(2*pi) )*sqrt((1/(2*const.epsilon0*const.c*const.mhe))...
%                     *outselect(2,@aprox_he_polz,f2wl(1550e-9))...
%                     *( (4)/(beam_waist^2) ) ...
%                     *  peak_intensity) 
                
