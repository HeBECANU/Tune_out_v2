addpath('./lib/Core_BEC_Analysis/lib/') %add the path to set_up_project_path, this will change if Core_BEC_Analysis is included as a submodule
                  % in this case it should be 
set_up_project_path('./../..')

hebec_constants %call the constants function that makes some globals




%% here we will manualy build a struct containing all the previous tune out measurements

to_lit_rev=[];

n=1;
to_lit_rev(n).to_freq_val=725736700e6;
to_lit_rev(n).to_freq_unc=sqrt(40e6^2+260e6^2);

he_trans=metastable_helium_transition_data;
% this is roughly right
to_lit_rev(n).transition_1=f2wl(1083e-9);
% weight by the oscillator strength
to_lit_rev(n).transition_1=wmean(arrayfun(@(n) a{n}.frequency_hz.val,1:3),arrayfun(@(n) a{n}.osc_strength.val,1:3));
to_lit_rev(n).transition_2=f2wl(389e-9);
to_lit_rev(n).transition_2=wmean(arrayfun(@(n) a{n}.frequency_hz.val,4:6),arrayfun(@(n) a{n}.osc_strength.val,4:6));


to_lit_rev(n).name='this work';
to_lit_rev(n).url='none';

n=2;
[f,df]=f2wl(790.032388e-9,0.000032e-9);
to_lit_rev(n).to_freq_val=f;
to_lit_rev(n).to_freq_unc=df;
% d1 and d2 line freq from
% https://steck.us/alkalidata/rubidium87numbers.1.6.pdf
to_lit_rev(n).transition_1=377.1074635e12;
to_lit_rev(n).transition_2=384.2304844685e12;
to_lit_rev(n).name='High-precision measurements of the 87Rb D-line tune-out wavelength';
to_lit_rev(n).url='https://doi.org/10.1103/PhysRevA.92.052501';


n=3;  
[f,df]=f2wl(768.9712e-9,0.0015e-9);
to_lit_rev(n).to_freq_val=f;
to_lit_rev(n).to_freq_unc=df;
% d1 and d2 line freq from
% https://tobiastiecke.nl/archive/PotassiumProperties.pdf
to_lit_rev(n).transition_1=389.286058716e12;
to_lit_rev(n).transition_2=391.01617003e12;
to_lit_rev(n).name='Measurement of a Wavelength of Light for Which the Energy Shift for an Atom Vanishes';
to_lit_rev(n).url='https://doi.org/10.1103/PhysRevLett.109.243004';

n=4;  
[f,df]=f2wl(790.01858e-9,0.00023e-9);
to_lit_rev(n).to_freq_val=f;
to_lit_rev(n).to_freq_unc=df;
% d1 and d2 line freq from
% https://steck.us/alkalidata/rubidium87numbers.1.6.pdf
to_lit_rev(n).transition_1=377.1074635e12;
to_lit_rev(n).transition_2=384.2304844685e12;
to_lit_rev(n).name='Precision measurement of the 87Rb tune-out wavelength in the hyperfine ground state F=1 at 790 nm';
to_lit_rev(n).url='https://doi.org/10.1103/PhysRevA.93.022507';

%n=5
%https://doi-org.virtual.anu.edu.au/10.1103/PhysRevLett.109.243003
% Precision Measurement of Transition Matrix Elements via Light Shift Cancellation
% only get to a ratio frac uncert of 0.002/1.617


%%

for ii=1:size(to_lit_rev,2)
    f_to_val=to_lit_rev(ii).to_freq_val;
    idiot_check_freq(f_to_val)
    f_to_unc=to_lit_rev(ii).to_freq_unc;
    f_trans_1=to_lit_rev(ii).transition_1;
    idiot_check_freq(f_trans_1)
    f_trans_2=to_lit_rev(ii).transition_2;
    idiot_check_freq(f_trans_2)
    
    if ~( f_to_val>min(f_trans_1,f_trans_2) && f_to_val<max(f_trans_1,f_trans_2) )
        error('tune out is not in the middle')
    end
    
    freq_interval=abs(f_trans_1-f_trans_2);
    to_lit_rev(ii).transtion_interval=freq_interval;
    
    to_lit_rev(ii).frac_transtion_interval=abs(frac_diff(f_trans_1,f_trans_2));
    
    E1=const.h*f_trans_1;
    E2=const.h*f_trans_2;
    
    % calulate the ratio of oscillator strengths
    x_val= (-E2^2 + (const.h^2)*(f_to_val^2))/(E1^2 - (const.h^2)*(f_to_val^2));
    
    pd_f_to_with_x=  ((E1-E2)* (E1+E2)) /... 
                    (2 *const.h * (x_val+1)^(3/2) * sqrt(E1^2*x_val+E2^2));
    pd_f_to_with_x=abs(pd_f_to_with_x) ;
    
    to_lit_rev(ii).pd_f_to_with_x=pd_f_to_with_x;
    
    x_unc=(1/pd_f_to_with_x)*f_to_unc;
    frac_x_unc_1= x_unc/x_val ;
    
    frac_x_unc_2=( 2 *(const.h^2) * f_to_val * (E1^2-E2^2)) /...
    ( (E1^2-(const.h^2)*f_to_val^2) * (-E2^2 + (const.h^2)* f_to_val^2 ) )*...
    f_to_unc;
    frac_x_unc_2=abs(frac_x_unc_2);
    
    if abs(frac_diff(frac_x_unc_2,frac_x_unc_1))>1e-3
        error('numeical problem')
    end
    to_lit_rev(ii).x_val=x_val;
    to_lit_rev(ii).frac_x_unc=frac_x_unc_1;
    
    
    % now compute the same for R, the ratio of matrix elements
    
    r_val= (f_trans_1*(f_to_val^2-f_trans_2^2) ) /...
        ( f_trans_2 * (f_trans_1^2-f_to_val^2));
    
    pd_f_to_with_r =  ( (f_trans_1-f_trans_2)*(f_trans_1+f_trans_2))/...
                      ( 2 * sqrt( (1/f_trans_1) +(r_val/f_trans_2)) ...
                            * ((f_trans_1+r_val*f_trans_2)^(3/2)) );
    pd_f_to_with_r=abs(pd_f_to_with_r);
    r_unc=(1/pd_f_to_with_r)*f_to_unc;
    frac_r_unc_1= r_unc/r_val ;
    
    
    
    frac_r_unc_2 = 2*(f_trans_1^2 - f_trans_2^2) * f_to_val* f_to_unc /...
                    ( (f_trans_1^2-f_to_val^2)*(-f_trans_2^2+f_to_val^2) );
    frac_r_unc_2=abs(frac_r_unc_2);
    
    if abs(frac_diff(frac_r_unc_1,frac_r_unc_2))>1e-3
        error('numeical problem')
    end
    
    to_lit_rev(ii).r_val=r_val;
    to_lit_rev(ii).frac_r_unc=frac_r_unc_1;
    
end

%% mitroy frac unc calc

frac_y= (to_lit_rev(1).to_freq_unc/to_lit_rev(1).to_freq_val) * (1/(2*0.0582))


%%
function idiot_check_freq(freq)
    if freq>f2wl(100e-9) || freq<f2wl(2000e-9)
        error('this is an insane freq')
    end
end


