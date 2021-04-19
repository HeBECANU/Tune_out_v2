%grav dechoherence calc
%consider a gravitational sensor consisting of an atom in the grounds tate of a harmonic trap
rdiff=1e-4;
rcent=1;
m1=const.mhe;
m_sens_atom=1.9e-29; %carbon
grav_const=6.67430e-11;
omega_trap=2*pi*1e12;

mass_of_sensor=6e24;%mass of earth%100e3;


uncert_force=omega_trap*const.hb/2
delt_force=rdiff*2*grav_const*m1*m_sens_atom/(rcent^3)
snr_single_sensor=delt_force/uncert_force
num_sensors=mass_of_sensor/m_sens_atom

snr_combined=snr_single_sensor*sqrt(num_sensors)

%%


uncert_force=omega_trap*const.hb/2
delt_force=rdiff*2*grav_const*m1*mass_of_sensor/(rcent^3)
snr_single_sensor=delt_force/uncert_force;


snr_combined=snr_single_sensor*sqrt(num_sensors)

%% lets attack this problem in a different way 
omega_year=2*pi/(360*24*60*60)

uncert_force=omega_year*const.hb/2
delt_force=rdiff*2*grav_const*m1*mass_of_sensor/(rcent^3)
snr_single_sensor=delt_force/uncert_force;


snr_combined=snr_single_sensor*sqrt(num_sensors)



%% eq 112 https://arxiv.org/pdf/1706.05677.pdf
mass_particle=const.me
a_c=(const.hb/grav_const)*((mass_particle^(-3)))
tau= mass_particle*a_c^2 /const.hb


tau=(const.hb/grav_const^2)*mass_particle^(-5)

