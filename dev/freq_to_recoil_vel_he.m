function atom_vel=freq_to_recoil_vel_he(optical_freq)
hebec_constants
global const 
photon_momentum=const.h*optical_freq/const.c;
atom_vel=photon_momentum/const.mhe;
end