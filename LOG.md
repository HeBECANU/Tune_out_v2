## Scanning cavity analysis

Analysis of light-on scans is just about done. However, background calibration is needed so let's get to that.

Done! But the peak ratios are really bonkers, and the undersampling is a killer.

Detour by the scanning driver [manual](https://www.thorlabs.com/drawings/94a934db0e39d559-89295696-026D-279A-77612E4304F2F0D2/SA201-EC-Manual.pdf):


# Fabry Perot theory
From wiki.

* $t_{RT} = 1/\Delta\nu_{FSR}$ - the FSR is the inverse of the return time $2l/c$$!
* Supported resonances have frequencies $q\Delta\nu_{FSR}$ where $q$ is an integer. These are the frequencies which have the same phase after one round trip $2\phi(\nu) = 2\pi\nut_{RT}$.
* Fourier transform the decaying, oscillating electric field to find the spectral field density
$$
E_{0}\frac{1}{(2\tau_c)^{-1} + 2\pi i (\nu-\nu_q)}
$$
where $\tau_c$ is the field decay time; the line profile is thus
$$
\frac{1}{\tau_c}\frac{1}{(2\tau_c)^{-2}+4\pi^2(\nu-\nu_q)^2}
$$
* In terms of the Lorentzian linewidth $\Delta\nu_c = \frac{1}{2\pi\tau_c}$;
$$
\gamma_q(\nu) = \frac{1}{\pi}\frac{\Delta\nu_c/2}{(\Delta\nu_c/2)^{2}+4\pi^2(\nu-\nu_q)^2}
$$
* The (Lorentzian) Finesse is $\mathcal{F}_c = \frac{\Delta\nu_{FSR}}{\Delta\nu_c} = \frac{2\pi}{-\textrm{ln}(R_1 R_2)}$, where $\Delta\nu_c$ is the Lorentzian linewidth,
* When the cavity is used as a scanning interferometer, one resolves a linear combination of several Airy distributions in the output.
* In wavelength terms, the phase difference between successive transmitted photons is $\delta = (2\pi/\lambda)2nl\cos\theta$, with $\theta$ the angle at which the light passes through the etalon. 
