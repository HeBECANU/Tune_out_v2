## Article title:
Potassium Tune-out-wavelength measurement using atom interferometry and a multipass optical cavity

### Authors:
Raisa Trubko, Maxwell D. Gregoire, William F. Holmgren, Alexander D. Cronin

#### Reference:
[Physical Review A **95**, 052507 (2017)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.95.052507)

## Purpose

The authors measure a tune-out wavelength in Potassium (surprise). I am reading this paper to understand how they characterize the systematic errors relevant to the laser light.

## Context

Tuneout wavelengths are constrained by the global structure of atoms, and so their measurements are tied to properties like oscillator strengths & ratios thereof, lifetimes, polarizabilities (scalar, vector & tensor), hyperpolarizabilities, electron correlations, and corrections from relativistic and QED effects. They are important for species- and state-specific optical traps.

## Approach

The authors illuminate twin atom beams with multiple passes of a tunable laser in a multipass optical cavity. The laser waist is chosen to create a uniform intensity gradient across the beam paths. The atom beams are split and combined by diffraction gratings, but their SNR is terrible with large beam separations so *both* beams are illuminated by the laser. They compensate for this by using the multipass cavity, which enhances the phase shift between the atom beams as they pass through the cavity. The phase shift is measured, I assume because it's not stated, by  the change in flux at the atom detector which monitors one output port of the interferometer. They set the probe beam to only two values either side of the TO and perform a linear fit between them.

## Contribution

The authors measure the longest TO in potassium. They calculate the ratios of line strengths and oscillator strengths for - [ ] he $D2/D1$ lines, correponding to the $4s-4p_{1/2}$ and $4s-4p_{3/2}$ transitions respectively.

## Relevance

* There is some good discussion on the dependence (and reduction of) the effects of broadband light, which is measured with a grating spectrometer (enhanced by a thin etalon).
* "In [22] we demonstrated that $\lambda_{zero,lab}$ is more sensitive to $\Omega_E$ when we use circularly polarized light,magnetic fields parallel to the light propagation (along ẑ),
and atom beams with broad velocity distributions."
Where $\Omega_E$ is the earth's rotation rate, manifest through the coriolis force... See figure 6.
* Otherwise, they do not delve into detail re: polarization and B field angle, which is what I was looking for. The references in **Further Reading** might be illuminating.
* Casual comment in the conclusion: "...such as synchronized pulsing of light on atoms in a time-orbiting potential trap to control $\vec{k}\cdot\vec{B}"[3]...$

## Quality

Discuss the quality of the writing and presentation of results, and try to assess how much you trust their assumptions, approach and conclusions

## Further reading
* Refs 1-5 include other measurements of TO values
* Ref 22 for discussion of the coriolis force...
## Comments

* The introduction has a good overview of the significance and application of TO data.

## Tags
