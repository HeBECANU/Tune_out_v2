
## Article title:
High-precision measurements of the $^{87}$Rb $D$-line tune-out wavelength

### Authors:
R. H. Leonard, A. J. Fallon, C. A. Sackett, M. S. Safranova

#### Reference:
[Physical Review A **92**, 051501 (2015)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.92.052501)


## Purpose

The authors measure a tuneout wavelength in the $F=2$ ground state of Rubidium. I am reading this paper to understand how they treat the relation between the lightfield vector and the magnetic field vector.

## Context

The precision of this measurement brings it into the realm where QED, the Breit interaction, hyperfine interactions, and core-valence electron correlations become relevant, as a test for atomic structure theory.

## Approach

Using a Bragg beam, the authors separate a trapped BEC into two clouds. While the clouds remain separated, they induce a phase shift by illuminating one with a Stokes beam near the tuneout. The accumulated phase is increased by increasing the power in the Stokes beam, and the slope of the phase accumulation increases with distance from the tuneout. By measuring this slope for varying wavelengths, the authors use a linear fit to determine where the slope goes to zero. This determines the TO, but this measurement is sensitive to the light polarization. By varying the polarization of the light and fitting the dependence of the TO with a calculated value of the tensor polarizability, the authors obtain a more accurate measurement of the TO through equation (3).


## Contribution

Describe their results and conclusions, assess how they fit into the field

## Relevance

Reproducing the expression provided by Li-Yan in page 2;
"The energy shift can be expressed as
$$
U = -\frac{\langle E^2\rangle}{2} \Big(\alpha^{(0)}-\mathcal{V}\cos\chi \frac{m_F}{F}\alpha^{(1)} + \Big(\frac{3\cos^2\xi -1}{2}\Big)\frac{3m_f^2-F(F-1)}{F(2F-1)}\alpha^{(2)}\Big)
$$

Where $\alpha^{(i)}$ are components of polarizabilsity (scalar, vector, tensor), the atom is in the state $|F,m_F\rangle$ relative to the field direction $\hat{b}=\textbf{B}/B$, the angle between the Stark beam k-vector and the magnetic field is $\chi$ (so $\cos\chi = \hat{k}\cdot\hat{b}$). Similarly $\cos\xi = \hat{\epsilon}\cdot\hat{b} $ is the projection of the light polarization vector onto the magnetic field, and $\mathcal{V}$ is the fourth Stokes parameter, which characterizes the degree of circular polarization, expressible as $\mathcal{V}\cos\chi = i(\hat{\epsilon}^*\times\hat{\epsilon})\cdot\hat{b}$ (nb $\epsilon$ is not a unit vector)  ... The tensor contribution is small but measurable, however the vector contribution can be quite large. For $\sigma_+$ light ($\mathcal{V} = -1$ and $\chi = 0$) the vector term completely eliminates the TO between D1 and D2 transitions, since the light does not couple our ground state to any states in the D1 manifold... It is necessary to keep $|\mathcal{V}\cos\chi<10^{-5}|$... Much below th elevel that can typically be maintained when a laser beam passes through a vacuum chamber window.
"
The authors address this by
* The rotating TOP field causes a sign change in $\chi$ which averages to zero
* By modulating the Stokes beam in phase with the TOP field, $\cos\chi$ can be forced close to $\pm 1$ and the polarization optics (zero-order HWP and QWP) adjusted to ensure the measured phase shift is the same in each case.
* Making this check after the data run allowed for an estimate of the drift

The authors note that polarization drift is their greatest source of error, about twice the average statistical error!

*The authors constrain the variation of $B$ by measureing the Zeeman linewidth by rf spectroscopy.

## Quality

This is a very readable paper. I think that being somewhat light on math helps. The experimental design is clean, with apparently few caveats, and the results are presented in a way that is easy to undertand. Much better than the Potassium paper, hah.
