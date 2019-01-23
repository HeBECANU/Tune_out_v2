# Notes on polarizability

From various sources, the dynamic dipole polarizability is the sum of the scalar, vector, and tensor terms

$$
\alpha_1(\omega) = \underbrace{\alpha_1^0(\omega)}_\textrm{scalar} + \underbrace{A \cos \theta_k \frac{M_J}{J}\alpha_{1}^{1}(\omega)}_\textrm{Vector} + \underbrace{\Big(\frac{3 \cos^2\theta_p -1}{2}\Big)\frac{3M_J^2 - J(J+1)}{J(2J-1)}\alpha^{2}_{1}(\omega)}_\textrm{tensor})
$$

(c.f. Kien et al, equation (17-18)

Where $\cos\theta_p = \hat{\epsilon}\cdot\hat{b}$ and $\cos\theta_k = \hat{k}\cdot\hat{b}$, with $\textbf{B} = B\hat{b}$, $\hat{k}$ is the Poynting vector of the light field, $\hat{b}$ is the direction of the magnetic field, and $A$ is the 'degree of circular polarization', or the fourth Stokes parameter,  $A\cos\theta_k = i(\hat{\epsilon}^*\times\hat{\epsilon})\cdot\hat{b}$

Safranova et al:
* It is necessary to keep $|A\cos\theta_k|<10^{-5}$... Much below the elevel that can typically be maintained when a laser beam passes through a vacuum chamber window.

We measured the polarization purity to be something about 1 part in 1000. How does this relate?



### Stokes parameters
The Stokes parameters give a parametrization of the positive Hermitian operators acting on a two-level system. In the wikipedia convention, the first parameter ($I$) specifies the trace (i.e. Hilbert-Schmidt norm) of the operator, which for $I=1$ specifies the set of density operators. Stokes parameters include mixed polarization, but the Jones vector is just pure polarizations.
In the 'increasing phase' convention (fixing the sign of $V$), they are related to the electric field by
$$
I = |E_x|^2 + |E_y|^2
$$
$$
Q = |E_x|^2 - |E_y|^2
$$
$$
U = 2\textrm{Re}(E_xE_{y}^*)
$$
$$
V=-2\textrm{Im}(E_xE_{y}^*)
$$
