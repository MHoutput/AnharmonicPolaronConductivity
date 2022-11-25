# AnharmonicPolaronConductivity

This Matlab code can be used to calculate the optical conductivity (real and imaginary part), relaxation time, and spectroscopic effective mass of an anharmonic large polaron gas, in the weak interaction approximation. The code accompanies the following article:
M. Houtput and J. Tempere, Optical conductivity of an anharmonic large polaron gas. arXiv:2210.10696 (2022).
This code was used to produce figures 4-6 in the above article.

## Basic use

The code can simply be downloaded and run in Matlab, no compilation is required.
The function conductivity.m calculates the conductivity as a function of frequency, given several dimensionless material parameters which are defined in the accompanying article. The complete list is:
- The Fröhlich coupling constant, $\alpha$
- The 3-phonon coupling constant, $T_0$
- The 1-electron-2-phonon coupling constant, $T_1$
- The dimensionless unit cell volume, $V_0$
- The ratio of the Rydberg energy and phonon energy of the material, $E_{Ry} / \hbar \omega_{\text{LO}}$
- The Wigner-Seitz radius, in units of the material's Bohr radius, $r_s$
- The temperature, $T$, in units of the phonon temperature $\hbar \omega_{\text{LO}}/ k_B$

The script MaterialParameters.m provides a way to calculate most of the above parameters from more familiar parameters (phonon frequency, electron band mass, low- and high-frequency dielectric constants, unit cell volume, and electron density).

The script main.m provides an example use, which calculates and plots the optical conductivity for the material parameters given in figure 5c. In this file, the conductivity is calculated as:
`sigma = conductivity(omega,alpha,T0,T1,V0,Eryw0,rs,T,model)`

## Further use

The function conductivity.m uses the functions implemented in the folder "Internal functions", which calculate the dynamical structure factor $S(k,\omega)$, the phonon spectral function $M(k,\omega)$, the memory function $\Sigma(k,\omega)$, and so on. These functions have some additional functionality.

The function MemFunc.m calculates the memory function $\Sigma(k,\omega)$, which is directly related to the conductivity. With this function, one can calculate the polaron inverse effective mass and inverse relaxation time:
`[Sigma,tauinv,minv] = MemFunc(omega,alpha,T0,T1,V0,Eryw0,rs,T,model)`

The function Sint.m implements four models for the dynamical structure factor are implemented: the single polaron model, the Hartree-Fock model, the Lindhard (RPA) model, and the Hubbard model, all for the isotropic homogeneous electron gas. If a different model for the dynamical structure factor is desired, it can be implemented in Sint.m, under the model 'UserDefined'.
