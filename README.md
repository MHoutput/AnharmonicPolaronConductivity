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

The script MaterialParameters.m provides a way to calculate most of the above parameters from more familiar parameters.

The script main.m provides an example use, which calculates and plots the optical conductivity for the material parameters given in figure 5c.
