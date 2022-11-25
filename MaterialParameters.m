% This script provides an easy way to calculate the required parameters of
% the function conductivity.m, given several experimentally accessible
% material parameters which may be found in the literature.

% Experimental material parameters (given for aluminium nitride, AlN):
w0 = 2.652e13; % LO phonon frequency, in Hz
mb = 0.285; %Band mass, in units of the electron mass
eps0 = 8.59; %Static permittivity (dimensionless)
epsInf = 4.62; %High-frequency permittivity (dimensionless)
V0 = 2.1317e-29; %Unit cell volume, in m^3
n = 1e24; %Electron density, in m^(-3) - can be tuned experimentally

% Nature constants:
hbar = 6.62607015e-34/(2*pi); %Planck constant, in J.s
me = 9.1093837015e-31; %Electron mass, in kg
e = 1.602176634e-19; %Conversion factor from J to eV, or alternatively, the electron charge in C
epsvac = 8.8541878128e-12; %Vacuum permittivity, in F/m
kB = 1.380649e-23; %Boltzmann constant, in J/K;

% Derived parameters in the model:
hw0 = hbar.*w0./e; %LO phonon frequency, in eV
EF_J = hbar.^2./(2*me*mb).*(3.*pi.^2.*n).^(2/3); %Fermi energy, in J;
EF = EF_J./e; %Fermi energy, in eV
hw_pl = hbar./e*sqrt(n*e.^2./(me*mb.*epsvac.*epsInf)); %Plasma frequency, in eV
wFw0 = EF./hw0; %Ratio of the Fermi and phonon energies
TF = EF_J/kB; %Fermi temperature, in K

% Input parameters of conductivity.m:
alpha = 1./(2*hbar*w0)*e^2/(4*pi*epsvac)*sqrt(2*mb*me*w0/hbar)*(1/epsInf-1/eps0); %Fröhlich coupling constant
V0dim = V0/(hbar/(2*mb*me*w0)).^(1.5); % Dimensionless unit cell volume
Eryw0 = (alpha/(1-epsInf/eps0))^2; %Ratio of the rydberg and phonon energies
rs = 3.*pi./16.*(9.*pi./4)^(1./3).*(hw_pl./EF).^2; %Wigner-Seitz radius