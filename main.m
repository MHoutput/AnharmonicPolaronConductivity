% Example script to plot the real and imaginary parts of the conductivity
% of the anharmonic polaron gas, for one specific set of values.
% This script reproduces the solid orange curve of figure 5c of the
% accompanying article (arXiv:2210.10696).

% Make sure to add the folder "Internal functions" to the MatLab path
% before running this script.

% Set the input values:
omega = 0:0.01:4;
alpha = 1;
T0 = 0;
T1 = 0.1;
V0 = 0.001;
Eryw0 = 8;
rs = 12;
T = 0;
model = 'Hubbard';

% Calculate the corresponding conductivity:
sigma = conductivity(omega,alpha,T0,T1,V0,Eryw0,rs,T,model);

% Plot the results
figure
plot(omega,real(sigma),'k-','LineWidth',2)
xlabel('$\omega/\omega_{LO}$','Interpreter','latex','FontSize',16)
ylabel('$\frac{\sigma_R(\omega)}{\frac{ne^2}{m\omega_0}}$',...
    'Interpreter','latex','FontSize',20,'Rotation',0,...
    'VerticalAlignment','middle', 'HorizontalAlignment','right')
title('Conductivity, real part','FontSize',18)
ylim([0,0.1])

figure
plot(omega,imag(sigma),'k-','LineWidth',2)
xlabel('$\omega/\omega_{LO}$','Interpreter','latex','FontSize',16)
ylabel('$\frac{\sigma_I(\omega)}{\frac{ne^2}{m\omega_0}}$',...
    'Interpreter','latex','FontSize',20,'Rotation',0,...
    'VerticalAlignment','middle', 'HorizontalAlignment','right')
title('Conductivity, imaginary part','FontSize',18)
ylim([0,2])
