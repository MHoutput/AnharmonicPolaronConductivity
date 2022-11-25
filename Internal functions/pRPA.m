function [P,dPdk] = pRPA(k,omega,rs,T,getDerivative)
% INTERNAL FUNCTION

% This function calculates the "polarisation bubble" of the electron gas,
% in the RPA approximation (G. Mahan, "Many Particle Physics", third
% edition, p.328, equations (5.151-5.152)).
% The RPA dielectric function is eps(k,omega) = 1-P(k,omega).
% All implemented models for the dynamical structure factor can be written
% in terms of P(k,omega) and a Hubbard factor G(k).

%Inputs:
%   -k: Wavevector, in units of the Fermi wavevector kF. It can be an array
%       of the same size as omega.
%   -omega: Frequency, in units of the Fermi frequency EF/hbar. It can be
%           an array of the same size as k.
%   -rs:    Wigner-Seitz radius, in units of the Bohr radius of the
%           material
%   -T: Temperature, in units of the Fermi temperature TF = EF/kB.
%       (optional, default 0)
%   -getDerivative: If this flag is true, this function also calculates the
%                   derivative dP/dk. The result is only an estimate.
%                   (Optional, default false)

% % Example use: Plot the Hartree-Fock structure factor
% k = 0:0.1:4;
% omega = 0:0.1:4;
% rs = 12;
% [K,Omega] = meshgrid(k,omega);
% P = pRPA(K,Omega,rs,0);
% surf(K,Omega,-0.7197./rs.*K.^2.*imag(P))
% xlabel('$k/k_F$','Interpreter','latex')
% ylabel('$\omega/\omega_F$','Interpreter','latex')
% zlabel('$\omega_F S(k,\omega)$','Interpreter','latex')

if nargin < 3
    error('plasmaK needs at least three input arguments.')
end
if nargin <= 3
    T = 0;
end
if nargin <= 4
    getDerivative = false;
end
if nargin > 5
    error('plasmaK cannot take more than four input arguments.')
end
if numel(rs) > 1
    error('rs must be a scalar value.')
end
if numel(T) > 1
    error('T must be a scalar value.')
end
if numel(k)>1 && numel(omega)>1 && any(size(k))~=any(size(omega))
    error('If k and omega are arrays, their size needs to be the same.')
end
if isempty(k) || isempty(omega)
    if all(size(omega) == [0,0]) || all(size(omega)==size(k))
        P = zeros(size(k));
        dPdk = zeros(size(k));
    else
        P = zeros(size(omega));
        dPdk = zeros(size(omega));
    end
    return
end

w = omega;
wplwF = sqrt(0.884581919475267*rs);
mu = ChemicalPotential(T);
Ep = (w+k.^2)./(2.*k);
Em = (w-k.^2)./(2.*k);
% These functions appear often in the required integrands
funcA = @(x,E) (0.5.*(E-x.^2).*(log((abs(x+sqrt(E))+1e-10))- ...
    log((abs(x-sqrt(E))+1e-10)))+x.*sqrt(E)).*(E>0);
logexp = @(x) x.*(x>0) + log(1+exp(-abs(x)));
reLU = @(x) x.*(x>0);

%Our program will calculate some NaN values and then overwrite them, but we
%don't want the NaN warning to show up every time
warning('off','MATLAB:integral:NonFiniteValue')
%Calculate the real and imaginary parts, P = A + B*i:
if T <= 1e-5
    B = 3*pi/16.*wplwF.^2./k.^3.*(reLU(mu-Ep.^2)-reLU(mu-Em.^2));
    A = -3/8.*wplwF.^2./k.^3.*(funcA(Ep,mu)-funcA(Em,mu));
else
    B = 3*pi/16.*wplwF.^2./k.^3.*T.*...
        (logexp((mu-Ep.^2)./T)-logexp((mu-Em.^2)./T));
    A = -3/8.*wplwF.^2./k.^3.*integral(...
        @(x) (funcA(Ep,mu+x.*T)-funcA(Em,mu+x.*T))./(4.*cosh(x/2).^2),...
        -Inf,+Inf,'ArrayValued',true);
end
warning('on','MATLAB:integral:NonFiniteValue')

iZeros = find(k<=1e-6);
if numel(w) == 1
    A(iZeros) = wplwF.^2./w.^2;
else
    A(iZeros) = wplwF.^2./w(iZeros).^2;
end
B(iZeros) = 0;

P = A+1i*B;

if getDerivative
    % The derivative dP/dk is used to estimate the width of the undamped
    % plasmon branch in the dynamical structure factor, in order to
    % determine a grid for the integral over k
    if T <= 1e-5
        nF = @(x) (x<mu);
        dfuncA = @(x,E) (-x.*(log(abs(x+sqrt(E))+1e-10)...
            -log(abs(x-sqrt(E))+1e-10))+2.*sqrt(E)).*(E>0);
        dA = -3.*A./k + ...
            3./8.*wplwF.^2./k.^4.*(dfuncA(Ep,1).*Em-dfuncA(Em,1).*Ep);
    else
        warning('off','MATLAB:integral:NonFiniteValue')
        nF = @(x) 1./(exp((x-mu)./T)+1);
        reSqrt = @(x) real(sqrt(x));
        dA = - 0.75.*wplwF.^2.*integral(...
            @(x) 1./k.^4.*(k.*reSqrt(mu+x.*T)+(Em.*Ep./T.*tanh(x/2)-1.5).*...
            (funcA(Ep,mu+x.*T)-funcA(Em,mu+x.*T)))./(4.*cosh(x/2).^2),...
            -Inf,+Inf,'ArrayValued',true);
        
        kSmalls = find((k<0.4.*(sqrt(1+w)-1)) & (k>0.2.*(sqrt(1+w)-1)));
        for i = kSmalls
            dA(i) = - 0.75.*wplwF.^2.*integral(...
                @(x) 1./k(i).^4.*(k(i).*reSqrt(mu+x.*T)+...
                (Em(i).*Ep(i)./T.*tanh(x/2)-1.5).*(funcA(Ep(i),mu+x.*T)...
                -funcA(Em(i),mu+x.*T)))./(4.*cosh(x/2).^2),-Inf,+Inf);
        end
        dAapprox = wplwF.^2.*(4.*k.^3./(w.^2-k.^4).^2 + ...
            24/5.*k.*(w.^4+6.*w.^2.*k.^4+k.^8)./(w.^2-k.^4).^4.*...
            integral(@(x) reSqrt(mu+x.*T).^5./(4.*cosh(x/2).^2),-Inf,+Inf));
        kSmallest = find(k<0.2.*(sqrt(1+w)-1));
        dA(kSmallest) = dAapprox(kSmallest);
        warning('on','MATLAB:integral:NonFiniteValue')
    end
    dB = -3.*B./k + 3*pi/8.*wplwF.^2./k.^4.*Ep.*Em.*(nF(Ep.^2)-nF(Em.^2));
    dA(iZeros) = 0;
    dB(iZeros) = 0;
    dPdk = dA+1i*dB;
else
    dPdk = nan(size(P));
end
end

