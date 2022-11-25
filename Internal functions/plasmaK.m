function [kpl1,kpl2,w2] = plasmaK(omega,rs,T,G)
% INTERNAL FUNCTION

% This function calculates the solution of the plasmon branch equation:
%    (1-G(k))A(k,omega) = 1
% where A(k,omega) is the real part of the polarization bubble, also
% calculated in pRPA.m.
% We are interested in the solutions k_pl for a given frequency omega. This
% is the inverse function of omega_pl(k), which is what is usually
% preferred.

% Inputs:
%   -omega: Frequency, in units of the Fermi frequency EF/hbar. Can be an
%           array.
%   -rs:    Wigner-Seitz radius, in units of the Bohr radius of the
%           material
%   -T: Temperature, in units of the Fermi temperature TF = EF/kB.
%       (optional, default 0)
%   -G: Function handle describing the Hubbard factor
%       (optional, default 0)
%       Apparently, new plasmon branches can pop up when including the
%       Hubbard factor if rs is large and we use approximate models. We
%       have to be careful that this never happens.

% The plasmon branch is multivalued: for every frequency, there are up to
% two values of k that satisfy the plasmon branch equation.
% Outputs:
%   -kpl1:  First solution of the plasmon branch, which corresponds with
%           the "undamped" plasmon branch
%   -kpl2:  Second solution of the plasmon branch
%   -w2: Highest frequency for which a solution exists

% % Example use: Plot the two plasmon branches in the RPA model
% omega=0:0.01:6;
% [kpl1,kpl2] = plasmaK(omega,12);
% plot(kpl1,omega,'b-',kpl2,omega,'r-')
% xlabel('$k/k_F$','Interpreter','latex')
% ylabel('$\omega/\omega_F$','Interpreter','latex')


if nargin < 2
    error('plasmaK needs at least two input arguments.')
end
if nargin <= 2
    T = 0;
end
if nargin <= 3
    G = @(k) zeros(size(k));
end
if nargin > 4
    error('plasmaK cannot take more than four input arguments.')
end
if numel(rs) > 1
    error('rs needs to be a scalar value.')
end
if numel(T) > 1
    error('T needs to be a scalar value.')
end
if ~isa(G,'function_handle')
    error('G must be a function handle of a single variable k.')
end
wplwF = sqrt(0.884581919475267*rs);
mu = ChemicalPotential(T);
kpl1 = nan(size(omega));
kpl2 = nan(size(omega));
Ep = @(k,w) (w+k.^2)./(2.*k);
Em = @(k,w) (w-k.^2)./(2.*k);

funcA = @(x,E) (0.5.*(E-x.^2).*log(abs(x+sqrt(E))./...
    (abs(x-sqrt(E))+1e-10))+x.*sqrt(E)).*(E>0);

% Define the function that we want to find the zeros of: (1-G(k))A(k,omega)
function y = zeroFunc(k,w)
    iZeros = find(k<=1e-6);
    iOthers = find(k>1e-6);
    y = nan(size(k));
    y(iZeros) = 1-wplwF.^2./w(iZeros).^2; % Zero-frequency limit
    kOthers = k(iOthers);
    wOthers = w(iOthers);
    y(iOthers) = 1-(1-G(kOthers)).*(-3/8.*wplwF.^2./kOthers.^3.*...
        integral(@(x) (funcA(Ep(kOthers,wOthers),mu+x.*T)-...
        funcA(Em(kOthers,wOthers),mu+x.*T))./(4.*cosh(x/2).^2),...
        -30,+30,'ArrayValued',true));
end

mins = nan(size(omega));
minVals = nan(size(omega));
[omegaSorted,P] = sort(omega(:));
currentSign = -1;
% The undamped plasmon branch only exists between w1 and w2, where:
w1 = wplwF; %Plasma frequency
w2 = omegaSorted(end)+0.5;
for i = 1:numel(omega)
    prevSign = currentSign;
    w = omegaSorted(i);
    if w <= 1e-10
        mins(P(i)) = 0;
        minVals(P(i)) = 0;
    else
        [mins(P(i)),minVals(P(i))]= ...
            fminsearch(@(x) zeroFunc(x,w),sqrt(1+w)-1);
        currentSign = sign(minVals(P(i)));
        if currentSign*prevSign < 0
            if numel(omega) > 1 && i > 1
                w2 = 0.5.*(omegaSorted(i)+omegaSorted(i-1));
            else
                w2 = omegaSorted(i)-0.5;
            end
                break %We can stop, since we don't need the minima above w2
        end
    end    
end

%Loop over the omegas smaller than w2
for i = find(omega(:) < w2)'
    w = omega(i);
    zeroFunc3 = @(k) zeroFunc(k,w);
    %Find the value of k that makes zeroFunc3 zero for this specific value
    %of k:
    if w <= w1
        if w == w1
            kpl1(i) = 0;
        end
        %If w<w1, there is no solution for kpl1, so we keep it at NaN
        if w < 1e-10
            kpl2(i) = 0;
        else
            kpl2(i) = fzero(zeroFunc3, [mins(i),5]);
            %In theory 5 should be +Inf, but in practice 5 is big enough
        end
    else
        kpl1(i) = fzero(zeroFunc3, [0,mins(i)]);
        kpl2(i) = fzero(zeroFunc3, [mins(i),5]);
    end
end

end

