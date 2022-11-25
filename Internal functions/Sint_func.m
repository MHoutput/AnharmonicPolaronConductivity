function [S_han,omegaMax] = Sint_func(model,rs,T)
% INTERNAL FUNCTION

% This function generates a function handle for the integrated structure
% factor, based on interpolation.

% The actual implementation of the integrated structure factor is in
% Sint.m; this function only interpolates the result.
% The possible input models implemented in Sint.m are:
%   'OnePolaron': Uses the structure factor of a single polaron (default)
%   'HartreeFock': Takes Fermi statistics into account
%   'RPA': Also takes interactions into account, to lowest order
%   'Hubbard': Hubbard model for the structure factor

% Units:
%   T is written in units of the Fermi temperature T_F = E_F/k_B.
%   rs is written in units of the Bohr radius of the material

% The output S_han is a function handle that represents
% S_int(omega/omega_F), where omega_F = E_F/hbar is the frequency
% corresponding to the Fermi energy
% The function also outputs omegaMax, which is the largest value that
% was used for interpolation. It can typically be ignored.

% %Example use:
% Sint_han = Sint_func('HartreeFock',12,0);
% omega = 0:0.01:4;
% plot(omega,Sint_han(omega))

possibleModels = {'OnePolaron','HartreeFock','RPA','Hubbard'};

if nargin < 1 %Model and temperature were not passed
    model = 'OnePolaron'; %Default model is OnePolaron
end
modelIndex = find(strcmpi(model,possibleModels),1,'first');
if nargin < 2
    rs = 1; %Default Wigner-Seitz radius is one
    if modelIndex > 2
        % The Wigner-Seitz radius is not necessary to use the OnePolaron
        % or HartreeFock models. Return an error otherwise:
        error(['The Wigner-Seitz radius rs is needed to use the ',model,...
            ' model']);
    end
end
if nargin < 3 %Temperature was not passed
    T = 0; %Default temperature is zero
end

if isempty(modelIndex)
    error(['Unknown model ',model,...
        '. Implemented models are: ',strjoin(possibleModels,', '),'.'])
end

if any(size(T)>1)
    error('The temperature T must be a scalar value.') 
end
if any(size(rs)>1)
    error('The Wigner-Seitz radius rs must be a scalar value.') 
end 

% There are two important frequencies:
wF = 1; %The Fermi frequency
wpl = sqrt(0.884581919475267*rs); %The plasma frequency
wScale = 2*max(wF,wpl); %We use a fine grid up to at least this frequency 

if modelIndex == 1
    %For the one-polaron model, the result is known exactly:
    S_han = @(x) sign(x).*0.5.*sqrt(abs(x));
    omegaMax = 100;
else
    %All other models can be handled in the same way
    N1 = 201;
    alpha = log(1.2);
    wMax = 100;
    h = wScale./(N1-1);
    N2 = ceil((log(1+(exp(alpha)-1).*wMax./h))./alpha);
    omegas = [linspace(0,wScale,N1),wScale+...
        h.*(exp(alpha.*(1:N2))-1)./(exp(alpha)-1)];
    omegaMax = omegas(end);
    S_vals = Sint(omegas,model,rs,T);
    omegaStep = omegas(2); % The step size of the first grid
    Sd0 = S_vals(2)./omegaStep; %Approximate derivative at 0

    mu = ChemicalPotential(T);
    reSqrt = @(x) real(sqrt(x));
    I1 = integral(@(x) reSqrt(mu+x.*T).^3./(4.*cosh(0.5.*x).^2), -30, 30);
    S_asy = @(x) sqrt(x)./2.*(1+0.3.*I1./(x+1e-5));
    S_vals_red = S_vals./S_asy(omegas);
    S_vals_red(1) = 0;
    S_pp = spline(omegas,S_vals_red);
    S_han = @(x) (abs(x)<omegaStep).*Sd0.*x + ...
        (abs(x)>=omegaStep).*sign(x).*S_asy(abs(x)).*...
        (1-(1-ppval(S_pp,abs(x))).*(abs(x)<omegaMax));
end
    
end