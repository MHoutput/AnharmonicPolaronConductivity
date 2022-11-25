function S = Sint(omega,model,rs,T)
% This function calculates and returns the structure factor of the
% electron gas, integrated over all wavevectors k. The integrated
% structure factor is a dimnensionless quantity, defined in figure 3
% of the accompanying article (arXiv:2210.10696).

% Inputs:
%   -omega: The frequencies where Sint(omega) must be calculated, in
%           units of the Fermi frequency E_F/hbar.
%   -model: Model for the structure factor. See below for a list of
%           implemented models.
%   -rs:    Wigner-Seitz radius of the electron gas, expressed in units of
%           the Bohr radius of the material.
%           (default: 1, does not need to be passed when using the
%           OnePolaron or HartreeFock models)
%   -T:     Temperature, in units of the Fermi temperature T_F = E_F/k_B.
%           (default: 0)

% Possible input models:
%    'OnePolaron': Uses the single polaron structure factor (default)
%    'HartreeFock': Takes Fermi statistics into account
%    'RPA': Also takes interactions into account, to lowest order
%    'Hubbard': Hubbard model for the structure factor
%    'UserDefined': Can be implemented by the user, currently returns zero

% % Example use: Calculate and plot Hartree-Fock structure factor
% omega = 0:0.01:4;
% S = Sint(omega,'HartreeFock');
% plot(omega,S)

possibleModels = {'OnePolaron','HartreeFock','RPA','Hubbard','UserDefined'};

if nargin < 1
    error('Function Sint was called with zero arguments.')
end
if nargin < 2 %Model and temperature were not passed
    model = 'OnePolaron'; %Default model is OnePolaron
end
modelIndex = find(strcmpi(model,possibleModels),1,'first');
if nargin < 3 %Wigner-Seitz radius was not passed
    rs = 1; %Default temperature is one
    if modelIndex > 2
        % The Wigner-Seitz radius is not necessary to use the OnePolaron
        % or HartreeFock models. Return an error otherwise:
        error(['The Wigner-Seitz radius rs is needed to use the ',...
            model,' model']);
    end
end
if nargin < 4 %Temperature was not passed
    T = 0; %Default temperature is zero
end

if isempty(modelIndex)
    error(['Unknown model ',model,'. Implemented models are: ',...
        strjoin(possibleModels,', '),'.'])
end

if any(size(T)>1)
    error('The temperature T must be a scalar value.') 
end
if any(size(rs)>1)
    error('The Wigner-Seitz radius rs must be a scalar value.') 
end

%We will explicitly exploit the symmetry Sint(-omega) = -Sint(omega):
signOmega = sign(omega);
omega = abs(omega);


wplwF = sqrt(0.884581919475267*rs); %The plasma frequency in our units

%For some models we will need Simpson integration. The following function
%creates a column vector of weights so that Simpson integration of a
%vector v is v*weights(length(v)).*h, where h is the distance between two
%grid points.
function ws = weights(N)
    if mod(N,2) == 0
        error('Simpson integration requires an odd number of points.')
    end
    ws = 4/3.*ones(N,1);
    ws(2:2:(N-1)) = 2/3;
    ws(1) = 1/3;
    ws(end) = 1/3;
end
Ngrid = 1000; %The typical amount of points a grid has.
%Double this to check convergence with respect to the grid



%Now, we go ahead to the different implemented models:
if modelIndex == 1
    %One polaron model, which is the omega -> ±Inf limit of other models
    %The integrated structure factor is known exactly:
    S = 0.5.*sqrt(omega);
elseif modelIndex == 2
    
    %Hartree-Fock model
    mu = ChemicalPotential(T);
    reSqrt = @(x) real(sqrt(x));
    S0 = 0.1875.*(mu>0).*(reSqrt(mu).*reSqrt(mu+omega).*(2.*mu+omega)...
        -omega.^2.*log(reSqrt(mu./(abs(omega)+1e-10))+...
        reSqrt(mu./(abs(omega)+1e-10)+sign(omega))+1e-10).*(mu+omega>0) ...
        - (reSqrt(mu).*reSqrt(mu-omega).*(2.*mu-omega)...
        -omega.^2.*log(reSqrt(mu./(abs(omega)+1e-10))+...
        reSqrt(mu./(abs(omega)+1e-10)-sign(omega))+1e-10).*(mu-omega>0)));
    if T==0
        S2 = 0; %At temperature zero, the above formula is exact
    else
        % At finite temperature, add another contribution
        fint = @(x,omega) reSqrt(x).*reSqrt(x+omega)...
            -reSqrt(x).*reSqrt(x-omega);
        S2 = 0.75*T.*integral(...
            @(x) (fint(mu+x.*T,omega)-fint(mu-x.*T,omega))./(exp(x)+1),...
            0,+Inf,'ArrayValued',true);
    end
    S = S0+S2;
elseif (modelIndex == 3 || modelIndex == 4)
    %RPA / Hubbard model
    if modelIndex == 3
        % RPA model, which is equivalent to the Hubbard model with no
        % Hubbard factor
        G = @(k) zeros(size(k));
        Gd = @(k) zeros(size(k));
    else
        %Hubbard model
        %The following Hubbard structure factor can be general, but we
        %need its derivative as well. We choose Hubbard's original
        %expression, equation (5.197) in G. Mahan, "Many particle physics",
        %third edition.
        G = @(k) 0.5.*k.^2./(k.^2+1+0.75*wplwF.^2);
        Gd = @(k) (1+0.75*wplwF.^2).*k./(k.^2+1+0.75*wplwF.^2).^2;
    end

    w = omega;
    g = @(k) k.^2./(1-G(k));
    kMax = @(w) 2.5.*(sqrt(1+w)+1);
    S = nan(size(omega));

    kpl1 = plasmaK(w,rs,T,G);
    %All indices where we have to be careful for a plasmon peak:
    plasmonWs = find(~isnan(kpl1));
    %All indices where we don't have to watch out for a plasmon peak:
    normalWs = find(isnan(kpl1));
    
    Ppl = nan(size(kpl1));
    dPpl = nan(size(kpl1));
    [Ppl(plasmonWs),dPpl(plasmonWs)] = pRPA(kpl1(plasmonWs),w(plasmonWs),rs,T,true);
    PHpl = (1-G(kpl1)).*Ppl;
    dPHpl = -Gd(kpl1).*Ppl+(1-G(kpl1)).*dPpl;
    
    %Corrected location of the plasmon peak:
    kpl = kpl1 + (1-real(PHpl))./(real(dPHpl)+1e-10);
    %Approximate width of the plasmon peak:
    kWidth = abs(imag(PHpl))./(abs(real(dPHpl))+1e-10);
    %We use the approximate width of the plasmon peak to determine whether
    %we can simply integrate it on the same grid as the remaining structure
    %factor, or whether we must integrate it on a separate grid or even
    %approximate it as a delta function.

    % First we form the k-integral for the frequencies that do not have a
    % plasmon peak: they can be done all at once, in vectorized form
    
    N = 3.*Ngrid+1; %Must be odd to apply Simpson's rule
    wNormals = w(normalWs);
    if size(wNormals,1) == 1
        wNormals = wNormals.'; % Make sure it is a column vector
    end
    nw = length(wNormals);
    kNormals = nan(nw,N);
    for i = 1:nw
        kNormals(i,:) = linspace(0,kMax(wNormals(i)),N);
    end
    wNormalsMesh = repmat(wNormals,1,N);
    hs = kMax(wNormals)./(N-1);
    PNormals = pRPA(kNormals,wNormalsMesh,rs,T);
    PHNormals = (1-G(kNormals)).*PNormals;
    SNormals = -2/pi./wplwF.^2.*g(kNormals).*imag(1./(1-PHNormals))...
        .*kNormals.^2;
    SNormalInts = SNormals*weights(N).*hs; %Simpson integration
    S(normalWs) = SNormalInts;

    % Next, we do the frequencies that do have a plasmon peak
    % involved. We have to do these one by one.
    for i = plasmonWs
        wNow = w(i);
        if abs(kWidth(i))/(abs(kpl(i))+1e-10) < 1e-4
            % First case: the plasmon peak is too narrow to
            % treat numerically. We approximate it as a delta
            % function.

            % Everything except the peak is just a continuous function,
            % which is integrated with the Simpson method:
            N = 3.*Ngrid+1;
            k1 = linspace(0,kMax(wNow),N);
            h = kMax(wNow)./(N-1);
            P1s = pRPA(k1,wNow,rs,T);
            PH1s = (1-G(k1)).*P1s;
            pow = 0.6;
            delta = 1e-10 + ...
                abs(imag(PHpl(i))).^pow.*max(abs(imag(PH1s))).^(1-pow);
            S1s = -2/pi./wplwF.^2.*g(k1).*imag(PH1s)./ ...
                ((1-real(PH1s)).^2+imag(PH1s).^2+delta.^2).*k1.^2;
            Sint2 = S1s*weights(N).*h;

            %The delta peak itself can be integrated analytically to give:
            Sint1 = 2./wplwF.^2.*kpl(i).^2.*g(kpl(i))./...
                (abs(real(dPHpl(i)))+1e-10).*...
                (1+abs(imag(Ppl(i)))./sqrt(imag(Ppl(i)).^2+delta.^2+1e-20));

            S(i) = Sint1 + Sint2;

        else
            %The plasmon peak is not a delta peak. We split the integral
            %into three parts: before the peak [0,kBorder1], the peak
            %itself [kBorder1,kBorder2], and after the peak
            %[kBorder2,+Inf].
            %Each subintegral is evaluated on a different grid.
            kBorder1 = kpl(i) - 20*kWidth(i);
            kBorder2 = kpl(i) + 20*kWidth(i);
            if kBorder1 > 0
                %Second case: the plasmon peak is not narrow
                %enough to approximate as a delta function, but
                %still quite narrow

                %Properties for the three grids:
                N1 = Ngrid+1;
                N2 = Ngrid+1;
                N3 = 2.*Ngrid+1;
                k1s = linspace(0,kBorder1,N1);
                k2s = linspace(kBorder1,kBorder2,N2);
                k3s = linspace(kBorder2,kMax(wNow),N3);
                h1 = kBorder1./(N1-1);
                h2 = (kBorder2-kBorder1)./(N2-1);
                h3 = (kMax(wNow)-kBorder2)./(N3-1);

                %Integral over the first grid:
                P1s = pRPA(k1s,wNow,rs,T);
                PH1s = (1-G(k1s)).*P1s;
                S1s = -2/pi./wplwF.^2.*g(k1s).*imag(1./(1-PH1s)).*k1s.^2;
                Sint1 = S1s*weights(N1).*h1;

                %Integral over the second grid:
                P2s = pRPA(k2s,wNow,rs,T);
                PH2s = (1-G(k2s)).*P2s;
                S2s = -2/pi./wplwF.^2.*g(k2s).*imag(1./(1-PH2s)).*k2s.^2;
                Sint2 = S2s*weights(N2).*h2;

                %Integral over the third grid:
                P3s = pRPA(k3s,wNow,rs,T);
                PH3s = (1-G(k3s)).*P3s;
                S3s = -2/pi./wplwF.^2.*g(k3s).*imag(1./(1-PH3s)).*k3s.^2;
                Sint3 = S3s*weights(N3).*h3;

                S(i) = Sint1 + Sint2 + Sint3;

            else
                %Third case: the plasmon peak is too broad, so
                %we should ignore the first grid altogether

                %Properties for the three grids:
                N2 = Ngrid+1;
                N3 = 2.*Ngrid+1;
                k2s = linspace(0,kBorder2,N2);
                k3s = linspace(kBorder2,kMax(wNow),N3);
                h2 = kBorder2./(N2-1);
                h3 = (kMax(wNow)-kBorder2)./(N3-1);

                %Integral over the 'second' grid:
                P2s = pRPA(k2s,wNow,rs,T);
                PH2s = (1-G(k2s)).*P2s;
                S2s = -2/pi./wplwF.^2.*g(k2s).*imag(1./(1-PH2s)).*k2s.^2;
                Sint2 = S2s*weights(N2).*h2;

                %Integral over the 'third' grid:
                P3s = pRPA(k3s,wNow,rs,T);
                PH3s = (1-G(k3s)).*P3s;
                S3s = -2/pi./wplwF.^2.*g(k3s).*imag(1./(1-PH3s)).*k3s.^2;
                Sint3 = S3s*weights(N3).*h3;

                S(i) = Sint2 + Sint3;

            end
        end
    end

elseif (modelIndex == 5)
    % User defined model: if you want to use a different model for the
    % dynamical structure factor, it can be implemented here. Choose a
    % model for the dynamical structure factor S(k,omega), and then
    % calculate:
    % Sint(omega) = E_F/(hbar*k_F^3)*integral(@(k) S(k,omega).*k.^2,0,+Inf)
    
    warning('This model is not yet implemented; zero was returned.')
    S = zeros(size(omega));
    
else
    error(['Unknown model ',model,'. Implemented models are: ',...
        strjoin(possibleModels,', '),'.'])
end

S = signOmega.*S;

end