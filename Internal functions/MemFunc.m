function [S,tauinv,minv] = MemFunc(omega,alpha,T0,T1,V0,Eryw0,rs,T,model)

%This function will calculate the memory function of an anharmonic polaron,
%to lowest order in alpha. The memory function will be expressed in units
%of the phonon frequency w_0. All inputs can be either scalars or vectors.
%Inputs:
% -omega:   The frequency at which we calculate the conductivity, in units 
%           of the phonon frequency w_0.
% -alpha:   The Fröhlich coupling constant.
% -T0:      The 3-phonon coupling constant.
% -T1:      The 1-electron-2-phonon coupling constant.
% -V0:      The dimensionless volume of the unit cell.
% -Eryw0:   The ratio of the Rydberg energy and the LO phonon energy.
% -rs:      The Wigner-Seitz radius, expressed in units of the Bohr radius.
% -T:       Temperature, expressed in units of the "phonon temperature"
%           hbar*w_0/k_B.
% -model:   A string with the model used for the structure factor. It can
%           be either a single string, or a 1D cell array containing
%           multiple models.
%           Possible input models:
%           'OnePolaron':   Uses the structure factor of a single polaron.
%                           Only makes sense at zero temperature.
%           'HartreeFock': Takes Fermi statistics into account (default).
%           'RPA': Also takes interactions into account, to lowest order.
%           'Hubbard': Hubbard model for the structure factor.      

%The outputs will be (up to) 8D cell arrays, where the first dimension is
%the dependence on alpha, the second is the dependence on T0, and so on.
%Each entry in the cell arrays is a vector, which is the memory function
%depending on omega.

%It is also allowed to pass all arguments together as a cell array.

%The implementation is based on the theory of the accompanying article
%(arXiv:2210.10696), especially equations (26) and (53). Besides these
%equations, it uses the integrated structure factor, which is implemented
%in the function Sint.m. This version of the function uses interpolation
%for the integrated structure factor (Sint_func.m), and is therefore
%well-suited for calculations with many values of omega at once.

% % Example use: Plot the memory function of a polaron gas
% omega = 0:0.01:4;
% S = MemFunc(omega,1,0,0.1,0.001,8,12,0,'HartreeFock');
% plot(omega,real(S{1}),'b-',omega,imag(S{1}),'r-')
% legend('Real part','Imaginary part')

if nargin == 1
    % If there is only one input argument, it is in the form of a cell
    % array containing all 9 input variables
    if iscell(omega) && length(omega) == 9
        omegas = omega{1};
        alphas = omega{2};
        T0s = omega{3};
        T1s = omega{4};
        V0s = omega{5};
        Eryw0s = omega{6};
        rss = omega{7};
        Ts = omega{8};
        models = omega{9};
    else
        error('If this function is called with one argument, it must be a cell array containing all 9 input variables.')
    end
else
    % Add "s" to the names of the array-valued variables; later in the
    % code, "omega", "alpha", ... will refer to a specific value in these
    % arrays
    omegas = omega;
    alphas = alpha;
    T0s = T0;
    T1s = T1;
    V0s = V0;
    Eryw0s = Eryw0;
    rss = rs;
    Ts = T;
    models = model;
end
%If models is just a string, we have to put it in a cell array first
if ~iscell(models)
    models = {models};
end

%Check if the inputs are what we expect. The last three get checked
%automatically in Sint.
if sum(size(omegas)>1)>1
    error('Input omega must be a scalar or a vector.')
end
if sum(size(alphas)>1)>1
    error('Input alpha must be a scalar or a vector.')
end
if sum(size(T0s)>1)>1
    error('Input T0 must be a scalar or a vector.')
end
if sum(size(T1s)>1)>1
    error('Input T1 must be a scalar or a vector.')
end
if sum(size(V0s)>1)>1
    error('Input V0 must be a scalar or a vector.')
end
if sum(size(Eryw0s)>1)>1
    error('Input wFw0 must be a scalar or a vector.')
end
if sum(size(rss)>1)>1
    error('Input rs must be a scalar or a vector.')
end
if sum(size(Ts)>1)>1
    error('Input T must be a scalar or a vector.')
end
if sum(size(models)>1)>1
    error('Input model must be a string or a 1D cell array of strings.')
end



n2 = length(alphas);
n3 = length(T0s);
n4 = length(T1s);
n5 = length(V0s);
n6 = length(Eryw0s);
n7 = length(rss);
n8 = length(Ts);
n9 = length(models);

S = cell(n2,n3,n4,n5,n6,n7,n8,n9);
minv = cell(n2,n3,n4,n5,n6,n7,n8,n9);
tauinv = cell(n2,n3,n4,n5,n6,n7,n8,n9);
w = omegas;

counter = 0;
goal = n6*n7*n8*n9;
make_waitbar = goal>1;
if make_waitbar
    han = waitbar(counter./goal,['Finished ',num2str(counter),...
        ' of ',num2str(goal),' memory functions']);
end

%Sometimes, we will need Simpson integration: the following function
%creates a column vector of weights so that Simpson integration of a
%vector v is v*weights(length(v)).*h
function ws = weights(N)
    if mod(N,2) == 0
        error('We need an odd number of points to perform Simpson integration.')
    end
    ws = 4/3.*ones(N,1);
    ws(2:2:(N-1)) = 2/3;
    ws(1) = 1/3;
    ws(end) = 1/3;
end
Ngrid = 1000; %The typical amount of points a grid has.
%Double this to check convergence with respect to the grid

zeroIndex = find(abs(w)<=1e-9);

% First, we loop over all variables that the integrated structure factor
% depends on: these are Eryw0, rs, T, and model
for i6 = 1:n6
for i7 = 1:n7
for i8 = 1:n8
for i9 = 1:n9
    Eryw0 = Eryw0s(i6);
    rs = rss(i7);
    T = Ts(i8);
    model = models{i9};

    wFw0 = Eryw0.*(9*pi/4).^(2/3)./rs.^2;
    TS = T./wFw0;

    %Get the function handle for the integrated structure factor:
    %Most of the runtime is in this fuction
    [Sf,wMaxwF] = Sint_func(model,rs,TS);

    if T > 0
        % Bose-Einstein distribution and its first three derivatives:
        nB = @(x) 1./(exp(x./T)-1);
        nBd = @(x) -1/(4.*T.*sinh(x./(2.*T)).^2);
        nBdd = @(x) coth(x./(2.*T))/(4.*T.^2.*sinh(x./(2.*T)).^2);
        nBddd = @(x) -(2+cosh(x./T))/(8.*T.^3.*sinh(x./(2.*T)).^4);
        % A coefficient that occurs often in our expressions
        cFunc = @(x) coth(x./(2.*T));
    else
        %Temperature zero limits of the above expressions
        nB = @(x) (x>0)-1;
        nBd = @(x) zeros(size(x));
        nBdd = @(x) zeros(size(x));
        nBddd = @(x) zeros(size(x));
        cFunc = @(x) sign(x);
    end
    cFI = cFunc(1); %coth(1/2T) , normally
    if strcmpi(model,'onepolaron')
        %There is one coefficient that we must explicitly set equal to 1
        cFR = 1;
    else
        cFR = cFunc(1);
    end

    % Loop over T0 and V0, and calculate the phonon spectral function
    % M(k,omega):
    for i3 = 1:n3
    for i5 = 1:n5
        T0 = T0s(i3);
        V0 = V0s(i5);

        kFap = sqrt(wFw0);

        %Poles of the phonon Green's function:
        sqrtT0 = sqrt(1+64./135.*T0.^2./V0.*cFR);
        x1 = sqrt(0.5.*(5-3.*sqrtT0)); %Approximately equal to 1
        x2 = sqrt(0.5.*(5+3.*sqrtT0)); %Approximately equal to 2

        if strcmpi(model,'onepolaron')
            % Our implementation doesn't work for the one polaron model
            % because of divergences caused by Sint'(0) = +Inf.
            % However, in this case, we can calculate the memory function
            % analytically.
            f = @(z) (2-sqrt(z+1)+1i.*sqrt(z-1))./z;
            S1 = 0.5./sqrt(wFw0.*x1).*f(w./x1);
            S2 = 0.5./sqrt(wFw0.*x2).*f(w./x2);
            I1_0 = -nBd(x1).*sqrt(x1./wFw0);
            I2_0 = -nBd(x2).*sqrt(x2./wFw0);
            intFunc = @(w,a) 0.25.*(sign(w-a).*cFunc(a)+cFunc(abs(w-a))) ...
                .*sqrt(abs(w-a))./w;
            intFuncdd0 = @(a) 0.5.*(nBd(a)./(4.*a.^1.5) ...
                - nBdd(a)./(2.*a.^0.5) - 1/3*nBddd(a).*a.^0.5);
            omegaMin = 1e-3;
            Rd1_0 = 2./pi./sqrt(wFw0).*...
                (omegaMin.*intFuncdd0(x1) + integral(...
                @(w) (intFunc(w,x1)-intFunc(w,-x1)+nBd(x1).*sqrt(x1))./w.^2,...
                omegaMin,+Inf));
            Rd2_0 = 2./pi./sqrt(wFw0).*(omegaMin.*intFuncdd0(x2)...
                + integral(@(w) ...
                (intFunc(w,x2)-intFunc(w,-x2)+nBd(x2).*sqrt(x2))./w.^2,...
                omegaMin,+Inf));

            nB = @(x) (x>0)-1;
            nBd = @(x) zeros(size(x));
            cFunc = @(x) sign(x);
        else
            %Calculate the memory function from equation (53) of the
            %accompanying article.
            %The core function that appears in (53) is:
            %intFunc = @(w,a) (1+nB(a)+nB(w-a)).*Sf((w-a)./wFw0);
            %but this doesn't behave well around w = a.
            %Instead, we use the following numerically stable form:
            intFunc = @(w,a) 0.5.*(sign(w-a).*cFunc(a)+...
                cFunc((abs(w-a)+1e-9))).*Sf((abs(w-a)+1e-9)./wFw0);
            %From this function, two parts appear in the memory
            %function:
            I1func = @(w) 1./(abs(w)+1e-9).*...
                (intFunc(abs(w)+1e-9,x1)-intFunc(abs(w)+1e-9,-x1));
            I2func = @(w) 1./(abs(w)+1e-9).*...
                (intFunc(abs(w)+1e-9,x2)-intFunc(abs(w)+1e-9,-x2));
            
            %Calculate the contributions to the imaginary part of the
            %approximate memory function:
            I1 = I1func(w);
            I2 = I2func(w);

            %We also calculate the values at omega=0 explicitly (useful for
            %the relaxation time):
            if isempty(zeroIndex)
                I1_0 = -2.*nBd(x1).*Sf(x1./wFw0);
                I2_0 = -2.*nBd(x2).*Sf(x2./wFw0);
            else
                I1_0 = I1(zeroIndex);
                I2_0 = I2(zeroIndex);
            end

            %Now, we calculate the real parts of the imaginary
            %contributions via the Kramers-Kronig relations
            
            %The integral is split into two parts: main features [0,wMax1],
            %and the high frequency region [wMax1,wMax2] up to the upper
            %interpolation cutoff.
            wMax1 = 4;
            wMax2 = max([wMaxwF*wFw0,wMax1+1,max(w)+1]);
            %For the region [wMax2,+Inf], it is assumed that the structure
            %factor is approximately the one polaron structure factor, and
            %the integral is performed analytically
            
            N1 = Ngrid+1;
            N2 = Ngrid+1;
            h1 = wMax1./(N1-1);
            h2 = (wMax2-wMax1)./(N2-1);
            ht = 1e-9; %A tiny step
            hd = 1e-3; %A step used to approximate the derivative
            t1 = linspace(0,wMax1,N1); %Integration variable on grid 1
            t2 = linspace(wMax1,wMax2,N2); %Integration variable on grid 2
            [T1,W1] = meshgrid(t1,w);
            [T2,W2] = meshgrid(t2,w);
            integrands = @(t,w,f) (f(t)-f(w+ht))./(t.^2-(w+ht).^2);
            integrand1 = integrands(T1,W1,I1func);
            integrand2 = integrands(T1,W1,I2func);
            lowIndices = (T1 < ht)&(W1 < hd);
            integrand1(lowIndices) = (I1func(hd)-I1func(0))./hd.^2;
            integrand2(lowIndices) = (I2func(hd)-I2func(0))./hd.^2;
            R1 = 2./pi.*w.*( (integrand1*weights(N1).*h1).' ... %First grid
                + (integrands(T2,W2,I1func)*weights(N2).*h2).' ) ... %Second grid
                + (1+2*nB(x1))./(pi*sqrt(w.*wFw0)).*...
                (atanh(sqrt(w/wMax2))-atan(sqrt(w/wMax2)))...
                - 2./pi.*atanh(w/wMax2).*I1; %Analytic contribution
            R2 = 2./pi.*w.*( (integrand2*weights(N1).*h1).' ... %First grid
                + (integrands(T2,W2,I2func)*weights(N2).*h2).' ) ... %Second grid
                + (1+2*nB(x2))./(pi*sqrt(w.*wFw0)).*...
                (atanh(sqrt(w/wMax2))-atan(sqrt(w/wMax2)))...
                - 2./pi.*atanh(w/wMax2).*I2; % Analytic contribution
            R1(zeroIndex) = 0;
            R2(zeroIndex) = 0; %Real part is zero when omega=0

            %Calculate the derivative of the real part at zero frequency
            %(useful for the effective mass):
            integrand01 = integrands(t1,0,I1func);
            integrand01(1) = (I1func(hd)-I1func(0))./hd.^2;
            integrand02 = integrands(t1,0,I2func);
            integrand02(1) = (I2func(hd)-I2func(0))./hd.^2;
            Rd1_0 = 2./pi.*( (integrand01*weights(N1).*h1).' ... % First grid
                + (integrands(t2,0,I1func)*weights(N2).*h2).' ) ... %Second grid
                + (1+2*nB(x1))./(1.5*pi*sqrt(wFw0)).*wMax2.^(-1.5) ...
                - 2./pi.*I1_0./wMax2; %Analytic contribution
            Rd2_0 = 2./pi.*( (integrand02*weights(N1).*h1).' ... % First grid
                + (integrands(t2,0,I2func)*weights(N2).*h2).' ) ... %Second grid
                + (1+2*nB(x2))./(1.5*pi*sqrt(wFw0)).*wMax2.^(-1.5) ...
                - 2./pi.*I2_0./wMax2; %Analytic contribution

            %Add real and imaginary parts:
            S1 = R1 + 1i.*I1;
            S2 = R2 + 1i.*I2;
        end




        %Given S1 and S2, there is an algebraic expression for the memory
        %function, in terms of the remaining variables alpha and T1:
        for i2 = 1:n2
        for i4 = 1:n4
            alpha = alphas(i2);
            T1 = T1s(i4);

            c1R = 0.5./x1.*(1+4/15*T1.^2./V0.*cFR...
                + (1+(8/3*T0.*T1-T1.^2).*4/15./V0.*cFR)./sqrtT0);
            c2R = 0.5./x2.*(1+4/15*T1.^2./V0.*cFR...
                - (1+(8/3*T0.*T1-T1.^2).*4/15./V0.*cFR)./sqrtT0);
            c1I = 0.5./x1.*(1+4/15*T1.^2./V0.*cFI...
                + (1+(8/3*T0.*T1-T1.^2).*4/15./V0.*cFI)./sqrtT0);
            c2I = 0.5./x2.*(1+4/15*T1.^2./V0.*cFI...
                - (1+(8/3*T0.*T1-T1.^2).*4/15./V0.*cFI)./sqrtT0);

            %Calculate the approximate memory function with equation (53):
            SUncorr = -4.*alpha./3.*kFap.*(c1R.*S1 + c2R.*S2);
            %Calculate the approximate memory function and its derivative
            %at omega=0:
            SUncorr_0 = -4i.*alpha./3.*kFap.*(c1I.*I1_0 + c2I.*I2_0);
            SUncorrd_0 = -4.*alpha./3.*kFap.*(c1I.*Rd1_0 + c2I.*Rd2_0);

            %Calculate the inverse relaxation time and effective mass from
            %equations (68)-(69):
            tauinv{i2,i3,i4,i5,i6,i7,i8,i9} = -imag(SUncorr_0);
            minv{i2,i3,i4,i5,i6,i7,i8,i9} = 1 + SUncorrd_0;

            if strcmpi(model,'onepolaron')
                SUncorr_0 = 0; %For the calculation of the memory function, we assume T=0
            end
            
            %Calculate the full memory function with equation (14):
            Scorr = SUncorr./(1+(SUncorr-SUncorr_0)./w);
            Scorr(zeroIndex) = SUncorr_0./(1+SUncorrd_0);
            S{i2,i3,i4,i5,i6,i7,i8,i9} = Scorr;

        end
        end
    end
    end

    counter = counter+1;
    if make_waitbar
        waitbar(counter./goal,han,['Finished ',num2str(counter),' of ',num2str(goal),' memory functions']);
    end
end
end
end
end
if make_waitbar
    close(han);
end

end