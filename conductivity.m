function sigma = cond(omega,alpha,T0,T1,V0,Eryw0,rs,T,model)

% This function will calculate the conductivity of an anharmonic polaron,
% to lowest order in alpha. The memory function will be expressed in units
% n*e²/m*w_0, where n is the electron density, e is the elementary charge,
% m is the electron band mass, and w_0 is the phonon frequency.

% Inputs: (can all be scalars or vectors)
%  -omega:   The frequency at which we calculate the conductivity, in units 
%            of the phonon frequency w_0.
%  -alpha:   The Fröhlich coupling constant.
%  -T0:      The 3-phonon coupling constant.
%  -T1:      The 1-electron-2-phonon coupling constant.
%  -V0:      The dimensionless volume of the unit cell.
%  -Eryw0:   The ratio of the Rydberg energy and the LO phonon energy.
%  -rs:      The Wigner-Seitz radius, expressed in units of the Bohr radius.
%  -T:       Temperature, expressed in units of the "phonon temperature"
%            hbar*w_0/k_B.
%  -model:   A string with the model used for the structure factor. It can
%            be either a single string, or a 1D cell array containing
%            multiple models.
%            Possible input models:
%            'OnePolaron':   Uses the structure factor of a single polaron.
%                            Only makes sense at zero temperature.
%            'HartreeFock': Takes Fermi statistics into account (default).
%            'RPA': Also takes interactions into account, to lowest order.
%            'Hubbard': Hubbard model for the structure factor.      

% If only omega is a vector and the other inputs are scalars, the output
% will be a single vector representing the conductivity depending on omega.
% If the inputs are vectors, the output will be (up to) 8D cell arrays,
% where the first dimension is the dependence on alpha, the second is the
% dependence on T0, and so on. Each entry in the cell arrays is a vector,
% which is the memory function depending on omega.

% It is also allowed to pass all arguments together as a cell array.

% This function simply calls MemFunc.m to calculate the memory function;
% the conductivity is related to the memory function by a simple algebraic
% formula.

% % Example use: Plot the optical absorption spectrum of a single polaron
% omega = 0:0.01:4;
% sigma = conductivity(omega,1,0,0.1,0.001,8,12,0,'OnePolaron');
% plot(omega,real(sigma),'b-')

if ~iscell(model)
    model = {model};
end
memoryFunction = MemFunc(omega,alpha,T0,T1,V0,Eryw0,rs,T,model);
if numel(memoryFunction) == 1
    sigma = 1i./(omega-memoryFunction{1});
else
    sigma = cell(length(alpha),length(T0),length(T1),length(V0),...
        length(Eryw0),length(rs),length(T),length(model));
    for i = 1:numel(memoryFunction)
        sigma{i} = 1i./(omega-memoryFunction{i});
    end
end