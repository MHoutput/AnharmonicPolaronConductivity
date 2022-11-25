function mu = ChemicalPotential(T,approx)
% INTERNAL FUNCTION

% This function calculates the reduced chemical potential mu/E_F, given the
% reduced temperature T/T_F. T can be a scalar or an array.

% Two different calculation methods are implemented: numerical solution of
% the defining equation which gives the exact result, or a fast approximate
% formula with an error less than 1%.
% If approx is False, the exact result is returned (default).
% If approx is True, the approximate formula is used.

% %Example use:
% T = 0:0.01:1;
% mu = ChemicalPotential(T)

if nargin == 1
    %By default, we do not use the approximate formula
    approx = false;
end

%This interpolation has an error of <1%
muApprox = @(T) (1-pi^2/12.*T.^2-pi^4/80.*T.^4)./(1+exp((T-0.5)/0.05))+ ...
    (T.*log(4/(3.*sqrt(pi)))-1.5.*T.*log(T+1e-10)+ ...
    sqrt(2./(9*pi.*(T+0.01))))./(1+exp(-(T-0.5)/0.05));

if approx
    mu = muApprox(T);
else
    %Solve the implicit equation to find the exact answer 
    mu = zeros(size(T));
    reSqrt = @(x) real(sqrt(x));
    intFunc = @(y,T) y.*reSqrt(y) + 3/2.*T*integral(@(x) ...
        (reSqrt(y+x.*T)-reSqrt(y-x.*T))./(exp(x)+1),0,+Inf)-1;
    for i = 1:numel(mu)
        mu(i) = fzero(@(x) intFunc(x,T(i)),muApprox(T(i)));
    end
end

end
