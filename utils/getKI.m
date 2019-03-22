function [K_star,SIF] = getKI(L, R, P, N)
% getKI calculates the SIF based on Nilson and Proffer (1984).

% Inputs:
% L: crack length starting from the rim of the sphere
% R: radius of the sphere
% P: net tensile stress (unit same as G)
% N: number of discretization steps

% Ouputs:
% K_star: dimensionless SIF
% KI: SIF with dimension

% Discretization
dx = L / N;
dxx = dx / N;
x = [0:dxx:dx, 2*dx:dx:L-dx, (L-dx+dx/N):dxx:L];

% KI at crack tip L
integrand = fFactor(x(1:end-1), L, R, 'axi', 'sphere')./ ... % removing x = L due to singularity
    sqrt(L^2 - x(1:end-1).^2);
K_star = 2/pi * trapz(x(1:end-1), integrand);
SIF = K_star * P * sqrt(pi*L);


function [f_factor, f_rad, f_notch] = fFactor(x, L, R, conf_rad, conf_notch)
% Nilson and Proffer (1984)

% f_rad
switch conf_rad
    case 'planar'
        f_rad = 1;
    case 'axi'       
        f_rad = (x/L+R/L)/(1+R/L); % axisymmetric
end

% f_notch
switch conf_notch
    case 'sphere'
        m = 3;
    case 'cylinder'
        m = 2;
end

f_notch = 1 + 0.3*(1-x/L)*(1/(1+L/R))^(2*m);

% f_factor = f_rad*f_notch
f_factor = f_rad.*f_notch;