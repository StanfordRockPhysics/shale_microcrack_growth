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


