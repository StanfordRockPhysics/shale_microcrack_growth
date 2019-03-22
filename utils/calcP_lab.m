function [P, delP] = calcP_lab(n2V, s, radius, Sw, T, K_sh, v_sh, conf, k, Kk, Kp)
% calcP_lab calculates the final pore pressure based on Eqn. 18 in Yang and
% Gary (2018).

T = T + 273;
R = 8.314; % cm3?MPa?K?1?mol?1
Vpi = 4/3*pi.*radius^3;
V0 = s*(1-Sw)*Vpi;
P0 = n2V*R*T; % MPa % abs pore pressure


A = Vpi*((1-s)*(1-Sw)/Kk + 1/Kp); % size(s)
if k ~= 0
    G_sh = 3*K_sh*(1-2*v_sh)/2/(1+v_sh);
    [~, ~, ~, ~, V_star] = COD(k*radius, radius, G_sh, v_sh, 10, 'axi', 'sphere', 50);
    coef = A + V_star;

else
    coef = A;
end

a = coef;
b = V0 + Sw*Vpi + coef*conf;
c = -V0*P0;

P = (-b+sqrt(b^2-4*a*c))/(2*a);
delP = P + conf;