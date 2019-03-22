function [P, delP] = calcP_geo(D, s, k, Sw, sigma_0, PH, Kk, Ko, Kw, Kp, K_sh, v_sh, radius)
% calcP_geo calculates the final pore pressure based on Eqn. 15 in Yang and
% Gary (2018).

if k == 0
    V_star = 0;
else
    G_sh = 3*K_sh*(1-2*v_sh)/2/(1+v_sh);
    [~, ~, ~, ~, V_star] = COD(k*radius, radius, G_sh, v_sh, 10, 'axi', 'sphere', 50);
end

V_star = V_star / (4/3*pi*radius^3); % Vr
A = Sw/Kw + s*D*(1-Sw)/Ko + (1-s)*(1-Sw)/Kk + 1/Kp;

LHS = A + V_star;
RHS = A*PH + (D-1)*s*(1-Sw) - V_star*sigma_0;
P = RHS/LHS;
delP = P - PH;
