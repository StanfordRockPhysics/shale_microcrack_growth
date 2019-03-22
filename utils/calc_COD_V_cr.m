function [w, V_cr] = calc_COD_V_cr(R, k, E_sh, v_sh, P)
% calc_COD_V_cr calculates the median of crack aperture and the crack
% volume

% Inputs:
% R: radius of the sphere in m
% k: ratio of pre-existing defect length to R (dimensionless)
% E_sh: Young's modulus of shale in MPa
% v_sh: Poisson's ratio of shale (dimensionless)
% P: net tensile stress (MPa)

% Outputs:
% w (1x1): median of crack aperture in m
% V_cr: microcrack volume in m^3

if min(size(P)) > 1
    disp('P has to be a vector');
    w = 0; V_cr = 0;
    return;
end

G = E_sh/2/(1+v_sh);

w = zeros(1, length(P)); V_cr = zeros(1,length(P));
for j = 1:length(P)
    [~, w_array, V_cr(j), ~, ~] = COD(k(j)*R, R, G, v_sh, P(j), 'axi', 'sphere', 50);
    w(j) = median(w_array);
end
