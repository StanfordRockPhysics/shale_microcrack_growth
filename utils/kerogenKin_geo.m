function TR = kerogenKin_geo(Ea, A, T0, H, t_vec)
% kerogenKin_geo calculates the transformation ratio (TR) based on hydrocarbon 
% generation kinetics for the geological settings. 

% Inputs:
% Ea: activation energy in kcal/mol 
% A: pre-exponential factor (frequency factor) (scale dependent on t_vec)
% T0: initial temperature in K or deg C
% H: heating rate in deg C per time unit
% t_vec: time steps where TR's are to be calculated

% Output:
% TR: transformation ratio at t_vec's. Ranges b/w 0 and 1

dt = t_vec(2) - t_vec(1);
R = 1.987e-3; % kcal/(K-mol)
T = (T0 + H * t_vec) + 273.15; % in K
rate = A * exp(-Ea ./(R*T));

TR = zeros(1, length(t_vec)); % TR 
for i = 2:length(t_vec)
    TR(i) = (TR(i-1) + rate(i) * dt)/(1 + rate(i) * dt);
end
