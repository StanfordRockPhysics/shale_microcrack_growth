function TR = kerogenKin_lab(Ea, A, T_vec, heatingRate)
% kerogenKin_lab calculates the transformation ratio (TR) based on hydrocarbon 
% generation kinetics for the lab settings. 

% Inputs:
% Ea: activation energy in kcal/mol 
% A: pre-exponential factor (frequency factor) (scale dependent on t_vec)
% T_vec: temp steps where TR's are to be calculated
% heatingRate: in deg C per min

% Output:
% TR: transformation ratio at T_vec's. Ranges b/w 0 and 1

dt = (T_vec(2) - T_vec(1)) * 60 * heatingRate; % s
R = 1.987e-3; % kcal/(K-mol)
T_vec = T_vec + 273;
rate = A * exp(-Ea ./(R*T_vec));

TR = zeros(1, length(T_vec)); % TR 
for i = 2:length(T_vec)
    TR(i) = (TR(i-1) + rate(i)*dt)/(1 + rate(i)*dt);
end
