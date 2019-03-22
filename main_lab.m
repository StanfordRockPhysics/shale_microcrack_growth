% main_lab.m runs the simulation in the laboratory settings (Yang and Mavko, 2018). 
% Iris Yang
% March, 2019
%% test case
% clear all
% close all
addpath('utils');
var = 'n2V'; % n2V, R or KIc
%% parameters
T_vec = 60:5:400; % temp range in deg C
Sw = 0.05;
TOC = 0.0992;

red = 0.1;
E_sh = red * 10e3;
v_sh = 0.3;
K_sh = E_sh / 3 / (1-2*v_sh); % MPa
Kk = red * 8.3e3; % MPa
Kw = red * 2.73e3; % MPa
Ks = red * 39e3; % MPa
phi = TOC * 1.20 * 2.0/(1-Sw);
Kp = phi./(1/K_sh - (1 + phi)./Ks); % Rock physics HB Table A4.1 Quartz with clay (Han)Kp = K_sh*phi/(b-phi);
conf = -0.1; % 1 atm in MPa

heatingRate = 1; % deg C per min
Ea = 49; % activation energy in kcal/mol
A = 5e13; % pre-exponential factor in 1/s

k0 = 0.2;
%% common parameters (some to be overidden by different cases)
radius = 40e-6; % sphere radius in m
Kc = 0.05; % SIF in MPa-sqrt(m)
n2V = 0.04; % volume expansion parameter in mol/cm3

switch var
    case 'n2V'
        vec = [0.02, 0.04, 0.06]; % n2V values to be tested
        str = {'n/V_{k,i} = 0.02 mol/cm^3', 'n/V_{k,i} = 0.04 mol/cm^3', ...
               'n/V_{k,i} = 0.06 mol/cm^3'}; % for legends in plotting
    case 'R'
        vec = 1e-6 * [20, 40, 60]; % R values to be tested
        str = {'R = 20 \mum','R = 40 \mum','R = 60 \mum'}; 
    case 'KIc'
        vec = 0.05:0.05:0.15; % KIc values to be tested
        str = {'K_{Ic} = 0.05 MPa-m^{1/2}', 'K_{Ic} = 0.10 MPa-m^{1/2}', 'K_{Ic} = 0.15 MPa-m^{1/2}'};
end
%% computing
P    = zeros(length(T_vec), length(vec));
delP = zeros(length(T_vec), length(vec));
SIF  = zeros(length(T_vec), length(vec));
L    = zeros(length(T_vec), length(vec));
w    = zeros(length(T_vec), length(vec)); 
V_cr = zeros(length(T_vec), length(vec));
area = zeros(length(T_vec), length(vec));

for i = 1:length(vec)
    switch var
        case 'n2V'
            n2V = vec(i);
        case 'R'
            radius = vec(i);
        case 'KIc'
            Kc = vec(i);
    end
    
    TR = kerogenKin_lab(Ea, A, T_vec, heatingRate);
    s = 0.7*TR;
    
    for j = 1:length(T_vec)
        T = T_vec(j);
        
        [P(j,i), delP(j,i)] = calcP_lab(n2V, s(j), radius, Sw, T, K_sh, v_sh, conf, k0, Kk, Kp);        
        [~, SIF(j,i)] = getKI(k0*radius, radius, delP(j,i), 50);
        
        if SIF(j,i) > Kc           
            kf = fzero(@(ktemp)...
                        myfun_lab(ktemp, n2V, Kc, radius, s(j), Sw, T,  conf, K_sh, v_sh, Kk, Kp), ...
                        [k0, 1000]); % looking for new k when SIF = Kc
            [P(j,i), delP(j,i)] = calcP_lab(n2V, s(j), radius, Sw, T, K_sh, v_sh, conf, kf, Kk, Kp);            
            [w_temp, V_temp] = calc_COD_V_cr(radius, kf, E_sh, v_sh, delP(j,i));
            w(j,i) = w_temp;
            V_cr(j,i) = V_temp / (4/3*pi * radius^3);
            L(j,i) = kf * radius;
        end
    end
    
    area(:,i) = 1e6*pi * ((L(:,i)+radius).^2 - radius^2); % mm^2
end

%% plotting
figure;
x_scale = [300, 400];
subplot(5, 2, 1:6);
plot(T_vec, area);
xlabel('T (^oC)');
ylabel('Fracture surface area (mm^2)');
xlim(x_scale);
legend(str, 'location', 'nw');

subplot(5,2,7);
plot(T_vec,delP);
ylabel('\DeltaP (MPa)');
xlim(x_scale);

subplot(5, 2, 8);
plot(T_vec, L*1e3);
ylabel('L (mm)');
xlim(x_scale);

subplot(5, 2, 9);
plot(T_vec, 1e6*w);
xlabel('T (^oC)'); 
ylabel('w_{med} (\mum)');
xlim(x_scale);

subplot(5, 2, 10);
plot(T_vec, V_cr);
xlabel('T (^oC)');
ylabel('V_{cr} / V_{p,i}');
xlim(x_scale);