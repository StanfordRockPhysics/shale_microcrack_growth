% main_geo.m runs the simulation in the geological settings (Yang and Mavko, 2018). 
% Iris Yang
% March, 2019
%% test case
% clear all
% close all
addpath('utils');
var = 'D'; % 'D', 'R' or 'KIc'
%% common parameters (some to be overidden by different cases)
radius = 40e-6; % sphere radius in m
Kc = 0.05; % SIF in MPa-sqrt(m)
D = 1.3; % ratio of kerogen to HC density

switch var
    case 'D'
        vec = 1.2:0.2:1.6; % D values to be tested
        str = {'D = 1.2','D = 1.4','D = 1.6'}; % for legends in plotting
    case 'R'
        vec = 1e-6 * [20, 40, 60]; % R values to be tested
        str = {'R = 20 \mum','R = 40 \mum','R = 60 \mum'}; 
    case 'KIc'
        vec = 0.05:0.05:0.15; % KIc values to be tested
        str = {'K_{Ic} = 0.05 MPa-m^{1/2}', 'K_{Ic} = 0.10 MPa-m^{1/2}', 'K_{Ic} = 0.15 MPa-m^{1/2}'};
end
%% more parameters
t_geo = 15:0.4:40; % geologic time; m.y.
S = 100; % sedimentation rate; m/m.y.
H = 2.5; % heating rate; C/m.y.
T0 = 25; % surface temp in deg C
Ea = 49; % activation energy in kcal/mol
A = 5e13*(3e13); % pre-exponential factor in 1/m.y.

TR = kerogenKin(Ea, A, T0, H, t_geo);
s = 0.7*TR;
OB = (-2200 * 9.81 * S * t_geo) / 1e6; % Overburden; MPa
PH = (1e3 * 9.81 * S * t_geo) / 1e6; % Hydrostatic; MPa

Sw = 0.05;
TOC = 0.0992;

red = 0.4;
E_sh = red*10e3;
v_sh = 0.3;
K_sh = E_sh / 3 / (1-2*v_sh); % MPa
Kk = red * 8.3e3; % MPa
Kw = 2.73e3; % MPa
Ko = 0.68e3; % MPa
Ks = red * 39e3; % MPa                                                       
phi = TOC * 1.20 * 2.0/(1-Sw);
Kp = phi ./ (1/K_sh - (1+phi)./Ks); % Rock physics HB Table A4.1 Quartz with clay (Han)Kp = K_sh*phi/(b-phi);

k0   = 0.2;
k    = k0*ones(length(s),length(vec));
P    = zeros(length(s), length(vec));
delP = zeros(length(s), length(vec));
SIF  = zeros(length(s), length(vec));
L    = zeros(length(s), length(vec));
w    = zeros(length(s), length(vec)); 
V_cr = zeros(length(s), length(vec));
kt   = linspace(k0, 50 * k0, 5);
%% computing
for i = 1:length(vec)
    switch var
        case 'D'
            D = vec(i);
        case 'R'
            radius = vec(i);
        case 'KIc'
            Kc = vec(i);
    end
    
    for j = 1:length(s)
        alpha = s(j);
        conf = OB(j);
        ph = PH(j);
   
        [P(j,i), delP(j,i)] = calcP_geo(D, alpha, k0, Sw, conf, ph, Kk, Ko, Kw, Kp, K_sh, v_sh, radius);
        [~, SIF(j,i)] = getKI(k0*radius, radius, P(j,i)+conf, 50);
    
        if SIF(j,i) > Kc
            if j > 1
                ktry = k0;
            else
                ktry = k0;
            end

            kf = fzero(@(ktemp)...
                myfun_geo(ktemp, D, Kc, radius, alpha, Sw, conf, ph, Kp, Ko, Kw, Kk, K_sh, v_sh), ...
                [ktry, 25]); % looking for new k when SIF = Kc
            [P(j,i), delP(j,i)] = calcP_geo(D, alpha, kf, Sw, conf, ph, Kk, Ko, Kw, Kp, K_sh, v_sh, radius);
            k(j,i) = kf;
            [w_temp, V_temp] = calc_COD_V_cr(radius, k(j,i), E_sh, v_sh, P(j,i)+conf);
            w(j,i) = w_temp;
            V_cr(j,i) = V_temp / (4/3*pi * radius^3);
        end
    end
    L(:,i) = k(:,i) * radius;
end

%% plotting
plotMat = [P'; -OB; PH];

figure;
subplot(5, 2, 1:6);
plot(t_geo, plotMat);
xlabel('Geologic time (m.y.)');
ylabel('Pressure (MPa)');
legend([str, 'Overburden', 'Hydrostatic pore p'], 'location', 'nw');

subplot(5, 2, 7);
plot(TR, delP);
ylabel('\DeltaP (MPa)');

subplot(5, 2, 8);
plot(TR, k*radius*1e6);
ylabel('L (\mum)');

subplot(5, 2, 9);
plot(TR, 1e6*w);
xlabel('TR'); 
ylabel('w_{med} (\mum)');

subplot(5, 2, 10);
plot(TR, V_cr);
xlabel('TR');
ylabel('V_{cr} / V_{p,i}');