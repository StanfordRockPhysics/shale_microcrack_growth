function [w_star, w, V, xi, V_star] = COD(L, R, G, v, P, opt1, opt2, N)
% COD calculates the crack aperture and volume based on Nilson and Proffer
% (1984).

% Outputs:
% w: crack aperture profile
% V: crack volume integrated from revolving around the axis of symmetry

dx = L/N;
da = dx/2;
x = 0:dx:L;
x_chopped = x(2:end-1)';

outerIntegral = zeros(length(x_chopped), 1);
for i = 1:length(x_chopped)
    x_pt = x_chopped(i);
    a = x_pt+da:da:L;
    
    if length(a) <= 1
        disp('not enough a');
        break;
    end
    
    outerIntegrand = zeros(1, length(a));
    for j = 1:length(a)
        a_pt = a(j);
        [K_star, ~] = getKI(a_pt, R, P, 10);                               % using N = 10
        innerIntegral = K_star*pi/2;
        
        outerIntegrand(j) = innerIntegral* ...
            fFactor(x_pt, a_pt, R, opt1, opt2)*...
            OuterIntegrandFactor(x_pt, a_pt, R, opt1)* ...
            a_pt./sqrt(a_pt^2-x_pt^2);
    end
    
    outerIntegral(i) = trapz(a, outerIntegrand);
    
end

w_star = outerIntegral/L;
w = 4/pi*P*(1-v)/G*outerIntegral;
V = trapz(x_chopped, 2*pi*(R+x_chopped).*w);
xi = x_chopped/L;
V_star = V/P;

% Helper function
function factor = OuterIntegrandFactor(x, a, R, opt)
switch opt
    case 'planar'
        factor = 1;
    case 'axi'
        factor = (a+R)/(x+R);
end