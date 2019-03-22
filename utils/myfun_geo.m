function diffSIF = myfun_geo(k, D, Kc, radius, s, Sw, conf, PH, Kp, Ko, Kw, Kk, K_sh, v_sh)
% find k such that SIF = Kc

p = calcP_geo(D, s, k, Sw, conf, PH, Kk, Ko, Kw, Kp, K_sh, v_sh, radius);
[~, SIF] = getKI(k*radius, radius, p+conf, 50);
diffSIF = SIF - Kc;