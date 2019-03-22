function diffSIF = myfun_lab(k, n2V, Kc, radius, s, Sw, T, conf, K_sh, v_sh, Kk, Kp)
% find k such that SIF = Kc

p = calcP_lab(n2V, s, radius, Sw, T, K_sh, v_sh, conf, k, Kk, Kp);
[~, SIF] = getKI(k*radius, radius, p+conf, 50);
diffSIF = SIF - Kc;