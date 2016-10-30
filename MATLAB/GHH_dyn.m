function out = GHH_dyn(in,param)

kp    = in(1);
hp    = in(2);
taup  = in(3);
k     = in(4);
h     = in(5);
tau   = in(6);
km    = in(7);
hm    = in(8);
taum  = in(9);
zp    = in(10);
z     = in(11);

alf = param(1);
bet = param(2);
sig = param(3);
theta = param(4);
del = param(5);
g = param(6);
psi   = param(7); 
omega = param(8);
rho = param(9);

[Y, r, w, c, i, Gam, Del, Pi] = ...
    GHH_defs(km, hm, taum, z, k, h, tau, param);
[Yp, rp, wp, cp, ip, Gamp, Delp, Pip] = ...
    GHH_defs(k, h, tau, zp, kp, hp, taup, param);

out1 = bet*(1+rp-del)*Gamp*Delp/(Gam*Del) - 1;
out2 = Del*w/Pi - 1;
out3 = tau - taum - 1;
out = [out1; out2; out3];