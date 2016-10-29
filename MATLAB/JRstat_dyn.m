function out = JRstat_dyn(in,param)

kp    = in(1);
hp    = in(2);
xp    = in(3);
k     = in(4);
h     = in(5);
x     = in(6);
km    = in(7);
hm    = in(8);
xm    = in(9);
zp    = in(10);
z     = in(11);

alf   = param(1);
bet   = param(2);
sig   = param(3);
theta = param(4);
gam   = param(5);
del   = param(6);
g     = param(7);
psi   = param(8); 
omega = param(9);
rho   = param(10);

[Y, r, w, c, i, Gam, Del, Phi, Pi, Lam] = ...
    JRstat_defs(km, hm, xm, z, k, h, x, param);
[Yp, rp, wp, cp, ip, Gamp, Delp, Phip, Pip, Lamp] = ...
    JRstat_defs(k, h, x, zp, kp, hp, xp, param);

out1 = bet*(Lam*Phip + Gamp*exp(g)*Delp*(1+rp-del))/(Gam*Del) - 1;
out2 = (Gam/Pi + bet*Lam*w*Gamp*exp(g)*Phip)/(Gam*Del*w) - 1;
out3 = c^gam*(xm*exp(-g))^(1-gam)/x - 1;
out = [out1; out2; out3];