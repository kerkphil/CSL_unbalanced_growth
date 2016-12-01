function [Y, r, w, c, i, Gam, Del, Phi, Pi, Lam] = ...
    JRstat_defs(km, hm, xm, z, k, h, x, param)

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

Y = km^alf*(h*exp(z))^(1-alf);
r = alf*Y/km;
w = (1-alf)*Y/h;
c = w*h + (1+r-del)*km - k*exp(g);
i = k*exp(g) - (1-del)*km;
Gam = (c - psi*h^theta*c^gam*(xm*exp(-g))^(1-gam))^(-sig);
Lam = gam*c^(gam-1)*(xm*exp(-g))^(1-gam);
Del = 1 - psi*h^theta*Lam;
Phi = psi*h^theta*c^(gam-1)*(1-gam)*(xm*exp(-g))^(-gam);
Pi = psi*theta*h^(theta-1)*c^gam*(xm*exp(-g))^(1-gam);
