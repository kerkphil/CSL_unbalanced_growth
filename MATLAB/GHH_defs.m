function [Y, r, w, c, i, Gam, Del, Pi] = ...
    GHH_defs(km, hm, sm, z, k, h, s, param)

alf = param(1);
bet = param(2);
sig = param(3);
theta = param(4);
del = param(5);
g = param(6);
psi   = param(7); 
omega = param(8);
rho = param(9);

Y = km^alf*(h*exp(g*s+z))^(1-alf);
r = alf*Y/km;
w = (1-alf)*Y/h;
c = w*h + (1+r-del)*km - k;
i = k - (1-del)*km;
Gam = (c - psi*h^theta)^(-sig);
Del = 1 - psi*h^theta;
Pi = psi*theta*h^(theta-1);