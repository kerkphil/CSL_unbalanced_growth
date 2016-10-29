% Jaimovich and Rebelo model - stationary version
clear

% set model parameters
alf = .35;
bet = .98;
sig = 2.5;
theta = 1.4;
gam = .1;
del = .025;
g   = .01;
psi = 1;
omega = .0013;
rho = .95;

% set up parameter vector to pass to DSGE function file
param = [alf bet sig theta gam del g psi omega rho];
% set numerical parameters
nx = 3;
ny = 0;
nz = 1;
nobs = 1000;
randomerr = 1;
logX = 0;
EE = 1;

% generate discret support for epsilon to be used in Euler error calcs
ne = 100;  %number of elements in support
Eps = zeros(ne,1);
Cum = -.5/ne;
Phi = ones(ne,1)/ne;
for e = 1:ne
    Cum = Cum + Phi(e);
    Eps(e) = norminv(Cum,0,sig);
end

% generate a history of Z's
Z = zeros(nobs+2,nz);
if randomerr
    eps = randn(nobs+2,1)*sig;
else
    load simulationshocks.mat
    eps = simeps;
    [nobs,~] = size(simeps)-2;
end
for t=1:nobs+1
    Z(t+1,:) = Z(t,:)*rho + eps(t+1,:);
end

%  steady state linearization
tic;
Zbar = 0;
%  find steady state
guess = [10; .5; .5];
XYbar = LinApp_FindSS(@JRstat_dyn,param,guess,Zbar,nx,ny)
if ny==0
    Xbar = XYbar;
    Ybar = [];
else
    Xbar = XYbar(1:nx);
    Ybar = XYbar(nx+1,nx+ny);
end

% set starting values
X0 = Xbar;

theta0 = [XYbar; XYbar; XYbar; 0; 0];
%  find derivatives
[AA, BB, CC, DD, FF, GG, HH, JJ, KK, LL, MM, WWW, TT] = ...
    LinApp_Deriv(@JRstat_dyn,param,theta0,nx,ny,nz,logX);
%  find policy function coefficients
[PP, QQ, UU, RR, SS, VVV] = ...
    LinApp_Solve(AA,BB,CC,DD,FF,GG,HH,JJ,KK,LL,MM,WWW,TT,rho,Zbar);
%  simulate
[XSSL, YSSL, EESSL] = LinApp_SSL_Euler(XYbar',Z,XYbar',logX,EE,Eps,Phi,...
                                      PP,QQ,UU,rho,@JRstat_dyn,param);
%  recover specific variables
kSSL  = XSSL(:,1);
hSSL  = XSSL(:,2);
xSSL  = XSSL(:,3);
YSSL  = zeros(nobs+2,1);
rSSL  = zeros(nobs+2,1);
wSSL  = zeros(nobs+2,1);
cSSL  = zeros(nobs+2,1);
iSSL  = zeros(nobs+2,1);
for t=2:nobs+1
    [YSSL(t+1), rSSL(t+1), wSSL(t+1), cSSL(t+1), iSSL(t+1)] = ...
     JRstat_defs(kSSL(t), hSSL(t), xSSL(t), Z(t), ...
     kSSL(t+1), hSSL(t+1), xSSL(t+1), param);
end
toc
% cacluate Euler error measures
AvgEESSL = mean(EESSL(2:nobs+1,:))
MaxAEESSL = max(abs(EESSL(2:nobs+1,:)))
RMSEESSL = sqrt(mean(EESSL(2:nobs+1,:).^2))
                                  
%  current state linarization
tic;
%  simulate
[XCSL, YCSL, EECSL] = LinApp_CSL_Euler(@JRstat_dyn,param,X0',Z,...
    rho,logX,EE,Eps,Phi);
%  recover specific variables
kCSL  = XCSL(:,1);
hCSL  = XCSL(:,2);
xCSL  = XCSL(:,3);
YCSL  = zeros(nobs+2,1);
rCSL  = zeros(nobs+2,1);
wCSL  = zeros(nobs+2,1);
cCSL  = zeros(nobs+2,1);
iCSL  = zeros(nobs+2,1);
for t=2:nobs+1
    [YCSL(t+1), rCSL(t+1), wCSL(t+1), cCSL(t+1), iCSL(t+1)] = ...
     JRstat_defs(kCSL(t), hCSL(t), xCSL(t), t, Z(t), ...
     kCSL(t+1), hCSL(t+1), xCSL(t+1), t+1, param);
end
toc
% cacluate Euler error measures
AvgEECSL = mean(EECSL(2:nobs+1,:))
MaxAEECSL = max(abs(EECSL(2:nobs+1,:)))
RMSEECSL = sqrt(mean(EECSL(2:nobs+1,:).^2))
