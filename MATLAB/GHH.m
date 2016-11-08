% Greenwood, Hercowitz and Huffman model
clear

% set model parameters
alf = .35;
bet = .98;
sig = 2.5;
theta = 1.4;
del = .025;
g   = .01;
psi = 2;
omega = .013;
rho = .95;

% set starting values
X0 = [.01; .025; 0];

% set up parameter vector to pass to DSGE function file
param = [alf bet sig theta del g psi omega rho];

% set numerical parameters
nx = 3;
ny = 0;
nz = 1;
nobs = 300;
burn = 50;
nsim = 1;
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
    Eps(e) = norminv(Cum,0,omega);
end

Mom = zeros(4,1);
% begin simulations loop
for i = 1:nsim
    % generate a history of Z's
    Z = zeros(nobs+burn+2,nz);
    if randomerr
        eps = randn(nobs+burn+2,1)*omega;
    else
        load simulationshocks.mat
        eps = simeps;
        [nobs,~] = size(simeps)-2;
    end
    for t=1:nobs+burn+1
        Z(t+1,:) = Z(t,:)*rho + eps(t+1,:);
    end

    %  current state linarization
    tic;
    [XCSL, temp, EulerErr] = LinApp_CSL(@GHH_dyn,param,X0',Z,...
        rho,logX,EE,Eps,Phi);
    k  = XCSL(burn:nobs+burn+2,1);
    h  = XCSL(burn:nobs+burn+2,2);
    Y = zeros(nobs+2,1);
    r  = zeros(nobs+2,1);
    w  = zeros(nobs+2,1);
    c = zeros(nobs+2,1);
    i  = zeros(nobs+2,1);

    for t=2:nobs+1
        [Y(t+1), r(t+1), w(t+1), c(t+1), i(t+1)] = ...
         GHH_defs(k(t), h(t), t, Z(t), k(t+1), h(t+1), t+1, param);
    end
    Time = toc;

    AvgEE = mean(mean(EulerErr(2:nobs+1,:)));
    MaxAEE = max(max(abs(EulerErr(2:nobs+1,:))));
    RMSEE = sqrt(mean(mean(EulerErr(2:nobs+1,:).^2)));

    Mom(1) = ((nsim-1)/nsim)* Mom(1) + (1/nsim)*AvgEE;
    Mom(2) = max(Mom(2), MaxAEE);
    Mom(3) = sqrt(((nsim-1)/nsim)* Mom(1)^2 + (1/nsim)*RMSEE); 
    Mom(4) = Mom(4) + Time;
end

Mom

%  plot capital and hours
JR_Plot1(k(2:nobs+1,:),h(2:nobs+1,:))

% plot Euler errors for capital and hours 
JR_Plot1(EulerErr(2:nobs+1,1),EulerErr(2:nobs+1,2))