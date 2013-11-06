% commands to run the tower example
addpath('../../matlab')

% read data
towerdat = readdata();

% initial set-up
pout = setupModel(towerdat.obsData, towerdat.simData);

% do burn-in and step size estimation
nburn = 100; 
nlev = 13;
pout = stepsize(pout, nburn, nlev);

% get some mcmc draws that we're actually going to use
nmcmc = 10000;
pout = gpmmcmc(pout, nmcmc, 'step', 1);
nmcmc = nmcmc + nburn * nlev;
save pout;
