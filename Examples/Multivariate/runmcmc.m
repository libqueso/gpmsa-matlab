% commands to run the Multivariate example...
addpath('../../matlab')

% read data
dat=fd(1:5);

% initial set-up
optParms.priors.lamOs.params=[1000 1000];
optParms.lamVzGroup=1:5;
params=setupModel(dat.obsData,dat.simData,optParms);

% modifications to defaults
params.priors.lamVz.params=repmat([1 0.0001],params.model.lamVzGnum,1);
params.priors.lamWs.params=repmat([1 0.0001],params.model.pu,1);

% modifications to initial values
params.model.lamOs=1;

% step size
nburn=100; nlev=21;
params=stepsize(params,nburn,nlev);

% mcmc
nmcmc=5000;
pout=gpmmcmc(params,nmcmc,'step',1);
save pout pout;

nmcmc=nmcmc+nburn*nlev;
pvec=floor(linspace(nburn*nlev+1,nmcmc,500));
pout.pvec500=pvec;

% calibration parameters
tmp=[pout.pvals.theta]';
save 'theta' tmp '-ascii';

% plots
%load pout;
fdPlots(pout,pvec,1:2);

holdout;
fdPlots(pout,pvec,3);

pvec=round(linspace(nburn*nlev+1,nmcmc,30));
pout.pvec30=pvec;
save pout pout

fdPlots(pout,pvec,4);

% cross-validation
pvec=floor(linspace(nburn*nlev+1,nmcmc,50));
fdPlots(pout,pvec,5);
