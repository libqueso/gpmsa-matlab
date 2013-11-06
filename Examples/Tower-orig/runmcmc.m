% commands to run the tower example...
addpath('../../matlab')

% read data
dat=sc(1);

% specify scalar output
optParms.scalarOutput=1;

% initial set-up
optParms.priors.lamOs.params=[1000 1000];
params=setupModel(dat.obsData,dat.simData,optParms);

% modifications to defaults
params.priors.lamVz.params=repmat([1 0.0001],params.model.lamVzGnum,1);
params.priors.lamWs.params=repmat([1 0.0001],params.model.pu,1);
params.priors.lamWOs.params=[2 1e-6];

% modifications to initial values
params.model.lamWOs=1e6;
params.model.lamOs=1;

% modifications to initial step widths
params.mcmc.lamWOswidth=0;

% step size
nburn=500; nlev=21;
params=stepsize(params,nburn,nlev);

% mcmc
nmcmc=10000;
pout=gpmmcmc(params,nmcmc,'step',1);
save pout pout;

nmcmc=nmcmc+nburn*nlev;
pvec=floor(linspace(nburn*nlev+1,nmcmc,500));
pout.pvec=pvec;

% plots
%load pout;
scPlots(pout,pvec,1:4);

% calibration parameters
tmp = [pout.pvals.theta]';
save 'theta' tmp '-ascii';

% sensitivity analysis
sens=gSens(pout,'pvec',pvec);

pout.sens=sens;
save pout pout;

tmp=[sens.smePm;sens.stePm]';
save 'sa_pm' tmp '-ascii';

scPlots(pout,pvec,5);
