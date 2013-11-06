% commands to run the smoke plume example...
addpath('../../matlab')

% read data
dat=sc(1);

% specify scalar output
optParms.scalarOutput=1;

% initial set-up
params=setupModel(dat.obsData,dat.simData,optParms);

% modifications to defaults
params.priors.lamWs.params=repmat([1 0.0001],params.model.pu,1);
params.priors.lamWOs.params=[2 1e-6];

% modifications to initial values
params.model.lamWOs=1e6;

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
scPlots(pout,pvec,1:3);

% sensitivity analysis
sens=gSens(pout,'pvec',pvec,'varlist','all');

pout.sens=sens;
save pout pout;

tmp=[sens.smePm;sens.stePm]';
save 'sa_pm' tmp '-ascii';
tmp=sens.siePm';
save 'sa_ie' tmp '-ascii'

scPlots(pout,pvec,4);

% cross-validation
pvec=floor(linspace(nburn*nlev+1,nmcmc,50));
scPlots(pout,pvec,5);
