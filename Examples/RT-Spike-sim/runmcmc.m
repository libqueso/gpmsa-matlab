% commands to run the RT-Spike mix example...
addpath('../../matlab')

% read data
dat=fd(1,'pcpct',0.99975);

% initial set-up
params=setupModel(dat.obsData,dat.simData);

% modifications to defaults
params.priors.lamWs.params=repmat([1 0.0001],params.model.pu,1);

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
fdPlots(pout,pvec,1:3);

% sensitivity analysis
sens=gSens(pout,'pvec',pvec,'varlist','all');

pout.sens=sens;
save pout pout;

tmp=[sens.smePm;sens.stePm]';
save 'sa_pm' tmp '-ascii';
tmp=sens.siePm';
save 'sa_ie' tmp '-ascii';

fdPlots(pout,pvec,4);
