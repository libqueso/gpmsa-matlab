% commands to run the Flyer plate example...
addpath('../../matlab')

% read data
dat=fd(1,'pcpct',0.95,'nkern',29);

% initial set-up
params=setupModel(dat.obsData,dat.simData);

% modifications to defaults
params.priors.lamVz.params=repmat([1 0.0001],params.model.lamVzGnum,1);
params.priors.lamWs.params=repmat([1 0.0001],params.model.pu,1);

% user-specified theta prior
params.priors.theta.fname = 'gThetaPrior';
tmp=params.priors.theta.params;
params.priors.theta.params=[];
params.priors.theta.params.default=tmp; tmp=[];
params.priors.theta.params.mean = textread('c.theta');
params.priors.theta.params.chCov = chol(textread('c.cov'));
params.priors.theta.params.ind = 2:8;
p=params.model.p;
params.priors.theta.params.xmin=params.simData.orig.xmin(p+ ...
                                (params.priors.theta.params.ind));
params.priors.theta.params.xrange=params.simData.orig.xrange(p+ ...
                                  (params.priors.theta.params.ind));

% step size
nburn=100; nlev=21;
params=stepsize(params,nburn,nlev);

% mcmc
nmcmc=5000;
pout=gpmmcmc(params,nmcmc,'step',1);
save pout pout;

nmcmc=nmcmc+nburn*nlev;
pvec=floor(linspace(nburn*nlev+1,nmcmc,500));
pout.pvec=pvec;

% plots
%load pout;
fdPlots(pout,pvec,1:4);

% calibration parameters
tmp = [pout.pvals.theta]';
save 'theta' tmp '-ascii';

% sensitivity analysis
q=pout.model.q;
rn=zeros(p+q,2); rn(:,2)=1; rn(1:p,:)=.5;
sens=gSens(pout,'pvec',pvec,'varlist','all','rg',rn);

pout.sens=sens;
save pout pout;

tmp=[sens.smePm;sens.stePm]';
save 'sa_pm' tmp '-ascii';
tmp=sens.siePm';
save 'sa_ie' tmp '-ascii';

fdPlots(pout,pvec,5);
