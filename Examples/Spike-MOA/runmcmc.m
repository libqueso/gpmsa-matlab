% commands to run the Spike mix example...
addpath('../../matlab')

% file names
fname={'RT','RM'};
ndat=length(fname);

% read data
for ii=1:ndat
  tmp(ii)=fd(1,'pcpct',0.99975,'fname',fname{ii});
end

% initial set-up
for ii=1:ndat
  params(ii)=setupModel(tmp(ii).obsData,tmp(ii).simData);
end

% initial values for MCMC from separate calibration runs
% run makeW.m for each separate analysis and move
% output to MOA run directory
for ii=1:ndat
  w=[]; load(['w' fname{ii}]);

  params(ii).priors.lamVz.params=repmat([1 0.0001],params(ii).model.lamVzGnum,1);
  params(ii).priors.lamWs.params=repmat([1 0.0001],params(ii).model.pu,1);

  params(ii).model.theta=w.model.theta;
  params(ii).model.betaU=w.model.betaU;
  params(ii).model.betaV=w.model.betaV;
  params(ii).model.lamVz=w.model.lamVz;
  params(ii).model.lamUz=w.model.lamUz;
  params(ii).model.lamWs=w.model.lamWs;
  params(ii).model.lamWOs=w.model.lamWOs;
  params(ii).model.lamOs=w.model.lamOs;

  params(ii).mcmc.thetawidth=w.mcmc.thetawidth;
  params(ii).mcmc.rhoUwidth=w.mcmc.rhoUwidth;
  params(ii).mcmc.rhoVwidth=w.mcmc.rhoVwidth;
  params(ii).mcmc.lamVzwidth=w.mcmc.lamVzwidth;
  params(ii).mcmc.lamUzwidth=w.mcmc.lamUzwidth;
  params(ii).mcmc.lamWswidth=w.mcmc.lamWswidth;
  params(ii).mcmc.lamWOswidth=w.mcmc.lamWOswidth;
  params(ii).mcmc.lamOswidth=w.mcmc.lamOswidth;
end

% clist for MOA
% specify only common thetas across models
clist=[1 1;2 2;3 3];

% step size
nburn=100; nlev=31;
params=stepsize(params,nburn,nlev,'clist',clist);

% mcmc
nmcmc=5000;
pout=gpmmcmc(params,nmcmc,'step',1,'clist',clist);
save pout pout;

nmcmc=nmcmc+nburn*nlev;
pvec=floor(linspace(nburn*nlev+1,nmcmc,500));
for ii=1:ndat, pout(ii).pvec=pvec; end
save pout pout;

% plots
%load pout;
for ii=1:ndat
  fdPlots(pout(ii),pvec,1:3,'fname',fname{ii});

  % calibration parameters
  tmp = [pout(ii).pvals.theta]';
  save(['theta' fname{ii}],'tmp','-ascii');
end
