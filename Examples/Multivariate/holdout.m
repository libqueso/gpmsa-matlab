% commands to run the Multivariate holdouts...
%addpath('../../matlab')
%load pout

% holdout subsets
nx=length(pout.obsData);
rep=4; isub=setxor(1:nx,rep);

ny=length(pout.obsData(1).orig.y);
a=[nx ny 7];
peta=zeros(a); pdelta=zeros(a); pzeta=zeros(a); pyhat=zeros(a);

% holdout subsets

for ii=isub
   paramsh=[];

   % read data
   dat=fd(1:5,'doPlot',0);
   
   % make variance large for holdouts
   dat.obsData(ii).Sigy=1e6*dat.obsData(ii).Sigy;
   if ii==3, dat.obsData(ii+1).Sigy=1e6*dat.obsData(ii+1).Sigy; end

   % initial set-up
   optParms.priors.lamOs.params=[1000 1000];
   optParms.lamVzGroup=1:5;
   paramsh=setupModel(dat.obsData,dat.simData,optParms);

   % modifications to defaults
   paramsh.priors.lamVz.params=repmat([1 0.0001],paramsh.model.lamVzGnum,1);
   paramsh.priors.lamWs.params=repmat([1 0.0001],paramsh.model.pu,1);

   % modifications to initial values
   paramsh.model.lamOs=1;

   % step size
   nburn=100; nlev=21;
   paramsh=stepsize(paramsh,nburn,nlev);

   % mcmc
   nmcmc=500;
   paramsh=gpmmcmc(paramsh,nmcmc,'step',1);

   % prediction
   ymean=paramsh.simData.orig.ymean; ysd=paramsh.simData.orig.ysd;
   nmcmc=nmcmc+nburn*nlev;
   pvec=nburn*nlev+1:nmcmc; npvec=length(pvec);
   pred=gPred(paramsh.obsData(ii).x,paramsh.pvals(pvec),paramsh.model,...
              paramsh.data,'uvpred');
   eta=paramsh.simData.Ksim*pred.u';
   eta=eta.*repmat(ysd,[1 npvec]);
   meanmat=repmat(ymean,[1 npvec]);
   eta=eta+meanmat;
   peta(ii,:,:)=quantile(eta,[.01 .05 .25 .5 .75 .95 .99],2);
   delta=paramsh.simData.orig.Dsim*pred.v';
   delta=delta.*repmat(ysd,[1 npvec]);
   pdelta(ii,:,:)=quantile(delta,[.01 .05 .25 .5 .75 .95 .99],2);
   zeta=eta+delta;
   pzeta(ii,:,:)=quantile(zeta,[.01 .05 .25 .5 .75 .95 .99],2);
   yhat=zeta+mvnrnd(zeros(npvec,ny),...
                    pout.obsData(ii).Sigy.*repmat(ysd.^2,1,ny))'./...
                    sqrt(repmat([paramsh.pvals(pvec).lamOs],[ny 1]));
   pyhat(ii,:,:)=quantile(yhat,[.01 .05 .25 .5 .75 .95 .99],2);
   if ii==3
     peta(ii+1,:,:)=peta(ii,:,:);
     pdelta(ii+1,:,:)=pdelta(ii,:,:);
     pzeta(ii+1,:,:)=pzeta(ii,:,:);
     pyhat(ii+1,:,:)=pyhat(ii,:,:);
   end
end

pout.holdout.eta=peta; pout.holdout.delta=pdelta;
pout.holdout.zeta=pzeta; pout.holdout.yhat=pyhat;

save pout pout
