function params=fd(doPlot,varargin);

 if ~exist('doPlot'); doPlot=0; else; doPlot=1; end

 pcpct=0.99;
 parseAssignVarargs({'pcpct'});

% read in design
 design=textread('design.txt'); m=size(design,1);
% scale design
 xmin=min(design); xrange=max(design)-xmin;
 design=(design-repmat(xmin,m,1))./repmat(xrange,m,1);

% read in sim data
 simdata=textread('sim_outputs');
 ysim=simdata(:,2:end);
 tsim=simdata(:,1);

 if(doPlot)
   figure(1);
   plot(tsim,ysim,'y');
   xlabel('Atwood number','FontSize',14);
   ylabel('spike \alpha','FontSize',14);
   figure(1); print -depsc2 fdDat; close;
 end
 
% summary stats from sims
 ysimmean=mean(ysim,2);
 ysimStd=ysim-repmat(ysimmean,1,m);
 ysimsd=std(ysimStd(:));
 ysimStd=ysimStd/ysimsd;

% K basis
% compute on simulations
 [U,S,V]=svd(ysimStd,0);
 lam=diag(S).^2/sum(diag(S).^2); lam=cumsum(lam);
 pu=sum(lam<pcpct)+1;
 Ksim=U(:,1:pu)*S(1:pu,1:pu)./sqrt(m);
 if(doPlot)
   figure(2);
   plot(tsim,Ksim);
   xlabel('Atwood number','FontSize',14);
   ylabel('spike \alpha','FontSize',14);
   figure(2); print -depsc2 fdPc; close;
 end

% record the data into the obsData and simData structures
% first simData
% required fields
 simData.x=design;
 simData.yStd=ysimStd;
 simData.Ksim=Ksim;
% extra fields: original data and transform stuff
 simData.orig.y=ysim;
 simData.orig.ymean=ysimmean;
 simData.orig.ysd=ysimsd;
 simData.orig.t=tsim;
 simData.orig.xmin=xmin;
 simData.orig.xrange=xrange;
% obsData
 obsData=[];

% pack up and leave
 params.simData=simData;
 params.obsData=obsData;

end
