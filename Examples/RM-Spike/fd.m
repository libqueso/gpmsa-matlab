function params=fd(doPlot,varargin);

 if ~exist('doPlot'); doPlot=0; else; doPlot=1; end

 pcpct=0.99; nkern=11;
 parseAssignVarargs({'pcpct','nkern'});

% read in design
 design=textread('design.txt'); m=size(design,1);
% scale design
 xmin=min(design); xrange=max(design)-xmin;
 design=(design-repmat(xmin,m,1))./repmat(xrange,m,1);

% read in sim data
 simdata=textread('sim_outputs');
 ysim=simdata(:,2:end);
 tsim=simdata(:,1);
 
% read in obs data
 n=1; % number of experiments
 for ii=1:n
   inf=['obs_outputs' int2str(ii)]; 
   obsdata{ii}=textread(inf);
   obsdata{ii}=obsdata{ii}(1:end-2,:);
   ydat{ii}=obsdata{ii}(:,2);
   tdat{ii}=obsdata{ii}(:,1);
   Sigy{ii}=1.0^2 .* eye(length(ydat{ii}));
 end

% summary stats from sims
 ysimmean=mean(ysim,2);
 ysimStd=ysim-repmat(ysimmean,1,m);
 ysimsd=std(ysimStd(:));
 ysimStd=ysimStd/ysimsd;

% structures for field data
% interpolate to data grid and standardize experimental data
 for ii=1:n
   yobs(ii).y=ydat{ii}; yobs(ii).t=tdat{ii};
   yobs(ii).ymean=interp1(tsim,ysimmean,yobs(ii).t,'linear','extrap');
   yobs(ii).yStd=(yobs(ii).y-yobs(ii).ymean)/ysimsd;
%   yobs(ii).Sigy=Sigy{ii}./(ysimsd.^2);
   yobs(ii).Sigy=Sigy{ii};
 end

% K basis
% compute on simulations
 [U,S,V]=svd(ysimStd,0);
 lam=diag(S).^2/sum(diag(S).^2); lam=cumsum(lam);
 pu=sum(lam<pcpct)+1;
 Ksim=U(:,1:pu)*S(1:pu,1:pu)./sqrt(m);
 if(doPlot)
   figure(1);
   plot(tsim,Ksim);
   xlabel('Atwood number','FontSize',14);
   ylabel('spike \theta','FontSize',14);
   figure(1); print -depsc2 fdPc; close;
 end

% interpolate K onto data grids
 for ii=1:n
   yobs(ii).Kobs=zeros(length(yobs(ii).yStd),pu);
   for jj=1:pu
     yobs(ii).Kobs(:,jj)=interp1(tsim,Ksim(:,jj),yobs(ii).t,'linear','extrap');
   end
 end

% D basis
% lay it out, and record decomposition on sim and data grids
% kernel centers and widths
 Dgrid=linspace(min(tsim),max(tsim),nkern); Dwidth=Dgrid(2)-Dgrid(1); 
 Dgrid=[Dgrid(1)-Dwidth Dgrid Dgrid(nkern)+Dwidth];
 pv=length(Dgrid);
% compute the kernel function map, for each kernel
 Dsim=zeros(size(ysimStd,1),pv);
 for ii=1:n, yobs(ii).Dobs=zeros(length(yobs(ii).yStd),pv); end
 for jj=1:pv
% first the obs
   for ii=1:n
     yobs(ii).Dobs(:,jj)=normpdf(yobs(ii).t,Dgrid(jj),Dwidth);
   end
% now the sim
   Dsim(:,jj)=normpdf(tsim,Dgrid(jj),Dwidth);
 end
% normalize the D maps
 Dmax=max(max(Dsim*Dsim'));
 Dsim=Dsim/sqrt(Dmax);
 for ii=1:n; yobs(ii).Dobs=yobs(ii).Dobs/sqrt(Dmax); end

% plot sims and data
 if (doPlot)
   Xp=[Ksim Dsim];
   for ii=1:n
     X=[yobs(ii).Kobs yobs(ii).Dobs]; Lamy=inv(yobs(ii).Sigy);
     uvhat=inv(X'*Lamy*X+1e-4*eye(pu+pv))*X'*Lamy*yobs(ii).yStd;
     yhat{ii}=Xp*uvhat; yhat{ii}=yhat{ii}*ysimsd+ysimmean;
   end
   figure(2);
   for ii=1:m, cysim(:,ii)=ysim(:,ii)-ysimmean; end
   plot(tsim,cysim,'y'); hold on;
   xlabel('Atwood number','FontSize',14);
   ylabel('spike \theta','FontSize',14);
   for ii=1:n
     plot(tdat{ii},yobs(ii).y-yobs(ii).ymean,'bo'); hold on;
     plot(tsim,yhat{ii}-ysimmean,'c.','MarkerSize',20); hold on;
   end
   xlim([min(tsim) max(tsim)]);
   ylim([min([cysim(:);[yobs.y]]) max([cysim(:);[yobs.y]])]);
   figure(2); print -depsc2 fdDat; close;
 end
      
% record the data into the obsData and simData structures
% first simData
% required fields
% must include dummy x if no x parameters in design
 simData.x=[.5*ones([m 1]) design];
 simData.yStd=ysimStd;
 simData.Ksim=Ksim;
% extra fields: original data and transform stuff
 simData.orig.y=ysim;
 simData.orig.ymean=ysimmean;
 simData.orig.ysd=ysimsd;
 simData.orig.Dsim=Dsim;
 simData.orig.t=tsim;
 simData.orig.xmin=[0.5 xmin];
 simData.orig.xrange=[0 xrange];
% obsData
% set the x values
% use dummy x if no x parameters in design
 x=[.5];
 for ii=1:n
% required fields
   obsData(ii).x=x(ii,:);
   obsData(ii).yStd=yobs(ii).yStd;
   obsData(ii).Kobs=yobs(ii).Kobs;
   obsData(ii).Dobs=yobs(ii).Dobs;
   obsData(ii).Sigy=yobs(ii).Sigy; 
% extra fields
   obsData(ii).orig.y=yobs(ii).y;
   obsData(ii).orig.ymean=yobs(ii).ymean;
   obsData(ii).orig.t=yobs(ii).t;
 end

% pack up and leave
 params.simData=simData;
 params.obsData=obsData;

end
