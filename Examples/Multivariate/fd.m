function params=fd(out,varargin);

 doPlot=1; pcpct=1;
 parseAssignVarargs({'doPlot','pcpct'});

% read in sim data
 nt=length(out); ysim=textread('sim_outputs');
 runs=ysim(:,1); ysim=ysim(:,out+1); ysim=ysim';

% read in design
 design=textread('design.txt');
 design=design(runs,:); m=size(design,1);
% scale design
 xmin=min(design); xrange=max(design)-xmin;
 design=(design-repmat(xmin,m,1))./repmat(xrange,m,1);

% read in obs data
 obsdata=textread('obs_outputs');
 n=size(obsdata,1); % number of experiments
% x in first nx columns, field data in rest
 nx=2; x=obsdata(:,1:nx);
 for ii=1:n, yobs(ii).y=obsdata(ii,1+nx:end)'; end
% specify assumed observation errors here
 for ii=1:n, yobs(ii).Sigy=diag((0.04^2)*ones(length(yobs(ii).y),1)); end
% scale obs x
 for ii=1:nx, x(:,ii)=(x(:,ii)-xmin(ii))./xrange(ii); end

% summary stats from sims
 ysimmean=mean(ysim,2);
 ysimStd=ysim-repmat(ysimmean,1,m);
 ysimsd=std(ysimStd,0,2);
 ysimStd=ysimStd./repmat(ysimsd,1,m);

% interpolate to data grid and standardize experimental data
 for ii=1:n
   yobs(ii).ymean=ysimmean;
   yobs(ii).yStd=(yobs(ii).y-yobs(ii).ymean)./ysimsd;
   yobs(ii).sd=sqrt(diag(yobs(ii).Sigy));
   yobs(ii).Sigy=yobs(ii).Sigy./(repmat(ysimsd.^2,1,nt));
 end

% K basis
% compute on simulations
 [U,S,V]=svd(ysimStd,0);
 lam=diag(S).^2/sum(diag(S).^2); lam=cumsum(lam);
 pu=sum(lam<pcpct)+1;
 Ksim=U(:,1:pu)*S(1:pu,1:pu)./sqrt(m);
 if(doPlot)
   figure(1);
   plot((1:pu)',Ksim);
   xlabel('Index','FontSize',14);
   ylabel('Output','FontSize',14);
   figure(1); print -depsc2 fdPc; close;
 end

% interpolate K onto data grids
 for ii=1:n, yobs(ii).Kobs=Ksim; end

% D basis
% compute on sim and data grids
 pv=nt;
 Dsim=eye(pv);
 for ii=1:n, yobs(ii).Dobs=eye(pv); end

% plot sims and data
 if (doPlot)
   Xp=[Ksim Dsim];
   for ii=1:n
     X=[yobs(ii).Kobs yobs(ii).Dobs]; Lamy=inv(yobs(ii).Sigy);
     uvhat=inv(X'*Lamy*X+1e-4*eye(pu+pv))*X'*Lamy*yobs(ii).yStd;
     yhat{ii}=Xp*uvhat; yhat{ii}=yhat{ii}.*ysimsd+ysimmean;
   end
   figure(2);
   for ii=1:m, cysim(:,ii)=ysim(:,ii)-ysimmean; end
   plot((1:pu)',cysim,'y*'); hold on;
   xlabel('Index','FontSize',14);
   ylabel('Output','FontSize',14);
   for ii=1:n
     plot((1:pu)',yobs(ii).y-yobs(ii).ymean,'b.','MarkerSize',20); hold on;
     plot((1:pu)',yhat{ii}-ysimmean,'c.','MarkerSize',10); hold on;
   end
   xlim([1 pu]); a=[yobs.y];
   ylim([min([cysim(:);a(:)]) max([cysim(:);a(:)])]);
   figure(2); print -depsc2 fdDat; close;
 end
      
% record the data into the obsData and simData structures
% first simData
% required fields
% must include dummy x if no x parameters in design
 simData.x=design;
 simData.yStd=ysimStd;
 simData.Ksim=Ksim;
% extra fields: original data and transform stuff
 simData.orig.y=ysim;
 simData.orig.ymean=ysimmean;
 simData.orig.ysd=ysimsd;
 simData.orig.Dsim=Dsim;
 simData.orig.xmin=xmin;
 simData.orig.xrange=xrange;
% obsData
% set the x values
% use dummy x if no x parameters in design
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
   obsData(ii).orig.sd=yobs(ii).sd;
 end

% pack up and leave
 params.simData=simData;
 params.obsData=obsData;

end
