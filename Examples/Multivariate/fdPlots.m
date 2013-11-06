function fdPlots(pout,pvec,plotnum,varargin)

% Set up parameter labels
 xlabs={'x_1','x_2'};
 thlabs={'\theta_1','\theta_2','\theta_3'};
 labs=[xlabs thlabs];
 nxlabs=length(xlabs); nlabs=length(labs);

 model=pout.model;
 data=pout.data;
 pvals=pout.pvals(pvec);
 nreal=length(pvals);

 pu=model.pu; pv=model.pv;
 p=model.p; q=model.q;

 doPlot(1:5)=0;
 if exist('plotnum'); doPlot(plotnum)=1; end

% process input arguments
 ngrid=21; subset=1:nlabs;
 parseAssignVarargs({'ngrid','subset'});
 lsub=length(subset);
 if nxlabs, thsub=subset(subset>p)-p; else thsub=subset; end

 if doPlot(1)
   % plot the betaU response in box plots
   figure(1); clf; colormap([0 0 0]);
   bu=[pvals.betaU]';
   ru=exp(-bu/4);
   complexity=[]; sparsity=[];
   for ii=1:pu;
     if nxlabs, b=bu(:,(ii-1)*(p+q)+1:ii*(p+q));
     else b=bu(:,(ii-1)*(p+q)+p+1:ii*(p+q)); end
     complexity=[complexity sum(b,2)]; sparsity=[sparsity sum(b.^2,2)];
     if nxlabs, r=ru(:,(ii-1)*(p+q)+1:ii*(p+q));
     else r=ru(:,(ii-1)*(p+q)+p+1:ii*(p+q)); end
     gPackSubplot(pu,1,ii,1);
     boxplot(r);
     ylab=['\rho_{w',num2str(ii),'k}'];
     ylabel(ylab,'FontSize',12); xlabel('k');
     if ii==pu, ylim([0 1.01]); else ylim([0.01 1.01]); end
     if ii==1
       text(1:nlabs,ones([1 nlabs])+0.08,labs,'FontSize',12,...
            'HorizontalAlignment','center');
     end
     text(0,0.5,['PC' int2str(ii)],'FontSize',12,...
          'HorizontalAlignment','center','Rotation',90);
   end
   figure(1); print -depsc2 fdRhoBox; close;
   save 'complexity' complexity '-ascii';
   save 'sparsity' sparsity '-ascii';
 end

 if doPlot(2)
   figure(2); clf;
   if length(pvec)>1000,
     pvec2=pvec(floor(linspace(1,length(pvec),1000)));
   else pvec2=pvec; end
   t=[pout.pvals(pvec2).theta]'; t=t(:,thsub);
   gPlotMatrix(t,'Pcontours',[0.5 0.9], ...
               'ustyle','imcont','lstyle','imcont', ...
               'ngrid',ngrid,'ksd',0.1,'labels',thlabs(thsub));
   figure(2); print -depsc2 fdPost; close;
 end    

 if doPlot(3)    
   h=[]; holdout=1;
   figure(3); clf;
   colormap jet;

   % predictions using all data
   ymean=pout.simData.orig.ymean; ysd=pout.simData.orig.ysd;
   a=[length(pout.obsData),length(pout.obsData(1).orig.y),nreal];
   peta=zeros(a); pdelta=zeros(a); pzeta=zeros(a); pyhat=zeros(a);
   xx=reshape([pout.obsData.x],2,[])';
   xx=reshape(repmat(xx,1,a(2))',2,[])';
   for ii=1:length(pout.obsData)
     ny=length(pout.obsData(ii).orig.y);
     pred=gPred(pout.obsData(ii).x,pvals,model,data,'uvpred');
     eta=pout.simData.Ksim*pred.u';
     eta=eta.*repmat(ysd,1,nreal);
     meanmat=repmat(ymean,1,nreal);
     peta(ii,:,:)=eta+meanmat;
     delta=pout.simData.orig.Dsim*pred.v';
     pdelta(ii,:,:)=delta.*repmat(ysd,1,nreal);
     pzeta(ii,:,:)=pdelta(ii,:,:)+peta(ii,:,:);
     pyhat(ii,:,:)=squeeze(pzeta(ii,:,:))+mvnrnd(zeros(nreal,ny),...
                   pout.obsData(ii).Sigy.*repmat(ysd.^2,1,ny))'./...
                   sqrt(repmat([pvals.lamOs],[ny 1]));
   end
   peta=quantile(peta,[.01 .05 .25 .5 .75 .95 .99],3); eta=peta;
   a=size(peta);
   peta=[xx reshape(permute(peta,[2 1 3]),a(1)*a(2),a(3))];
   save 'etapred' peta '-ascii';
   pdelta=quantile(pdelta,[.01 .05 .25 .5 .75 .95 .99],3);
   pdelta=[xx reshape(permute(pdelta,[2 1 3]),a(1)*a(2),a(3))];
   save 'deltapred' pdelta '-ascii';
   pzeta=quantile(pzeta,[.01 .05 .25 .5 .75 .95 .99],3); zeta=pzeta;
   pzeta=[xx reshape(permute(pzeta,[2 1 3]),a(1)*a(2),a(3))];
   save 'zetapred' pzeta '-ascii';
   pyhat=quantile(pyhat,[.01 .05 .25 .5 .75 .95 .99],3); yhat=pyhat;
   pyhat=[xx reshape(permute(pyhat,[2 1 3]),a(1)*a(2),a(3))];
   save 'yhatpred' pyhat '-ascii';

   % make multivariate prediction plot
   if holdout, heta=pout.holdout.eta; hzeta=pout.holdout.zeta;
               hyhat=pout.holdout.yhat; end
   ny=length(pout.obsData(1).orig.y); nx=length(pout.obsData);
   if holdout, xx=[0.05 0.25 0.45 0.7 0.95]; else xx=0.1:0.2:0.9; end
   tc=norminv(.95);
   kk=0;
   for ii=1:ny
     for jj=1:nx
       kk=kk+1; 
       h(kk)=gPackSubplot(ny,nx,ii,jj);
       xxx=[xx(1)-.05;xx(1)+.05];
       for ll=1:model.m
         yyy=pout.simData.orig.y(ii,ll).*ones(2,1);
         plot(xxx,yyy,'g-'); hold on;
       end
       xxx=[xx(2)-.05;xx(2)+.05];
       yyy=pout.obsData(jj).orig.y(ii).*ones(2,1);
       plot(xxx,yyy,'k-'); hold on;
       xxx=xx(2).*ones(2,1);
       yyy(1)=yyy(1)-tc*pout.obsData(jj).orig.sd(ii);
       yyy(2)=yyy(2)+tc*pout.obsData(jj).orig.sd(ii);
       plot(xxx,yyy,'k-','LineWidth',6); hold on;
       xxx=[xx(3)-.05;xx(3)+.05];
       yyy=eta(jj,ii,4).*ones(2,1);
       plot(xxx,yyy,'b-'); hold on;
       xxx=xx(3).*ones(2,1);
       yyy(1)=eta(jj,ii,2); yyy(2)=eta(jj,ii,6);
       plot(xxx,yyy,'b-','LineWidth',6); hold on;
       xxx=[xx(4)-.05;xx(4)+.05];
       if holdout, xxx=xxx-.05; end
       yyy=zeta(jj,ii,4).*ones(2,1);
       plot(xxx,yyy,'r-'); hold on;
       xxx=xx(4).*ones(2,1); if holdout, xxx=xxx-.05; end
       yyy(1)=zeta(jj,ii,2); yyy(2)=zeta(jj,ii,6);
       plot(xxx,yyy,'r-','LineWidth',6); hold on;
       if holdout
         xxx=[xx(4)-.05;xx(4)+.05]+.05;
         yyy=hzeta(jj,ii,4).*ones(2,1);
         plot(xxx,yyy,'-','Color',[0.5 0 0]); hold on;
         xxx=xx(4).*ones(2,1)+.05;
         yyy(1)=hzeta(jj,ii,2); yyy(2)=hzeta(jj,ii,6);
         plot(xxx,yyy,'-','Color',[0.5 0 0],'LineWidth',6); hold on;
       end
       xxx=[xx(5)-.05;xx(5)+.05];
       if holdout, xxx=xxx-.05; end
       yyy=yhat(jj,ii,4).*ones(2,1);
       plot(xxx,yyy,'y-'); hold on;
       xxx=xx(5).*ones(2,1); if holdout, xxx=xxx-.05; end
       yyy(1)=yhat(jj,ii,2); yyy(2)=yhat(jj,ii,6);
       plot(xxx,yyy,'y-','LineWidth',6);
       if holdout
         hold on;
         xxx=[xx(5)-.05;xx(5)+.05]+.05;
         yyy=hyhat(jj,ii,4).*ones(2,1);
         plot(xxx,yyy,'-','Color',[1 .5 0]); hold on;
         xxx=xx(5).*ones(2,1)+.05;
         yyy(1)=hyhat(jj,ii,2); yyy(2)=hyhat(jj,ii,6);
         plot(xxx,yyy,'-','Color',[1 .5 0],'LineWidth',6);
       end
       v=axis;
       ylb=linspace(v(3),v(4),5); ylb=ylb(2:4);
       if (jj==1),
          ylabel(strcat('y_',int2str(ii)));
          set(gca,'Ytick',ylb);
       else
          set(gca,'YtickLabel','');
       end
       if (ii==1)
          xlabel(strcat('Expt. ',int2str(jj)));
          set(gca,'XaxisLocation','top');
       end
       if holdout, set(gca,'XLim',[-0.05 1.1]); else
                   set(gca,'XLim',[0 1]); end
       set(gca,'XtickLabel','');
     end
   end
   kk=1:nx*ny; kk=reshape(kk,nx,[])';
   for ii=1:ny, axisNorm(h(kk(ii,:)),'xymax'); end
   figure(3); print -depsc2 fdPreds; close;
 end

 if doPlot(4)
   h=[];
   ymean=pout.simData.orig.ymean; ysd=pout.simData.orig.ysd;
   % uncalibrated response surface for the simulator
   ngrid=11;
   grid=linspace(0,1,ngrid);
   [g1 g2]=meshgrid(grid,grid);
   g1v=g1(:); g2v=g2(:);
   gvlen=length(g1v);

   % do some theta_1, theta_3 surfaces
   theta=0.5*ones(gvlen,q);
   theta(:,1)=g1v; theta(:,3)=g2v;
   pred=gPred(0.5*ones(gvlen,p),pvals,model,data,'wpred',theta);
   pm=squeeze(mean(pred.w,1));
   r=(pout.simData.Ksim*pm)'.*repmat(ysd',gvlen,1);
   v=r+repmat(ymean',gvlen,1);

   % make response surfaces
   tg1=g1*pout.simData.orig.xrange(p+1)+pout.simData.orig.xmin(p+1);
   tg2=g2*pout.simData.orig.xrange(p+3)+pout.simData.orig.xmin(p+3);
   for ii=1:pu
     figure(30+ii);
     surf(tg1,tg2,reshape(v(:,ii),ngrid,ngrid)); colormap('summer');
     xlabel('\theta_1'); ylabel('\theta_3'); title(['y_' num2str(ii)]);
     zlim([min(pout.simData.orig.y(ii,:)) max(pout.simData.orig.y(ii,:))]);
   end

   % do some theta_1, theta_2 surfaces
   theta=0.5*ones(gvlen,q);
   theta(:,1)=g1v; theta(:,2)=g2v;
   pred=gPred(0.5*ones(gvlen,p),pvals,model,data,'wpred',theta);
   pm=squeeze(mean(pred.w,1));
   r=(pout.simData.Ksim*pm)'.*repmat(ysd',gvlen,1);
   v=r+repmat(ymean',gvlen,1);

   % make response surfaces
   tg1=g1*pout.simData.orig.xrange(p+1)+pout.simData.orig.xmin(p+1);
   tg2=g2*pout.simData.orig.xrange(p+2)+pout.simData.orig.xmin(p+2);
   for ii=1:pu
     figure(30+ii+pu);
     surf(tg1,tg2,reshape(v(:,ii),ngrid,ngrid)); colormap('summer');
     xlabel('\theta_1'); ylabel('\theta_2'); title(['y_' num2str(ii)]);
     zlim([min(pout.simData.orig.y(ii,:)) max(pout.simData.orig.y(ii,:))]);
   end

   % save some plots
   figure(37); print -depsc2 surf212.eps;
   figure(36); print -depsc2 surf112.eps;
   figure(33); print -depsc2 surf313.eps;
   figure(34); print -depsc2 surf413.eps;
   figure(35); print -depsc2 surf513.eps;
 end

 if doPlot(5)
   [h,cv]=gXval(pout,pvec,'figNum',5);
   figure(5); print -depsc2 fdXval; close;
   yobs=pout.simData.orig.y; Ksim=pout.simData.Ksim;
   ymean=pout.simData.orig.ymean; ysd=pout.simData.orig.ysd;
   ypred=zeros(size(yobs)); sdpred=zeros(size(yobs));
   for ii=1:model.m
     ypred(:,ii)=diag(ysd)*Ksim*cv.pred(ii,:)'+ymean;
     Sigma=diag(ysd)*Ksim*squeeze(cv.cov(ii,:,:))*Ksim'*diag(ysd);
     sdpred(:,ii)=sqrt(diag(Sigma));
   end
   abserr=100*abs(yobs-ypred)./abs(yobs);
   sderr=100*sdpred./abs(yobs);
   cverr=[];
   for ii=1:size(ypred,1)
     cverr(3*(ii-1)+(1:3),:)=[yobs(ii,:);abserr(ii,:);sderr(ii,:)];
   end
   save 'cverr' cverr '-ascii'
 end

end
