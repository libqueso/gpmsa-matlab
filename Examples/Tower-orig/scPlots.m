function scPlots(pout,pvec,plotnum,varargin)

% Set up parameter labels
 xlabs={'x'};
 thlabs={'\theta'};
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

% Set up prediction grid
 grid=linspace(0,1,ngrid);

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
   end
   figure(1); print -depsc2 scRhoBox; close;
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
   figure(2); print -depsc2 scPost; close;
 end    

 if doPlot(3)    
   h=[];
   figure(3); clf;
   nx=length(pout.obsData);
   zc=norminv(.95);
   tr=[pout.simData.orig.xmin(1) pout.simData.orig.xmin(1)+ ...
                                 pout.simData.orig.xrange(1)];
   rtr=pout.simData.orig.xrange(1);
   zz=rtr*pout.data.z+tr(1);
   % calibrated eta and delta
   for ii=1:nx
     pred=gPred(pout.obsData(ii).x,pvals,model,data,'uvpred');
     eta(:,ii)=pout.simData.Ksim*pred.u'.*pout.simData.orig.ysd;
     deltaR(:,ii)=pout.simData.orig.Dsim*pred.v'.*pout.simData.orig.ysd;
   end
   save 'etapred' eta '-ascii';
   save 'deltapred' deltaR '-ascii';
   etabounds=prctile(eta,[5 95]);
   meanmat=repmat(pout.simData.orig.ymean,size(etabounds));
   peta=etabounds+meanmat;

   % now delta
   deltaRbounds=prctile(deltaR,[5 95]);

   % now a discrepancy-adjusted prediction
   yhat=deltaR+eta;
   save 'zetapred' yhat '-ascii';
   yhatbounds=prctile(yhat,[5 95]);
   pyhat=yhatbounds+meanmat;

 % plot
   h(1)=gPackSubplot(1,3,1,1);
   plot(zz,pout.simData.orig.y,'y*','MarkerSize',4); hold on;
   for ii=1:nx
     xx=rtr*pout.obsData(ii).x+tr(1);
     sd=pout.obsData(ii).orig.sd;
     plot(xx,pout.obsData(ii).orig.y,'b.','MarkerSize',4); hold on;
     plot([xx xx],pout.obsData(ii).orig.y+zc*[-sd sd],'k','LineWidth',2);
     hold on;
     plot([xx xx],peta(:,ii)','g');
     if ii<nx, hold on; end
   end
   ylabel('time','FontSize',12);
   title('calibrated simulator','Fontsize',12);

   h(2)=gPackSubplot(1,3,1,2);
   plot(zz,pout.simData.orig.y,'y*','MarkerSize',4); hold on;
   for ii=1:nx
     xx=rtr*pout.obsData(ii).x+tr(1);
     sd=pout.obsData(ii).orig.sd;
     plot(xx,pout.obsData(ii).orig.y,'b.','MarkerSize',4); hold on;
     plot([xx xx],pout.obsData(ii).orig.y+zc*[-sd sd],'k','LineWidth',2);
     hold on;
     plot([xx xx],pyhat(:,ii)','m'); 
     if ii<nx, hold on; end
   end
   set(gca,'YtickLabel','');
   xlabel('x','FontSize',12);
   title('discrepancy-adjusted','Fontsize',12);

   h(3)=gPackSubplot(1,3,1,3);
   for ii=1:nx
     xx=rtr*pout.obsData(ii).x+tr(1);
     plot([xx xx],deltaRbounds(:,ii)','c','LineWidth',2);
     if ii<nx, hold on; end
   end
   line(tr(1):rtr/100:tr(2),0,'LineStyle','--');
   set(gca,'YaxisLocation','right');
   title('discrepancy','Fontsize',12);

   axisNorm(h(1:2),'ymax'); axisNorm(h(3),'ymax');
   axisNorm(h,'x',[tr(1)-0.1*rtr tr(2)+0.1*rtr]);
   figure(3); print -depsc2 scPreds; close;
 end

 if doPlot(4)
   h=[];
   % Now the realizations over each parameter
   figure(4); clf; colormap('copper');
   Ksim=pout.simData.Ksim; ysd=pout.simData.orig.ysd;
   npred=ngrid;
   np=min(lsub,4);
   for ii=1:lsub
     jj=subset(ii);  
     des=ones(npred,p+q)*0.5;
     if nxlabs
       des(:,jj)=grid';
       tt=pout.simData.orig.xrange(jj)*grid+pout.simData.orig.xmin(jj);
     else
       des(:,p+jj)=grid';
       tt=pout.simData.orig.xrange(p+jj)*grid+pout.simData.orig.xmin(p+jj);
     end
     xpred=des(:,1:p); theta=des(:,1+p:end);
     pred=gPred(xpred,pvals,model,data,'wpred',theta);
     pw=zeros(length(pvals),pu,npred);
     for kk=1:pu
       pw(:,kk,:)=pred.Myhat(:,(kk-1)*npred+1:kk*npred);
     end
     my=squeeze(pw); pm=mean(my,1);
     r=diag(ysd)*(Ksim*pm);
     for kk=1:npred
       r(:,kk)=r(:,kk)+pout.simData.orig.ymean;
     end
     pvar=zeros(length(pvals),npred,pu,pu);
     for kk=1:length(pvals)
       for ll=1:npred
         pvar(kk,ll,:,:)=pred.Syhat{kk}(ll:npred:(pu-1)*npred+ll,...
                                        ll:npred:(pu-1)*npred+ll);
       end
     end
     for kk=1:length(pvals)
       for ll=1:npred
         svar=diag(diag(ysd)*Ksim*squeeze(pvar(kk,ll,:,:))*Ksim'*diag(ysd));
         sy(kk,ll)=sqrt(svar(1,1));
       end
     end
     v=mean(sy.^2+my.^2,1)-pm.^2; s=sqrt(v);

     h(ii)=gPackSubplot(ceil(lsub/np),np,ceil(ii/np),mod(ii-1,np)+1,0.5);
     plot(tt,r,'b'); hold on;
     plot(tt,r-2*s,'b--'); hold on; plot(tt,r+2*s,'b--');
     xlim([min(tt) max(tt)]);
     xlabel(labs(jj)); if ~mod(ii-1,np), ylabel('time'); end
     set(gca,'Xgrid','on','Ygrid','on');
   end
   axisNorm(h,'ymax');
   figure(4); print -depsc2 scNomSens; close;
 end

 if doPlot(5)
   h=[];
   % Now the marginalizations over each parameter
   me=pout.sens.tmef.m; npred=size(me,2);
   sd=pout.sens.tmef.sd;
   figure(5); clf; colormap('copper');
   tdat=0:1.0/(npred-1):1.0;
   np=min(lsub,4);
   for ii=1:lsub
     jj=subset(ii);
     if nxlabs
       tt=pout.simData.orig.xrange(jj)*tdat+pout.simData.orig.xmin(jj); 
     else
       tt=pout.simData.orig.xrange(p+jj)*tdat+pout.simData.orig.xmin(p+jj);
     end
     r=squeeze(me(jj,:)); s=squeeze(sd(jj,:));
     h(ii)=gPackSubplot(ceil(lsub/np),np,ceil(ii/np),mod(ii-1,np)+1,0.5);
     plot(tt,r,'b'); hold on;
     plot(tt,r-2*s,'b--'); hold on; plot(tt,r+2*s,'b--');
     xlim([min(tt) max(tt)]);
     xlabel(labs(jj)); if ~mod(ii-1,np), ylabel('time'); end
     set(gca,'Xgrid','on','Ygrid','on');
   end
   axisNorm(h,'ymax');
   figure(5); print -depsc2 scMeSens; close;
 end

end
