function fdPlots(pout,pvec,plotnum,varargin)

% Set up parameter labels
 xlabs={};
 thlabs={'CT','CB','CD'};
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
 ngrid=21; subset=1:nlabs; fname=[];
 parseAssignVarargs({'ngrid','subset','fname'});
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
     text(0.2,0.5,['PC' int2str(ii)],'FontSize',12,...
          'HorizontalAlignment','center','Rotation',90);
   end
   figure(1); print('-depsc2',['fdRhoBox' fname]); close;
   save(['complexity' fname],'complexity','-ascii');
   save(['sparsity' fname],'sparsity','-ascii');
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
   figure(2); print('-depsc2',['fdPost' fname]); close;
 end    

 if doPlot(3)    
   h=[];
   ctr=0;
   if ctr,
     for ii=1:pout.model.m
       ysim(:,ii)=pout.simData.orig.y(:,ii)-pout.simData.orig.ymean;
     end
     for ii=1:pout.model.n
       yobs{ii}=pout.obsData(ii).orig.y-pout.obsData(ii).orig.ymean;
     end
   else
     ysim=pout.simData.orig.y;
     for ii=1:pout.model.n, yobs{ii}=pout.obsData(ii).orig.y; end
   end
   reps=0;
   if reps
     ind=1;
     nind=length(ind);
   else ind=1:pout.model.n; nind=pout.model.n; end
   tr=[min(pout.simData.orig.t) max(pout.simData.orig.t)]; 
   rtr=range(tr);
   for ii=1:nind
     jj=ind(ii);
     figure(3); clf;
     outF=strcat('fdPreds',int2str(ii),fname);

     % calibrated eta
     pred=gPred(pout.obsData(jj).x,pvals,model,data,'uvpred');
     eta=pout.simData.Ksim*pred.u'.*pout.simData.orig.ysd;
     save(['etapred' fname],'eta','-ascii');
     etabounds=prctile(eta,[5 95],2);
     if ctr, meanmat=0;
     else meanmat=repmat(pout.simData.orig.ymean,[1 2]); end
 
     % now delta
     deltaR=pout.simData.orig.Dsim*pred.v'.*pout.simData.orig.ysd;
     save(['deltapred' fname],'deltaR','-ascii');
     deltaRbounds=prctile(deltaR,[5 95],2);
 
     % now a discrepancy-adjusted prediction
     yhat=deltaR+eta;
     save(['zetapred' fname],'yhat','-ascii');
     yhatbounds=prctile(yhat,[5 95],2);
 
     % plot
     h(1)=gPackSubplot(1,3,1,1);
     plot(pout.simData.orig.t,ysim,'y');
     hold on;
     if ii==nind, ll=pout.model.n; else ll=ind(ii+1)-1; end
     for kk=jj:ll, plot(pout.obsData(kk).orig.t,yobs{kk},'bo'); hold on; end
     plot(pout.simData.orig.t,etabounds+meanmat,'g','LineWidth',1);
     ylabel('spike \alpha','FontSize',12);
     title('calibrated simulator','Fontsize',12);
 
     h(2)=gPackSubplot(1,3,1,2);
     plot(pout.simData.orig.t,ysim,'y');
     hold on;
     for kk=jj:ll, plot(pout.obsData(kk).orig.t,yobs{kk},'bo'); hold on; end
     plot(pout.simData.orig.t,yhatbounds+meanmat,'k','LineWidth',1);
     set(gca,'YtickLabel','');
     xlabel('Atwood number','FontSize',12);
     title('discrepancy-adjusted','Fontsize',12);

     h(3)=gPackSubplot(1,3,1,3);
     plot(pout.simData.orig.t,deltaRbounds,'c','LineWidth',1); hold on;
     line(tr(1):rtr/100:tr(2),0,'LineStyle','--');
     set(gca,'YaxisLocation','right');
     title('discrepancy','Fontsize',12);

     axisNorm(h(1:2),'ymax'); axisNorm(h(3),'ymax');
     axisNorm(h(1:3),'x',[tr(1)-0.1*rtr tr(2)+0.1*rtr]);
     figure(3); print('-depsc2',outF); close;
   end
 end

 if doPlot(4)
   h=[];
   ctr=0;
   % Now the realizations over each theta
   figure(4); clf; colormap('copper');
   npred=ngrid;
   AzEl=[45 55];
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
     pm=reshape(squeeze(mean(pw,1)),pu,[]);
     r=(pout.simData.Ksim*pm)'*pout.simData.orig.ysd;
     if ~ctr
       for kk=1:npred
         r(kk,:)=r(kk,:)+pout.simData.orig.ymean';
       end
     end
     h(ii)=gPackSubplot(ceil(lsub/np),np,ceil(ii/np),mod(ii-1,np)+1,0.6);
     plot3(repmat(pout.simData.orig.t,size(grid)),...
           repmat(tt,size(pout.simData.orig.t)),r','b'); view(AzEl);
     xlabel('Atwood #'); ylabel(labs(jj)); zlabel('spike \alpha');
     alpha(0.25);
     axis tight;
     set(gca,'Xgrid','on','Ygrid','on','Zgrid','on');
   end
   axisNorm(h,'xzmax');
   figure(4); print('-depsc2',['fdNomSens' fname]); close;
 end

 if doPlot(5)
   h=[];
   ctr=0;
   % Now the marginalizations over each theta
   me=pout.sens.tmef.m; npred=size(me,3);
   if ctr, meanmat=repmat(pout.simData.orig.ymean,[1 npred]);
   else meanmat=0; end
   figure(5); clf; colormap('copper');
   tdat=0:1.0/(npred-1):1.0;
   AzEl=[45 55];
   np=min(lsub,4); 
   for ii=1:lsub
     jj=subset(ii);
     if nxlabs
       tt=pout.simData.orig.xrange(jj)*tdat+pout.simData.orig.xmin(jj); 
     else
       tt=pout.simData.orig.xrange(p+jj)*tdat+pout.simData.orig.xmin(p+jj);
     end
     r=squeeze(me(jj,:,:))-meanmat;
     h(ii)=gPackSubplot(ceil(lsub/np),np,ceil(ii/np),mod(ii-1,np)+1,0.6);
     plot3(repmat(pout.simData.orig.t,size(tdat)),...
           repmat(tt,size(pout.simData.orig.t)),r,'b'); view(AzEl);
     xlabel('Atwood #'); ylabel(labs(jj)); zlabel('spike \alpha');
     alpha(0.25);
     axis tight;
     set(gca,'Xgrid','on','Ygrid','on','Zgrid','on');
   end
   axisNorm(h,'xzmax');
   figure(5); print('-depsc2',['fdMeSens' fname]); close;
 end

end
