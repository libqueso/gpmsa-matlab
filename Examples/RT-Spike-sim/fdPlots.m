function fdPlots(pout,pvec,plotnum,varargin)

% Set up parameter labels
 labs={'CT','CB','CD'};

 model=pout.model;
 data=pout.data;
 pvals=pout.pvals(pvec);
 nreal=length(pvals);

 pu=model.pu; p=model.p;

 doPlot(1:4)=0;
 if exist('plotnum'); doPlot(plotnum)=1; end

% process input arguments
 ngrid=21; subset=1:p;
 parseAssignVarargs({'ngrid','subset'});
 lsub=length(subset);

% Set up prediction grid
 grid=linspace(0,1,ngrid);

 if doPlot(1)
   % plot the betaU response in box plots
   figure(1); clf; colormap([0 0 0]);
   bu=[pvals.betaU]';
   ru=exp(-bu/4);
   complexity=[]; sparsity=[];
   for ii=1:pu;
     b=bu(:,(ii-1)*p+1:ii*p);
     complexity=[complexity sum(b,2)]; sparsity=[sparsity sum(b.^2,2)];
     r=ru(:,(ii-1)*p+1:ii*p);
     gPackSubplot(pu,1,ii,1);
     boxplot(r);
     ylab=['\rho_{w',num2str(ii),'k}'];
     ylabel(ylab,'FontSize',12); xlabel('k');
     if ii==pu, ylim([0 1.01]); else ylim([0.01 1.01]); end
     if ii==1
       text(1:p,ones([1 p])+0.08,labs,'FontSize',12,...
            'HorizontalAlignment','center');
     end
     text(0.2,0.5,['PC' int2str(ii)],'FontSize',12,...
          'HorizontalAlignment','center','Rotation',90);
   end
   figure(1); print -depsc2 fdRhoBox; close;
   save 'complexity' complexity '-ascii';
   save 'sparsity' sparsity '-ascii';
 end

 if doPlot(2)    
   h=[];
   ymean=pout.simData.orig.ymean; ysd=pout.simData.orig.ysd;
   ctr=0;
   if ctr,
     for ii=1:pout.model.m
       ysim(:,ii)=pout.simData.orig.y(:,ii)-ymean;
     end
   else
     ysim=pout.simData.orig.y;
   end
   tr=[min(pout.simData.orig.t) max(pout.simData.orig.t)]; 
   rtr=range(tr);

   xpred=rand(20,3); npred=size(xpred,1);
   pred=gPred(xpred,pvals,model,data,'etamod');
   pw=zeros(length(pvals),pu,npred);
   for ii=1:pu
     pw(:,ii,:)=pred.Myhat(:,(ii-1)*npred+1:ii*npred);
   end
   pm=reshape(squeeze(mean(pw,1)),pu,[]);
   r=(pout.simData.Ksim*pm)*ysd;

   etapm=zeros(size(pout.simData.Ksim,1),npred);
   etalb=zeros(size(pout.simData.Ksim,1),npred);
   etaub=zeros(size(pout.simData.Ksim,1),npred);
   for ii=1:npred
     ps=squeeze(pred.w(:,:,ii))';
     eta=(pout.simData.Ksim*ps)*ysd;
     etalb(:,ii)=prctile(eta,5,2); etaub(:,ii)=prctile(eta,95,2);
     if ~ctr, etapm(:,ii)=etapm(:,ii)+ymean;
              etalb(:,ii)=etalb(:,ii)+ymean;
              etaub(:,ii)=etaub(:,ii)+ymean;
     end
   end
   save 'etapredPM' etapm '-ascii';
   save 'etapredLB' etalb '-ascii';
   save 'etapredUB' etaub '-ascii';

   figure(2); clf;
   for ii=1:4
     for jj=1:5
       % plot
       h(5*(ii-1)+jj)=gPackSubplot(4,5,ii,jj);
       plot(pout.simData.orig.t,ysim,'y'); hold on;
       plot(pout.simData.orig.t,etalb(:,5*(ii-1)+jj),'g','LineWidth',1);
       hold on;
       plot(pout.simData.orig.t,etaub(:,5*(ii-1)+jj),'g','LineWidth',1);
       xlim([tr(1)-0.1*rtr tr(2)+0.1*rtr]);
       if jj==1, ylabel('spike \alpha','FontSize',12);
       else set(gca,'YtickLabel',''); end
       if ii==4 && jj==3, xlabel('Atwood number','FontSize',12); end
     end
   end
   axisNorm(h,'ymax');
   figure(2); print -depsc2 fdPreds; close;
 end

 if doPlot(3)
   h=[];
   ctr=0;
   % Now the realizations over each theta
   figure(3); clf; colormap('copper');
   npred=ngrid;
   AzEl=[45 55];
   np=min(lsub,4);
   for ii=1:lsub
     jj=subset(ii);
     des=ones(npred,p)*0.5;
     des(:,jj)=grid';
     tt=pout.simData.orig.xrange(jj)*grid+pout.simData.orig.xmin(jj);
     pred=gPred(des,pvals,model,data,'etamod');
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
   figure(3); print -depsc2 fdNomSens; close;
 end

 if doPlot(4)
   h=[];
   ctr=0;
   % Now the marginalizations over each theta
   me=pout.sens.tmef.m; npred=size(me,3);
   if ctr, meanmat=repmat(pout.simData.orig.ymean,[1 npred]);
   else meanmat=0; end
   figure(4); clf; colormap('copper');
   tdat=0:1.0/(npred-1):1.0;
   AzEl=[45 55];
   np=min(lsub,4); 
   for ii=1:lsub
     jj=subset(ii);
     tt=pout.simData.orig.xrange(jj)*tdat+pout.simData.orig.xmin(jj); 
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
   figure(4); print -depsc2 fdMeSens; close;
 end

end
