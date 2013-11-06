function scPlots(pout,pvec,plotnum,varargin)

% Set up parameter labels
 labs={'h','f','r','a'};

 model=pout.model;
 data=pout.data;
 pvals=pout.pvals(pvec);
 nreal=length(pvals);

 pu=model.pu; p=model.p;

 doPlot(1:5)=0;
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
   end
   figure(1); print -depsc2 scRhoBox; close;
   save 'complexity' complexity '-ascii';
   save 'sparsity' sparsity '-ascii';
 end

 if doPlot(2)    
   h=[];
   ymean=pout.simData.orig.ymean; ysd=pout.simData.orig.ysd;
   % uncalibrated response surface for the simulator
   grid=linspace(0,1,ngrid);
   [g1 g2]=meshgrid(grid,grid);
   g1v=g1(:); g2v=g2(:);
   gvlen=length(g1v);

   % do x_2, x_4 surface
   des=0.5*ones(gvlen,p);
   des(:,2)=g1v; des(:,4)=g2v;
   pred=gPred(des,pvals,model,data,'etamod');
   pm=squeeze(mean(pred.w,1));
   r=(pout.simData.Ksim*pm)*ysd;
   v=r+ymean;

   % make response surface
   tg1=g1*pout.simData.orig.xrange(2)+pout.simData.orig.xmin(2);
   tg2=g2*pout.simData.orig.xrange(4)+pout.simData.orig.xmin(4);
   figure(2);
   for ii=1:pu
     surf(tg1,tg2,reshape(v(:,ii),ngrid,ngrid)); colormap('summer');
     xlabel('x_2'); ylabel('x_4'); zlabel('time');
   end
   figure(2); print -depsc2 scPreds; close;
 end

 if doPlot(3)
   h=[];
   % Now the realizations over each parameter
   figure(3); clf; colormap('copper');
   Ksim=pout.simData.Ksim; ysd=pout.simData.orig.ysd;
   npred=ngrid;
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
   figure(3); print -depsc2 scNomSens; close;
 end

 if doPlot(4)
   h=[];
   % Now the marginalizations over each parameter
   me=pout.sens.tmef.m; npred=size(me,2);
   sd=pout.sens.tmef.sd;
   figure(4); clf; colormap('copper');
   tdat=0:1.0/(npred-1):1.0;
   np=min(lsub,4);
   for ii=1:lsub
     jj=subset(ii);
     tt=pout.simData.orig.xrange(jj)*tdat+pout.simData.orig.xmin(jj); 
     r=squeeze(me(jj,:)); s=squeeze(sd(jj,:));
     h(ii)=gPackSubplot(ceil(lsub/np),np,ceil(ii/np),mod(ii-1,np)+1,0.5);
     plot(tt,r,'b'); hold on;
     plot(tt,r-2*s,'b--'); hold on; plot(tt,r+2*s,'b--');
     xlim([min(tt) max(tt)]);
     xlabel(labs(jj)); if ~mod(ii-1,np), ylabel('time'); end
     set(gca,'Xgrid','on','Ygrid','on');
   end
   axisNorm(h,'ymax');
   figure(4); print -depsc2 scMeSens; close;
 end

 if doPlot(5)
   [h,cv]=gXval(pout,pvec,'figNum',5);
   figure(5); print -depsc2 scXval; close;
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
