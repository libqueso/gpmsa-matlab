function neddPlots(pout,pvec, plotnum, varargin)

newFig=1;
holdOutMods=[];
plotStyle='2d';
pout1=[];
parseAssignVarargs({'newFig','holdOutMods','plotStyle','pout1'});

simData=pout.simData; obsData=pout.obsData;
model=pout.model;
data=pout.data;
pvals=pout.pvals(pvec);
nreal=length(pvals);

pu=model.pu; pv=model.pv;
p=model.p; q=model.q;

doPlot(10)=0;
if exist('plotnum'); doPlot(plotnum)=1;
else
end


if doPlot(1) 
    % Plot the betaU response in box plots
    if newFig; figure(1); end
    clf; colormap([0 0 0]);
    %just for the thetas, so discard the (dummy) x columns
    bu=[pvals.betaU]';
    ru=exp(-0.25*bu);
    for ii=1:p+q;
        r=ru(:,(ii-1)*(p+q)+1:ii*(p+q));
        subplot(p+q,1,ii);
        boxplot(r);
        xlabel('[x \theta]');
        ylabel(['PC' num2str(ii)]);
        axisNorm(gca,'y',[0 0 0 1]);
    end
end 

if doPlot(2)
    % First, theta posteriors.
    if newFig; figure(2); end
    clf; %colormap('gray');
    t=[pout.pvals(min(pvec):max(pvec)).theta]';
    C=gPlotMatrix(t,'ngrid',50,'ustyle','imcont','ksd',0.04);

end    

if doPlot(3)
    % Now the realizations over each x,theta for sensitivity

    % Set up prediction grid
    switch(plotStyle)
      case '3d'
        pgrid=linspace(0,1,15);
      case '2d'
        pgrid=linspace(0,1,5); 
    end
    glen=length(pgrid);  mlen=glen*glen;  
    [gridx,gridy]=meshgrid(pgrid,pgrid);
    
    if newFig; figure(3); end;
    clf; colormap('copper');
    npred=glen;
    AzEl=[45 55];    
    for ii=1:p+q
        tdat=ones(npred,p+q)*0.5;
        tdat(:,ii)=pgrid';
        xpred=tdat(:,1:p);
        theta=tdat(:,p+1:end);
        pred=gPred(xpred,pvals,model,data,'wpred',theta);
        wrm=squeeze(mean(pred.w))';
        
        %r=(pout.simData.Keta*wrm')';
        %v=r*pout.simData.yesd + repmat(pout.simData.yemean(:)',npred,1);
        r=(pout.simData.Ksim*wrm')'*pout.simData.orig.ysd;
        v=r + repmat(pout.simData.orig.ymean(:)',npred,1);
        r=r(:,1:22);
        v=v(:,1:22);

        switch(plotStyle)
        case '3d'
          h1(ii)=subplot(p+q,3,ii*3-2);
            mesh(pout.simData.orig.time,pgrid,r,'MeshStyle','row'); view(AzEl);
            alpha(0.25)
          h2(ii)=subplot(p+q,3,ii*3-1);
            mesh(pout.simData.orig.time,pgrid,v,'MeshStyle','row'); view(AzEl);
            alpha(0.25)

        case '2d'
          h1(ii)=subplot(p+q,3,ii*3-2);
            hold on;
            for kk=length(pgrid):-1:1
              plot(pout.simData.orig.time,r(kk,:), ...
                   'color',[1 1 1]*0.75*(1 - kk/length(pgrid)),'linewidth',1.25 );
            end
          h2(ii)=subplot(p+q,3,ii*3-1);
            hold on;
            for kk=length(pgrid):-1:1
              plot(pout.simData.orig.time,v(kk,:), ...
                   'color',[1 1 1]*0.75*(1 - kk/length(pgrid)),'linewidth',1.25 );
            end
        end
        subplot(p+q,3,ii*3-2);
          if ii<=p; var=['x_' num2str(ii)];else; var=['\theta_' num2str(ii-p)];end
          xlabel('time'); ylabel(var); zlabel('std. radius');
          axis tight
        subplot(p+q,3,ii*3-1);
          xlabel('time'); ylabel(var); zlabel('radius');
          axis tight
       drawnow
    end
    
    switch(plotStyle)
      case '3d'
        axisNorm(h1,'xyzmax');
        axisNorm(h2,'xyzmax');
      case '2d'
        axisNorm(h1,'xymax');
        axisNorm(h2,'xymax');
    end
    drawnow
end

if doPlot(4)    
    AzEl=[140 30];
    if newFig; figure(4); end
    clf
    for ii=1:3
      h(ii)=subplot(3,3,ii);
        npred=1;
        pred=gPred(pout.obsData(ii).x(1),pvals,model,data,'uvpred');
        urm=mean(pred.u);
        vrm=mean(pred.v);
        
        yvmhat=(pout.obsData(1).orig.Dsim*vrm')';
        yumhat=(pout.simData.Ksim*urm')';
        rm=yvmhat+yumhat;
        ysd=pout.simData.orig.ysd; ymean=pout.simData.orig.ymean';
        vm=(yvmhat+yumhat)*ysd + ymean;
        
        plot3(pout.obsData(ii).orig.phi,pout.obsData(ii).orig.time, ...
              pout.obsData(ii).orig.y,...
              'bo','MarkerFaceColor','b','MarkerSize',3);
        view(AzEl);
        xlabel('\phi');ylabel('time');zlabel('r');
        set(gca,'XGrid','on','YGrid','on','ZGrid','on');
        hold on;
        
        colormap([0 0 0])
        mesh(pout.simData.orig.phimat,pout.simData.orig.timemat,...
             reshape(vm,size(pout.simData.orig.timemat)));
        hidden off;
        title(['Experiment ' num2str(ii)]);
        axis tight;
        drawnow;
    end
    axisNorm(h,'xyzmax');
end

if doPlot(5)
  if newFig; figure(5); end;
  clf
  % Plot the bounded preds in polar co-ords
  for jj=1:length(pout.obsData);
    npred=1;
    pred=gPred(pout.obsData(jj).x(1),pvals,model,data,'uvpred');
    predv=pred.v;
    predu=pred.u;
    
    pyd=(pout.obsData(jj).orig.Dsim*predv')';
    pye=(pout.simData.Ksim*predu')';
    py=pyd+pye;        
    
    pym=mean(py);
    pyh=prctile(py,95);
    pyl=prctile(py,5);
    
    ysd=pout.simData.orig.ysd; ymean=pout.simData.orig.ymean';
    vm=pym*ysd + ymean;
    vl=pyl*ysd + ymean;
    vh=pyh*ysd + ymean;
    
    vm=reshape(vm,size(pout.simData.orig.timemat));
    vl=reshape(vl,size(pout.simData.orig.timemat));
    vh=reshape(vh,size(pout.simData.orig.timemat));
    
    d1u=unique(pout.obsData(jj).orig.time);
    for ii=1:length(d1u)
      XTickOn=1;
      YTickOn=1;
      switch(jj)
      case 1
       switch(ii)
       case 1
        D2subplot(1,3,5,[.35 .08 .05 .02],[0 0 0 0])
        XTickOn=0;
        title('experiment 1')
       case 2
        D2subplot(6,3,5,[.35 .08 .05 .02],[0 0 0 0])
        XTickOn=0;
       case 3
        D2subplot(11,3,5,[.35 .08 .05 .02],[0 0 0 0])
       end
      case 2
        D2subplot(12,3,5,[.35 .08 .05 .02],[0 0 0 0])
        YTickOn=0;
        title('experiment 2')
      case 3
       YTickOn=0;
       switch(ii)
       case 1
        D2subplot(8,3,5,[.35 .08 .05 .02],[0 0 0 0])
        XTickOn=0;
        title('experiment 3')
       case 2
        D2subplot(13,3,5,[.35 .08 .05 .02],[0 0 0 0])
       end
      end
   
       dinds=pout.obsData(jj).orig.time==d1u(ii);
        hold on;    

        [x,y]=pol2cart(pout.obsData(jj).orig.phi(dinds),pout.obsData(jj).orig.y(dinds));
          plot(x,y,'o','MarkerFaceColor','b','MarkerSize',3);
          axis equal; axis([-1 1 -1 1]*3.8); 
        timeline=pout.simData.orig.timemat(:,1);
        tdist=timeline-d1u(ii);
        tdist1=tdist; tdist1(tdist1>0)=-Inf; [y1 jj1]=max(tdist1); t1=timeline(jj1);
        tdist2=tdist; tdist2(tdist2<0)= Inf; [y2 jj2]=min(tdist2); t2=timeline(jj2);
        c=(d1u(ii)-timeline(jj1))/(timeline(jj2)-timeline(jj1));
        vmp = (1-c)* vm(jj1,:) + c*vm(jj2,:); 
        vhp = (1-c)* vh(jj1,:) + c*vh(jj2,:); 
        vlp = (1-c)* vl(jj1,:) + c*vl(jj2,:); 
        vlp=max(vlp,0);
        %[x,y]=pol2cart(pout.simData.phimat(1,:),vmp); plot(x,y);     %mean
        [x,y]=pol2cart(pout.simData.orig.phimat(1,:),vlp); plot(x,y,'k'); %lower
        [x,y]=pol2cart(pout.simData.orig.phimat(1,:),vhp); plot(x,y,'k'); %upper
      if ~XTickOn; set(gca,'XTick',[]); else; xlabel('cm'); end
      if ~YTickOn; set(gca,'YTick',[]); else; ylabel('cm'); end
    end
    drawnow
  end
  axes('Position',[0 0 1 1],'Visible','off');
  text(.65,.85,'t = 15 \mus','FontSize',14);
  text(.65,.65,'t = 25 \mus','FontSize',14);
  text(.65,.45,'t = 45 \mus','FontSize',14);

end

if doPlot(6) % Holdout Plots
  if ~exist('poutHO'); error('poutHO must be passed in for holdout plots'); end
  if newFig; figure(6); end;
  clf
  % Plot the bounded preds in polar co-ords
  for jj=1:length(pout.obsData);
    npred=1;
    pred=gPred(pout.obsData(jj).x(1),poutHO(jj).pvals(pvec), ...
               poutHO(jj).model,poutHO(jj).data,'uvpred');
    predv=pred(:,1:pv);
    predu=pred(:,pv+1:end);
    
    pyd=(pout.obsData(jj).Dsim*predv')';
    pye=(pout.simData.Keta*predu')';
    py=pyd+pye;        
    
    pym=mean(py);
    pyh=prctile(py,95);
    pyl=prctile(py,5);
    
    ysd=pout.simData.yesd; ymean=pout.simData.yemean';
    vm=pym*ysd + ymean;
    vl=pyl*ysd + ymean;
    vh=pyh*ysd + ymean;
    
    vm=reshape(vm,size(pout.simData.timemat));
    vl=reshape(vl,size(pout.simData.timemat));
    vh=reshape(vh,size(pout.simData.timemat));
    
    d1u=unique(pout.obsData(jj).time);
    for ii=1:length(d1u)
        subplot(3,3,jj+3*((3-length(d1u))+ii-1))
          hold on;
        dinds=pout.obsData(jj).time==d1u(ii);
    
        [x,y]=pol2cart(pout.obsData(jj).phi(dinds),pout.obsData(jj).y(dinds));
          plot(x,y,'o','MarkerFaceColor','b','MarkerSize',3);
          axis([-2.5 2.5 -2.5 2.5]); axis equal;
          ylabel(['Time ' num2str(floor(d1u(ii)*1e6)) ' \mus']);
        timeline=pout.simData.timemat(:,1);
        tdist=timeline-d1u(ii);
        tdist1=tdist; tdist1(tdist1>0)=-Inf; [y1 jj1]=max(tdist1); t1=timeline(jj1);
        tdist2=tdist; tdist2(tdist2<0)= Inf; [y2 jj2]=min(tdist2); t2=timeline(jj2);
        c=(d1u(ii)-timeline(jj1))/(timeline(jj2)-timeline(jj1));
        vmp = (1-c)* vm(jj1,:) + c*vm(jj2,:); 
        vhp = (1-c)* vh(jj1,:) + c*vh(jj2,:); 
        vlp = (1-c)* vl(jj1,:) + c*vl(jj2,:); 
        vlp=max(vlp,0);
        %[x,y]=pol2cart(pout.simData.phimat(1,:),vmp); plot(x,y);     %mean
        [x,y]=pol2cart(pout.simData.phimat(1,:),vlp); plot(x,y,'k'); %lower
        [x,y]=pol2cart(pout.simData.phimat(1,:),vhp); plot(x,y,'k'); %upper
    end
    drawnow
  end
  axes('Position',[0 0 1 1],'Visible','off');
  text(0.17,0.96,'Experiment 1','FontSize',14)
  text(0.44,0.96,'Experiment 2','FontSize',14)
  text(0.73,0.96,'Experiment 3','FontSize',14)
end

if doPlot(7)
    figure(7);
    showPvals(pvals);
end

if doPlot(8) % misc plots for paper
  if newFig; figure(8); end
  clf;
  AzEl=[140 30];

  % first frame is a summary of the data and the sims with the full surface
  h1=subplot(2,2,1); cla
    for ii=1:36
      plot3(repmat(simData.x(ii,1),22,1),simData.orig.time,...
            simData.orig.y(1:22,ii))
      hold on
    end
    grid on;
    xlabel('x'); ylabel('time');zlabel('radius')  
    axis([0 1 0 5e-5 0 3])
    view(AzEl);
  
    plot3(pout1.simData.x(:,1), ...
        pout1.simData.x(:,2) * pout1.simData.orig.trange, ...
        pout1.simData.orig.y,'.')

    for ii=1:length(pout1.obsData)
      y=pout1.obsData(ii).orig.y;
      ymm(1)=min(y); ymm(2)=max(y);
      plot3([1 1]*pout1.obsData(ii).x(1), ...
            [1 1]*pout1.obsData(ii).x(2) * pout1.simData.orig.trange, ...
            ymm,'m','lineWidth',6)
      plot3([1 1]*pout1.obsData(ii).x(1), ...
            [1 1]*pout1.obsData(ii).x(2) * pout1.simData.orig.trange, ...
            [0 ymm(2)],'m-');
    end
  
  subplot(2,2,2); cla
   if(1)
    t=[pout1.pvals(ilinspace(1000,10000,1000)).theta]'; 
    plot(pout1.simData.x(:,3),pout1.simData.x(:,4),'mx')
    hold on
    xlabel('\theta_1'); ylabel('\theta_2')
    % plot the estimated prob contours of the calibration
      %plot(t(:,1),t(:,2),'.'); not the data
      ksd=0.05; ngrid=20; Pcontours=[0.1 0.5 0.9];
      % Generate a grid and supporting data structures
        gridvals = linspace(0,1,ngrid);
        [g1 g2] = meshgrid(gridvals,gridvals);
        g1v = g1(:); g2v = g2(:);
        gvlen = length(g1v);
        dens = zeros(gvlen,1);        
         for i=1:gvlen
           f = normpdf(t(:,1),g1v(i),ksd).*normpdf(t(:,2),g2v(i),ksd);
           dens(i) = sum(f);
         end
       % normalize dens
         dens = dens/sum(dens);
       % get the contours
         for j=1:length(Pcontours)
          hlevels(j) = fzero(@(x) getp(x,dens)-Pcontours(j),[0 max(dens)]);
         end
         %imagesc(g1v,g2v,reshape(dens,ngrid,ngrid)); axis xy; hold on;
         colormap(repmat([.9:-.02:.3]',[1 3]));
         contour(g1,g2,reshape(dens,ngrid,ngrid),hlevels,'LineWidth',1.0,'Color','b'); 
    end

  % predictions for further plots      
    %for ii=1:length(pout1.obsData); xpred(ii,:)=pout1.obsData(ii).x; end
    pgrid=linspace(0,1,10);
    glen=length(pgrid);  mlen=glen*glen;  
    [gridx,gridy]=meshgrid(pgrid,pgrid);
    pred=gPred([gridx(:) gridy(:)],pout1.pvals(pvec),pout1.model,pout1.data,'uvpred');

  subplot(2,2,3)
  % third frame is estimated discrepancy 
     ysd=pout1.simData.orig.ysd;
      du=squeeze(prctile(pred.v,95)) * ysd;
      dl=squeeze(prctile(pred.v,5)) * ysd;
      polyy=gridy(1:glen,1); polyy=[polyy; polyy(end:-1:1)] ;
      polyy=polyy * pout1.simData.orig.trange;
      for ii=1:glen
         polyx=repmat(gridy(ii,1),[glen*2 1]);
         duseg=du((ii-1)*glen+1:ii*glen); dlseg=dl((ii-1)*glen+1:ii*glen);
         polyz=[duseg; dlseg(end:-1:1)];
         fill3(polyx,polyy,polyz,'k'); hold on
       end

      %surf(gridx,gridy * pout1.simData.orig.trange, ...
      %     reshape(du,glen,glen),'facecolor','b','facealpha',0.5);
      %hold on
      %surf(gridx,gridy * pout1.simData.orig.trange, ...
      %     reshape(dl,glen,glen),'facecolor','b','facealpha',0.5);
      axis([0 1 0 5e-5 -1 1])
      grid on;
        xlabel('x'); ylabel('time');zlabel('discrepancy (cm)')  
       view(AzEl);

  subplot(2,2,4); cla
  % fourth frame is estimated eta 
     ymean=pout1.simData.orig.ymean;
     ysd=pout1.simData.orig.ysd;
      eu=squeeze(prctile(pred.u,95))*ysd + ymean;
      el=squeeze(prctile(pred.u,5))*ysd + ymean;
      polyy=gridy(1:glen,1); polyy=[polyy; polyy(end:-1:1)] ;
      polyy=polyy * pout1.simData.orig.trange;
      for ii=1:glen
         polyx=repmat(gridy(ii,1),[glen*2 1]);
         euseg=eu((ii-1)*glen+1:ii*glen); elseg=el((ii-1)*glen+1:ii*glen);
         polyz=[euseg; elseg(end:-1:1)];
         fill3(polyx,polyy,polyz,'k'); hold on
       end
      %surf(gridx,gridy * pout1.simData.orig.trange, ...
      %     reshape(eu,glen,glen),'facecolor','b','facealpha',0.5);
      %hold on
      %surf(gridx,gridy * pout1.simData.orig.trange, ...
      %     reshape(el,glen,glen),'facecolor','b','facealpha',0.5);
      axis([0 1 0 5e-5 0 3])
      grid on;
        xlabel('x'); ylabel('time');zlabel('calibrated pred. (cm)')  
       view(AzEl);

    for ii=1:length(pout1.obsData)
      y=pout1.obsData(ii).orig.y;
      ymm(1)=min(y); ymm(2)=max(y);
      plot3([1 1]*pout1.obsData(ii).x(1), ...
            [1 1]*pout1.obsData(ii).x(2) * pout1.simData.orig.trange, ...
            ymm,'m','lineWidth',6)
      plot3([1 1]*pout1.obsData(ii).x(1), ...
            [1 1]*pout1.obsData(ii).x(2) * pout1.simData.orig.trange, ...
            [0 ymm(2)],'m-');
    end

end

return

if ~exist('neddfigs'); mkdir('neddfigs'); end
cd neddfigs
figure(1); print -deps2c rhobox.eps
figure(2); print -deps2c thetascatter.eps
figure(3); print -deps2c thetapreds.eps
figure(4); print -deps2c datapreds.eps
figure(5); print -deps2c polarBounds.eps
figure(6); print -deps2c polarHoldoutBounds.eps
figure(7); print -deps2c pvals.eps
cd ..

return

end

% function to get probability of a given level h
function pout = getp(h,d);
    iabove = (d >= h);
    pout = sum(d(iabove));
end

