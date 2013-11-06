function nedd1DPlots(pout,pvec, plotnum, varargin)

newFig=1;
for ii=1:2:length(varargin)
  switch varargin{ii}
  case 'newFig'
    newFig=varargin{ii+1};
  case 'holdOutMods'
    poutHO=varargin{ii+1};
  otherwise
    error('invalid extended argument passed to neddPlots');
  end
end

simData=pout.simData;
obsData=pout.obsData;

model=pout.model;
data=pout.data;
pvals=pout.pvals(pvec);

pu=model.pu; pv=model.pv;
p=model.p; q=model.q;

doPlot(10)=0;
if exist('plotnum','var'); doPlot(plotnum)=1;
else
end

% Set up prediction grid
grid=linspace(0,1,11); glen=length(grid);
[gridx,gridy]=meshgrid(grid,grid);

if doPlot(1) 
  figure(1); clf;
  xpred=[gridx(:) gridy(:)];
  pred=gPred(xpred,pvals,model,data,'uvpred'); 
  eta=squeeze(mean(pred.u));
  mesh(gridx,gridy,reshape(eta,size(gridx)),...
       'edgecolor',[0 0 0],'faceAlpha',0)
  xlabel('x'); ylabel('time'); zlabel('radius');
  
  hold on;
  for ii=1:size(simData.x,1)
    plot3(simData.x(ii,1),simData.x(ii,2),simData.yStd(ii),'.');
  end
  
  a=axis; bot=a(5)
  
  for ii=1:length(obsData);
    plot3(obsData(ii).x(1),obsData(ii).x(2),obsData(ii).yStd,'g*');
    plot3([1 1]*obsData(ii).x(1),[1 1]*obsData(ii).x(2),[bot obsData(ii).yStd],'g');
  end
  
  figure(2); clf;
  mesh(gridx,gridy,reshape(squeeze(mean(pred.v)),size(gridx)), ...
       'edgecolor',[0 0 0],'faceAlpha',0)
  xlabel('x'); ylabel('time'); zlabel('radius');

end

if doPlot(3)
  figure(3); clf
  
  xpred=reshape([pout.obsData.x],2,6)';
  pred=gPred(xpred,pvals,model,data,'uvpred');
  eta=squeeze(mean(pred.u + pred.v));
  etal=squeeze(prctile(pred.u + pred.v,5));
  etah=squeeze(prctile(pred.u + pred.v,95));
  
  obsix={[1 2 3],[4],[5 6]}; obsix2=[1 1 1 2 3 3];
  for ii=1:length(obsix)
    h(ii)=subplot(length(obsix),1,ii); cla; hold on;
    for jj=1:length(obsix{ii})
      ix=obsix{ii}(jj);
      plot(xpred(ix,2),eta(ix),'.');
      plot([1 1]*xpred(ix,2),[etal(ix) etah(ix)]);
      plot(obsData(ix).x(2),obsData(ix).yStd,'g.');
    end
  end
  axisNorm(h,'xy',[0 1 -2.5 1]);
  xlabel(h(3),'time');
  ylabel(h(1),'Exp. 1')
  ylabel(h(2),'Exp. 2')
  ylabel(h(3),'Exp. 3')

  % do some holdout predictions
  for ii=1:length(obsix2)
    predObsData=obsData(setxor(obsix{obsix2(ii)},1:6));
    predMod=setupModel(predObsData,simData);
    pred=gPred(pout.obsData(ii).x,pvals,predMod.model,predMod.data,'uvpred');
    etap(ii)=squeeze(mean(pred.u + pred.v));
    etapl(ii)=squeeze(prctile(pred.u + pred.v,5));
    etaph(ii)=squeeze(prctile(pred.u + pred.v,95));
  end
  
  
  obsStr={'x','o','^'}; 
  figure(4); clf; hold on
  for ii=1:length(obsData)
    plot(ii,obsData(ii).yStd-eta(ii),['g' obsStr{obsix2(ii)}]);
    plot([1 1]*ii,obsData(ii).yStd-[etal(ii) etah(ii)],'g');
    %plot(ii,0,['b' obsStr{obsix2(ii)}]);
    plot(ii+0.1,      obsData(ii).yStd-etap(ii),['m' obsStr{obsix2(ii)}]);
    plot([1 1]*ii+0.1,obsData(ii).yStd-[etapl(ii) etaph(ii)],'m');
  end
  plot([0.7 6.3],[0 0],'k:');
  
  axisNorm(gca,'x',[0.5 6.5]);
  set(gca,'XTick',[1:6]);
  xlabel('Observation number');
  ylabel('Residual observation error mean and 5/95% bounds');
  
end


%%%%%%%%%%%%%%
% Junk
%%%%%%%%%%%%%%
if doPlot(5)
    % Now the realizations over each theta
    if newFig; figure(3); end;
    clf; colormap('copper');
    npred=glen;
    AzEl=[45 55];    
    for ii=1:p+q
        tdat=ones(npred,p+q)*0.5;
        tdat(:,ii)=grid';
        xpred=tdat(:,1:p);
        theta=tdat(:,p+1:end);
        pred=gPred(xpred,pvals,model,data,'wpred',theta);
        pm=squeeze(mean(pred.w));
        
        wrm=reshape(pm,npred,pu);
        %r=(pout.simData.Keta*wrm')';
        %v=r*pout.simData.yesd + repmat(pout.simData.yemean(:)',npred,1);
        r=(pout.simData.Ksim*wrm')'*pout.simData.orig.ysd;
        v=r + repmat(pout.simData.orig.ymean(:)',npred,1);
        r=r(:,1:22);
        v=v(:,1:22);
        
        h1(ii)=subplot(p+q,3,ii*3-2);
        mesh(pout.simData.orig.time,grid,r,'MeshStyle','row'); view(AzEl);
        if ii<=p; var=['x_' num2str(ii)];else; var=['\theta_' num2str(ii-p)];end
        xlabel('time'); ylabel(var); zlabel('std. radius');
        alpha(0.25)
        axis tight

        h2(ii)=subplot(p+q,3,ii*3-1);
        mesh(pout.simData.orig.time,grid,v,'MeshStyle','row'); view(AzEl);
        xlabel('time'); ylabel(var); zlabel('radius');
        axis tight
        alpha(0.25)
        drawnow
    end
    axisNorm(h1,'xyzmax');
    axisNorm(h2,'xyzmax');
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
