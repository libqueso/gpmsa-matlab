function params=sc(doPlot);

 if ~exist('doPlot'); doPlot=0; else; doPlot=1; end

% read in design
 design=textread('design.txt'); m=size(design,1);
% scale design
 xmin=min(design); xrange=max(design)-xmin;
 design=(design-repmat(xmin,m,1))./repmat(xrange,m,1);

% read in sim data
 ysim=textread('sim_outputs'); ysim=ysim';
 
% plot sims
 if (doPlot)
   figure(1);
   for ii=1:size(design,2)
     xxs=design(:,ii)'; xxs=xxs*xrange(ii)+xmin(ii);
     subplot(2,2,ii);
     plot(xxs,ysim,'y*'); hold on;
     xlabel(['x_' int2str(ii)],'FontSize',10);
     ylabel('time','FontSize',10);
     xlim([xmin(ii)-0.1*xrange(ii) xmin(ii)+1.1*xrange(ii)]);
   end
   figure(1); print -depsc2 scDat; close;
 end

% standardize sims
 ysimmean=mean(ysim);
 ysimStd=ysim-ysimmean;
 ysimsd=std(ysimStd);
 ysimStd=ysimStd/ysimsd;

% K basis
% compute on sim grids
 Ksim=1;

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
 simData.orig.xmin=xmin;
 simData.orig.xrange=xrange;
% obsData
 obsData=[];

% pack up and leave
 params.simData=simData;
 params.obsData=obsData;

end
