function params=sc(doPlot);

 if ~exist('doPlot'); doPlot=0; else; doPlot=1; end

% read in design
 design=textread('design.txt'); m=size(design,1);
% scale design
 xmin=min(design); xrange=max(design)-xmin;
 design=(design-repmat(xmin,m,1))./repmat(xrange,m,1);

% read in sim data
 ysim=textread('sim_outputs'); ysim=ysim';
 
% read in obs data
 obsdata=textread('obs_outputs'); 
 n=size(obsdata,1); % number of experiments
% x in first nx columns, field data in last
 x=obsdata(:,1:end-1); nx=size(x,2);
 for ii=1:n, yobs(ii).y=obsdata(ii,end); end
% specify assumed observation errors here
 for ii=1:n, yobs(ii).Sigy=0.075^2; end
% scale obs x
 for ii=1:nx, x(:,ii)=(x(:,ii)-xmin(ii))./xrange(ii); end

% plot sims and data
 if (doPlot)
   xxs=design(:,1); xxs=xxs*xrange(1)+xmin(1);
   xxo=x; xxo=xxo*xrange(1)+xmin(1);
   yo=[yobs.y]'; sdo=sqrt([yobs.Sigy])';
   zc=norminv(.95);
   figure(1);
   plot(xxs,ysim,'y*'); hold on;
   xlabel('x','FontSize',14);
   ylabel('time','FontSize',14);
   plot(xxo,yo,'b.','MarkerSize',4);
   xlim([xmin(1)-0.1*xrange(1) xmin(1)+1.1*xrange(1)]);
   for ii=1:n
     hold on;
     plot([xxo(ii) xxo(ii)],yo(ii)+zc*[-sdo(ii) sdo(ii)],'k','LineWidth',2);
   end  
   figure(1); print -depsc2 scDat; close;
 end

% standardize sims
 ysimmean=mean(ysim);
 ysimStd=ysim-ysimmean;
 ysimsd=std(ysimStd);
 ysimStd=ysimStd/ysimsd;

% standardize experimental data
 for ii=1:n
   yobs(ii).ymean=ysimmean;
   yobs(ii).yStd=(yobs(ii).y-yobs(ii).ymean)/ysimsd;
   yobs(ii).sd=sqrt(yobs(ii).Sigy);
   yobs(ii).Sigy=yobs(ii).Sigy/(ysimsd.^2);
 end

% K basis
% compute on sim and data grids
 Ksim=1;
 for ii=1:n, yobs(ii).Kobs=1; end

% D basis
% compute on sim and data grids
 Dsim=1;
 for ii=1:n; yobs(ii).Dobs=1; end

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
 simData.orig.Dsim=Dsim;
 simData.orig.xmin=xmin;
 simData.orig.xrange=xrange;
% obsData
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
