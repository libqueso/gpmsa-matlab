function dataStruct=neddeg(doPlot);

if ~exist('doPlot','var'); doPlot=0; else doPlot=1; end

new=0;
if new

% we'll be using the design below...
Xlhs = [    0.7714    0.4286    0.0286;
    0.3714    0.1143    0.7143;
    0.1714    0.4571    0.8857;
    0.3429    0.6000    0.8000;
    0.8000    0.6286    0.4000;
    0.7429    0.5429         0;
    0.6571    1.0000    0.6286;
    0.2857         0    0.4571;
    0.5143    0.9429    0.2286;
    0.6857    0.3143    0.6571;
    0.8286    0.2000    1.0000;
    0.9714    0.3429    0.6000;
    0.4000    0.8000    0.2000;
    0.5429    0.2857    0.2857;
    0.9143    0.8857    0.2571;
    0.0571    0.0286    0.0857;
    0.1143    0.5714    0.7429;
    0.2000    0.2286    0.3714;
    0.4571    0.9143    0.3429;
    0.6286    0.7143    0.6857;
         0    0.8286    0.9429;
    0.8857    0.0857    0.9714;
    0.2286    0.0571    0.5714;
    0.7143    0.1714    0.8571;
    0.2571    0.4857    0.1429;
    0.5714    0.4000    0.8286;
    0.9429    0.6857    0.4857;
    0.4857    0.1429    0.1143;
    1.0000    0.8571    0.9143;
    0.6000    0.6571    0.5143;
    0.1429    0.7429    0.5429;
    0.8571    0.2571    0.0571;
    0.3143    0.3714    0.4286;
    0.4286    0.7714    0.7714;
    0.0286    0.9714    0.3143;
    0.0857    0.5143    0.1714];

%%%%%%% Neddermeyer stuff below %%%%%%%%%
% these parameter values simulate expt 10 in Neddermeyer '43
% params = [R1 lam s rho mratio u0]
params10 = [1.5*2.54 2/3 3e5*6.84e4 7.5 .32 1.65e10]';
paramslhs = zeros([36 6]);
for i=1:length(Xlhs(:,1))
    paramslhs(i,:) = params10'.*[1 1 Xlhs(i,2)*.2+.9 1 ... 
        Xlhs(i,1)*.65+.5   Xlhs(i,3)*.2+.9];
end
% the computer runs will vary mratio from .32*[.5 to 1.15]
%                             s      from s0*[.9 to 1.1]
%                             u0     from u0*[.9 to 1.1]

tend = 5.0e-5; nt = 21;
time = (0:tend/nt:tend)';
rinner = nedderimp(time,params10);
lam = params10(2); R1 = params10(1);
router = sqrt(rinner.^2 + 1 - lam^2);
phi = (0:.04:1) .* 2*pi;
xcirc = cos(phi);
ycirc = sin(phi);
r = nedderimp(time,params10);
%%%% a plot of the implosion over time
% plot of cylinder for 10 frames
if doPlot
figure(1)
tframes=round(1:((nt-1)/10):nt);
for i=1:10
    it = tframes(i);
    D2subplot(i,2,5,[.5 .08 .05 .02],[0 0 0 0])
%    plot(xcirc*rinner(it)*R1,ycirc*rinner(it)*R1,'r',...
%         xcirc*router(it)*R1,ycirc*router(it)*R1,'b');
    patch(xcirc*router(it)*R1,ycirc*router(it)*R1,[.6 .6 .6]);
    patch(xcirc*rinner(it)*R1,ycirc*rinner(it)*R1,[1 1 1]);

         set(gca,'Xtick',[]); set(gca,'Ytick',[]); ylabel('');
         if(i==1) || (i==6)
             set(gca,'Ytick',[-2 0  2 ]);
             ylabel('cm');
         end
         axis equal; axis([-1 1 -1 1]*R1*1.2); 
         xlabel(''); hold on;
    %plotdome([2.2 1.9 1.6 1.3 ]',6,rinner(it)*R1); 
    plot(time(1:it)/(6e-5)*6-4,repmat(-4.2,[it 1]),'g','LineWidth',4);
    text(-4:1:1,repmat(-3.9,[1 6]),{'0','1','2','3','4','X 10^{-5} s'},'FontSize',7)
    hold off;
end
end

% generate simulation data according to the design Xlhs and paramslhs
m = size(Xlhs,1); 
for i=1:m
    params = paramslhs(i,:)';
    yr(:,i) = params(1).*nedderimp(time,params);
end
if doPlot
    figure(2)
    plot(time,yr,'g');
    axis tight; xlabel('time (s)'); ylabel('inner radius (cm)');
end

% define the time,phi values at which the sim will output...
nt = length(time); np = length(phi);
ye = repmat(yr,[np 1]);
yearr = reshape(ye,[nt np m]);
if doPlot
    figure(3)
    mesh(time,phi,yearr(:,:,2)'); hold on;
    mesh(time,phi,yearr(:,:,1)');
    mesh(time,phi,yearr(:,:,29)');
    set(gca,'YTick',[0 pi 2*pi],'YTickLabel',{'0','pi','2pi'});
    axis([0 max(time) 0 2*pi 0 2.6]);
    xlabel('time (s)'); ylabel('angle (radians)');zlabel('inner radius (cm)');
end

% generate discrepancy for experimental data
phidat = 0:(2*pi/16):(2*pi-.1); 
ndat1 = length(phidat);
xcircd = cos(phidat); ycircd = sin(phidat);
phiknots = 0:(2*pi/8):(2*pi-.1); pphiknots = length(phiknots);
Ddelt = zeros([ndat1 pphiknots]);
Dsim  = zeros([length(phi) pphiknots]);
for k=1:pphiknots
    Ddelt(:,k) = dnorm(dist2pi(phidat,phiknots(k)*ones(size(phidat))),0,pi/8)';
    Dsim(:,k) = dnorm(dist2pi(phi,phiknots(k)*ones(size(phi))),0,pi/8)';
end
dknots = [.04 -.03 .03 -.03 .02 -.03 .03 -.03]'*2.5;
datadelt = Ddelt*dknots;
simdelt = Dsim*dknots;

% fakesol = yearr(:,:,1) + repmat(0:1/21:1,[26 1])' .* repmat(simdelt,[1 22])';
%     mesh(time,phi,fakesol');
%     set(gca,'YTick',[0 pi 2*pi],'YTickLabel',{'0','pi','2pi'});
%     axis([0 max(time) 0 2*pi 0 2.6]);
%     xlabel('time (s)'); ylabel('angle (radians)');zlabel('inner radius (cm)');

% generate experimental data
% experiment 1 (expt 10 from Neddermeyer '43)
paramsdat1 = params10;
timedat1 = [1.5e-5 2.7e-5 4.5e-5]';
rdat1 = paramsdat1(1)*nedderimp(timedat1,paramsdat1);
rout1 = sqrt(rdat1.^2./(paramsdat1(1)^2) + 1 - lam^2).*paramsdat1(1);
ydat1 = repmat(rdat1',[ndat1 1]) + repmat(datadelt,[1 3]);
ydat1 = ydat1 + randn(size(ydat1)).*.01;
rsimout1 = repmat(rout1',[length(phi) 1]);
rsim1 = repmat(rdat1',[length(phi) 1]) + repmat(simdelt,[1 3]);
% experiment 2 (expt 2 from Neddermeyer '43)
paramsdat2 = params10; paramsdat2(5) = .17;
timedat2 = 4.5e-5;
rdat2 = paramsdat2(1)*nedderimp(timedat2,paramsdat2);
rout2 = sqrt(rdat2.^2./(paramsdat2(1)^2) + 1 - lam^2).*paramsdat2(1);
ydat2 = repmat(rdat2',[ndat1 1]) + repmat(datadelt,[1 1]);
ydat2 = ydat2 + randn(size(ydat2)).*.01;
rsimout2 = repmat(rout2',[length(phi) 1]);
rsim2 = repmat(rdat2',[length(phi) 1]) + repmat(simdelt,[1 1]);
% experiment 3 (expt 11 from Neddermeyer '43)
paramsdat3 = params10; paramsdat3(5) = .36;
timedat3 = [2.5e-5 4.5e-5]';
rdat3 = paramsdat3(1)*nedderimp(timedat3,paramsdat3);
rout3 = sqrt(rdat3.^2./(paramsdat3(1)^2) + 1 - lam^2).*paramsdat3(1);
ydat3 = repmat(rdat3',[ndat1 1]) + repmat(datadelt,[1 2]);
ydat3 = ydat3 + randn(size(ydat3)).*.01;
rsimout3 = repmat(rout3',[length(phi) 1]);
rsim3 = repmat(rdat3',[length(phi) 1]) + repmat(simdelt,[1 2]);
xdat = ([params10(5) .17 .36]/.32-.5)./.65;

if doPlot
figure(4)
D2subplot(1,3,5,[.35 .08 .05 .02],[0 0 0 0])
patch(xcirc'.*rsimout1(:,1),ycirc'.*rsimout1(:,1),[.6 .6 .6]); hold on;
patch(xcirc'.*rsim1(:,1),ycirc'.*rsim1(:,1),[1 1 1]);
         plot(xcircd'.*ydat1(:,1),ycircd'.*ydat1(:,1),'b.','MarkerSize',11);
         axis equal; axis([-1 1 -1 1]*3.8); title('experiment 1'); ylabel('cm');
         set(gca,'Xtick',[]); %set(gca,'Ytick',[]);

D2subplot(6,3,5,[.35 .08 .05 .02],[0 0 0 0])
patch(xcirc'.*rsimout1(:,2),ycirc'.*rsimout1(:,2),[.6 .6 .6]); hold on;
patch(xcirc'.*rsim1(:,2),ycirc'.*rsim1(:,2),[1 1 1]);
         plot(xcircd'.*ydat1(:,2),ycircd'.*ydat1(:,2),'b.','MarkerSize',11);
         axis equal; axis([-1 1 -1 1]*3.8); ylabel('cm');
         set(gca,'Xtick',[]); %set(gca,'Ytick',[]);
      
D2subplot(11,3,5,[.35 .08 .05 .02],[0 0 0 0])
patch(xcirc'.*rsimout1(:,3),ycirc'.*rsimout1(:,3),[.6 .6 .6]); hold on;
patch(xcirc'.*rsim1(:,3),ycirc'.*rsim1(:,3),[1 1 1]);
         plot(xcircd'.*ydat1(:,3),ycircd'.*ydat1(:,3),'b.','MarkerSize',11);
         axis equal; axis([-1 1 -1 1]*3.8); ylabel('cm');
         xlabel('cm');
         %set(gca,'Xtick',[]); %set(gca,'Ytick',[]);
 
D2subplot(12,3,5,[.35 .08 .05 .02],[0 0 0 0])
patch(xcirc'.*rsimout2(:,1),ycirc'.*rsimout2(:,1),[.6 .6 .6]); hold on;
patch(xcirc'.*rsim2(:,1),ycirc'.*rsim2(:,1),[1 1 1]);
         plot(xcircd'.*ydat2(:,1),ycircd'.*ydat2(:,1),'b.','MarkerSize',11);
         axis equal; axis([-1 1 -1 1]*3.8);  title('experiment 2');
         set(gca,'Ytick',[]); %set(gca,'Xtick',[]);
         xlabel('cm');

D2subplot(8,3,5,[.35 .08 .05 .02],[0 0 0 0])
patch(xcirc'.*rsimout3(:,1),ycirc'.*rsimout3(:,1),[.6 .6 .6]); hold on;
patch(xcirc'.*rsim3(:,1),ycirc'.*rsim3(:,1),[1 1 1]);
         plot(xcircd'.*ydat3(:,1),ycircd'.*ydat3(:,1),'b.','MarkerSize',11);
         axis equal; axis([-1 1 -1 1]*3.8);  title('experiment 3');
         set(gca,'Xtick',[]); set(gca,'Ytick',[]); 
      
D2subplot(13,3,5,[.35 .08 .05 .02],[0 0 0 0])
patch(xcirc'.*rsimout3(:,2),ycirc'.*rsimout3(:,2),[.6 .6 .6]); hold on;
patch(xcirc'.*rsim3(:,2),ycirc'.*rsim3(:,2),[1 1 1]);
         plot(xcircd'.*ydat3(:,2),ycircd'.*ydat3(:,2),'b.','MarkerSize',8);
         axis equal; axis([-1 1 -1 1]*3.8); 
         set(gca,'Ytick',[]); %set(gca,'Ytick',[]);
         xlabel('cm');
axes('Position',[0 0 1 1],'Visible','off');
text(.65,.85,'t = 15 \mus','FontSize',14);
text(.65,.65,'t = 25 \mus','FontSize',14);
text(.65,.45,'t = 45 \mus','FontSize',14);
end

save neddTemp
return

else
load neddTemp
end

doPlot=1;

% now the data come in 3 expts with ydat1 ydat2 ydat3 holding the
% data and timedat1 timedat2 timedat3 holding the times at which
% the r's were measured, and phidat holds the angles for each
% time slice of each dataset

%
%
% Data is now complete, build up the structures for GPM modeling.
%
%

% standardize the simulations...
ysimmean = mean(ye,2);
ysimsd = std(ye(:));
ysimStd = (ye - repmat(ysimmean,[1 m]))./ysimsd;

% now represent ye images using eof's
[U,S,V]=svd(ysimStd,0);

% we'll record these.
timearr = repmat(time,[1 length(phi)]);
phiarr = repmat(phi,[length(time) 1]);

% plot the 3-d eof representation
if doPlot
  % 3-d representation
  y3 = (U(:,1:3)*S(1:3,1:3)*V(:,1:3)')*ysimsd+repmat(ysimmean,[1 m]);
  y3e = reshape(y3,size(yearr));

  timearr = repmat(time,[1 length(phi)]);
  phiarr = repmat(phi,[length(time) 1]);

  figure(5);
  plot3(timearr,phiarr,yearr(:,:,k),'w');
  set(gca,'YTick',[0 pi 2*pi],'YTickLabel',{'0','pi','2pi'});
  xlabel('time'); ylabel('angle'); zlabel('radius');
  hold on;
  for k=1:m
    plot3(timearr,phiarr,yearr(:,:,k),'g');
    plot3(timearr,phiarr,y3e(:,:,k),'.c');
    %plot3(timearr,phiarr,y3e(:,:,k),'.g');
  end
  grid on; axis tight;
end

% fill up the simulator object - simData - with stuff.
simData.orig.y = ye;      % matrix of neta x m simulation output; each column is a different sim run.
simData.x = Xlhs;  % m x p vector of input settings (x,theta)
simData.orig.ymean = ysimmean; % neta-vector that's the mean of the m sim outputs.
simData.orig.ysd = ysimsd;  % sd of sims
simData.yStd=ysimStd;

numPC = 3; % number of eof bases used to represent the sim output vectors
% compute the matrix Keta which is neta x peta, each column corresponding
% to an eof
simData.Ksim = U(:,1:numPC)*S(1:numPC,1:numPC)./sqrt(m);
simData.orig.time = time; % the values of the time for the sim output in the image matrix format
simData.orig.phi = phi; % the values of the angle for the sim output in the image matrix format
simData.orig.timemat = timearr;
simData.orig.phimat = phiarr;

% now a plot to make sure the eof representation is approximating the the
% actual dataset.  The plot is over radius and phi here.
% 2-d and 1-d plots of the pc's
if doPlot
  figure(6)
  for k=1:numPC
   D2subplot(k,2,3,[.2 .02 0 0], [.08 .08 .02 .02]);
   mesh(phi,time,reshape(simData.Ksim(:,k),[length(time) length(phi)]))
   set(gca,'YTick',[0 pi 2*pi],'YTickLabel',{'0','pi','2pi'});
   ylabel('time'); xlabel('angle'); zlabel('r');
   axis([0 2*pi 0 max(time) -1.5 .4]);
  end
  for k=1:numPC
   D2subplot(k+3,2,3,[.2 .02 0 0], [.08 .08 .02 .02]);
   plot(time,simData.Ksim(1:nt,k))
   xlabel('time'); ylabel('r');
   axis([0 max(time) -1.5 .4]);
  end
end

% fill up the obsData() list - each component corresponds to data from a
% particular experiment. n = # of experiments (here n=3).  Everything  must
% be input here.  radius and angle are used later to interpolate the eof.
n = 3; % # of experiments
obsData(1).orig.y = ydat1(:);   % observed radii at various time slices
obsData(2).orig.y = ydat2(:);
obsData(3).orig.y = ydat3(:);
for k=1:n obsData(k).x = [xdat(k)]; end    % observed x,theta conditions (note theta values are irrelevant here - they're just place holders)
obsData(1).orig.time = [ones(size(phidat'))*timedat1(1); ones(size(phidat'))*timedat1(2);...
                   ones(size(phidat'))*timedat1(3)];
obsData(2).orig.time = [ones(size(phidat'))*timedat2(1)];
obsData(3).orig.time = [ones(size(phidat'))*timedat3(1); ones(size(phidat'))*timedat3(2)];
obsData(1).orig.phi =  repmat(phidat',[3 1]);  % phi angle values for the measurements.
obsData(2).orig.phi = repmat(phidat',[1 1]);
obsData(3).orig.phi = repmat(phidat',[2 1]);
% compute simulator mean values simdat.ymean interpolated to the data values...
for k=1:n 
    ymk = interp2(simData.orig.phimat,simData.orig.timemat, reshape(simData.orig.ymean,[length(simData.orig.time) length(simData.orig.phi)]),...
        unique(obsData(k).orig.phi'),unique(obsData(k).orig.time));
    ymk = ymk';
    obsData(k).orig.ymean = ymk(:);
end

% a plot to check things...
if doPlot
  figure(7)
  mesh(simData.orig.phimat, simData.orig.timemat, reshape(simData.orig.ymean,[22 26])); hold on;
  for k=1:n
    plot3(obsData(k).orig.phi,obsData(k).orig.time,obsData(k).orig.ymean,'o');
    plot3(obsData(k).orig.phi,obsData(k).orig.time,obsData(k).orig.y,'.');
  end
   set(gca,'XTick',[0 pi 2*pi],'XTickLabel',{'0','pi','2pi'});
   ylabel('time'); xlabel('angle'); zlabel('r');
   axis([0 2*pi 0 max(time) 0 2.8]);

  hold off;
end

% now compute the centered, scaled observed arrival times yStd
for k=1:n 
    obsData(k).yStd = (obsData(k).orig.y - obsData(k).orig.ymean)./simData.orig.ysd; 
end  

% now compute the interpolated eof's for each dataset - held in the 
% obsData(k).n x peta matrix Kdat
for k=1:n
    obsData(k).Kobs = zeros([length(obsData(k).yStd) numPC]);
    for j=1:numPC
        obsData(k).Kobs(:,j)= interp2(simData.orig.phimat, simData.orig.timemat, reshape(simData.Ksim(:,j),size(simData.orig.phimat)),...
            obsData(k).orig.phi,obsData(k).orig.time);
    end
end

% a plot to check things...
if doPlot
  figure(8);
  for j=1:numPC
    D2subplot(j,1,3,[.5 0 0 0],[.1 .1 .02 .02]);
    mesh(simData.orig.phimat, simData.orig.timemat, reshape(simData.Ksim(:,j),[22 26])); hold on;
    set(gca,'XTick',[0 pi 2*pi],'XTickLabel',{'0','pi','2pi'});
    ylabel('time'); xlabel('angle'); zlabel('r'); axis tight;
    %axis([0 2*pi 0 max(time) 0 2.8]);
    for k=1:n
        plot3(obsData(k).orig.phi,obsData(k).orig.time,obsData(k).Kobs(:,j),'k.','MarkerSize',14);
    end
  end
  hold off;
end

% compute the basis functions for the discrepancy function.  Each data set
% will get a obsData(k).n x pdelta matrix of basis functions.  
phiknots = 0:(2*pi/8):(2*pi-.1); pphiknots = length(phiknots);
Ddelt = zeros([ndat1 pphiknots]);
Dsim  = zeros([length(phi) pphiknots]);
for k=1:pphiknots
    Ddelt(:,k) = dnorm(dist2pi(phidat,phiknots(k)*ones(size(phidat))),0,pi/8)';
    Dsim(:,k) = dnorm(dist2pi(phi,phiknots(k)*ones(size(phi))),0,pi/8)';
end
dknots = [.04 -.03 .03 -.03 .02 -.03 .03 -.03]'*2.5;
datadelt = Ddelt*dknots;
simdelt = Dsim*dknots;

phiknots = 0:(2*pi/8):(2*pi-.1); pphiknots = length(phiknots);
timeknots = [0:.25:.5]*1e-4; ptimeknots=length(timeknots);
for k=1:n
    obsData(k).orig.knotlocstime = reshape(repmat(timeknots,[pphiknots 1]),[ptimeknots*pphiknots 1]);
    obsData(k).orig.knotlocsphi = repmat(phiknots',[ptimeknots 1]);
    pv = length(obsData(k).orig.knotlocstime);
    obsData(k).Dobs = zeros([length(obsData(k).yStd) pv]);
    obsData(k).orig.Dsim = zeros([size(simData.yStd,1) pv]);
    for j=1:pv
        obsData(k).Dobs(:,j) = dnorm(obsData(k).orig.time,obsData(k).orig.knotlocstime(j),.25*1e-4).*...
            dnorm(dist2pi(obsData(k).orig.phi,obsData(k).orig.knotlocsphi(j)*ones(size(obsData(k).orig.phi))),0,pi/8);
        obsData(k).orig.Dsim(:,j) = dnorm(simData.orig.timemat(:),obsData(k).orig.knotlocstime(j),.25*1e-4).*...
            dnorm(dist2pi(simData.orig.phimat(:),obsData(k).orig.knotlocsphi(j)*ones(size(simData.orig.phimat(:)))),0,pi/8);
    end
end
%plot3(simData.timemat',simData.phimat',[reshape(obsData(1).Dsim(:,1),[22 26])]')
% now normalize Ddat and Dsim so that it gives a var=1 process
Dsim = obsData(1).orig.Dsim;
dmax = max(max(Dsim*Dsim'));
for k=1:n
    obsData(k).Dobs = obsData(k).Dobs/sqrt(dmax);
    obsData(k).orig.Dsim = obsData(k).orig.Dsim/sqrt(dmax);
end

% now some images to check things out...
if doPlot
  figure(5);
  imagesc(reshape(obsData(3).orig.Dsim,[22 26*24])); % the basis functions
  for k=1:n
    vtest = randn([size(obsData(k).Dobs,2) 1]);
    D2subplot(k,2,2,[0 0 0 0],[.1 .1 .02 .02]);
    mesh(simData.orig.timemat,simData.orig.phimat,reshape(obsData(k).orig.Dsim*vtest,size(simData.orig.timemat)));
    hold on;
    plot3(obsData(k).orig.time,obsData(k).orig.phi,obsData(k).Dobs*vtest,'.r');
  end
  hold off; % this checks out well
end
% plot the discrepancy bases
if doPlot
    figure(9)
    for k=1:size(obsData(1).Dobs,2)
        D2subplot(k,3,8,[.5 .1 0 0],[.01 .01 .01 .01]);
        imagesc(reshape(obsData(3).orig.Dsim(:,k),[22 26]));
        axis equal; axis tight;
        set(gca,'Ytick',[]); set(gca,'Xtick',[]);
        if k==9 ylabel('time','FontSize',13); end
    end
    axes('Position',[0 0 1 1],'Visible','off');
    text(.5,.45,'angle \phi','FontSize',13);
end
        
% now a figure to check - 1st, the fit to the simulator eofs
if doPlot

  % compute the uystar's for each experiment
  for k=1:n 
    DKtmp = [obsData(k).Dobs obsData(k).Kobs];
    pv = size(obsData(k).Dobs,2);
    Pobs = diag(ones([pv+numPC 1]))*0.00001;
    DKtmp2=DKtmp'*DKtmp + Pobs;
    DKtmp2inv=inv(DKtmp2);% + eye(size(DKtmp2))*1e-8);
    vu = DKtmp2inv*DKtmp'*obsData(k).yStd;
    %% sqrt(mean((obsData(k).y0-DKtmp*vu).^2))
    sdCC = sqrt(diag(DKtmp2inv));
    VC = sdCC*sdCC';
    %% CC./VC
    vudat(k).vstar = vu(1:pv);
    vudat(k).ustar = vu(pv+1:end);
  end

  figure(9);
  for k=1:n
    D2subplot(k,2,2,[0 0 0 0],[.1 .1 .02 .02]);
    mesh(simData.orig.phimat, simData.orig.timemat, reshape(simData.Ksim*vudat(k).ustar,[22 26])); hold on;
    plot3(obsData(k).orig.phi,obsData(k).orig.time,obsData(k).yStd,'o'); 
    plot3(obsData(k).orig.phi,obsData(k).orig.time,obsData(k).Kobs*vudat(k).ustar,'.r'); hidden off;
  end
%keyboard
  % 2nd, the fit to the simulator + delta
  figure(10);
  for k=1:n
    D2subplot(k,2,2,[0 0 0 0],[.1 .1 .02 .02]);
    mesh(simData.orig.phimat, simData.orig.timemat, reshape([obsData(k).orig.Dsim simData.Ksim]*[vudat(k).vstar; vudat(k).ustar],[22 26])); hold on;
    plot3(obsData(k).orig.phi,obsData(k).orig.time,obsData(k).yStd,'o'); 
    plot3(obsData(k).orig.phi,obsData(k).orig.time,[obsData(k).Dobs obsData(k).Kobs]*[vudat(k).vstar; vudat(k).ustar],'.r'); hidden off;
  end
end

dataStruct.obsData=obsData;
dataStruct.simData=simData;
