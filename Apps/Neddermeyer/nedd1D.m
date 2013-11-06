function e=nedd1D(d)

simData=d.simData;
obsData=d.obsData;

% get the outer (sim) time range, for normalization (time is now an "x")
  simtime=simData.orig.time;
  stmin=simtime(1); strange=simtime(end) - simtime(1);

%OK, first, we need to get the times of the data mean at each point
ind=1;
for ii=1:length(obsData); 
  times=unique(obsData(ii).orig.time);
  for jj=1:length(times)
    obs(ind).x=[obsData(ii).x (times(jj)-stmin)/strange];
    obs(ind).y=obsData(ii).orig.y(obsData(ii).orig.time==times(jj));
    ind=ind+1;
  end
end

% do the some sim points
simtimeN=(simtime-stmin)/strange;
newtimes=linspace(0,1,6);
ind=1;
for ii=1:size(simData.orig.y,2)
  y=simData.orig.y(1:length(simtime),ii);
  newy=interp1(simtimeN,y,newtimes,'spline');
  for jj=1:length(newy)
    sim.x(ind,:)=[simData.x(ii,1) newtimes(jj)];
    sim.t(ind,:)=simData.x(ii,2:end);
    sim.y(:,ind)=newy(jj);
    ind=ind+1;
  end
end

% Organize into gpm struct
clear simData obsData

ymean=mean(sim.y); ysd=std(sim.y);

for ii=1:length(obs); 
  obsData(ii).x=obs(ii).x;
  obsData(ii).yStd=(mean(obs(ii).y) - ymean)/ysd;
  obsData(ii).Kobs=1;
  obsData(ii).Dobs=1;
  obsData(ii).orig.y=obs(ii).y;
end

simData.x=[sim.x sim.t];
simData.yStd=(sim.y-ymean)/ysd;
simData.Ksim=1;
simData.orig.ymean=ymean;
simData.orig.ysd=ysd;
simData.orig.y=sim.y;
simData.orig.tmin=stmin;
simData.orig.trange=strange;

e.obsData=obsData;
e.simData=simData;

