function data = readdata();

% we want to read in the inputs/outputs of the simuator and of the field
% experiments and get them into the format required by the GPM/SA code.

dirstr = './'; % where the files are

% read in the simulated data and the design
tsim = textread([dirstr 'sim.dat']); % times
hsim = textread([dirstr 'sim.height']); % heights
[Rsim Csim] = textread([dirstr 'sim.design']); % design (x=R and theta=C)

% read in the field data
tfield = textread([dirstr 'field.dat'])'; % times
hfield = textread([dirstr 'field.height']); % heights
Rfield = textread([dirstr 'field.radii']); % radii (x=R)
%
% NEED A SIGY - check with Dave whether to generate the field data with
% error or just to define a sigy here and pretend it's the right thing.
% Ok, Dave says make something up that matches the box size on the plots of the
% generated data, but generate the data without error.
% Sigy = .15; % this is box height...but does this need to be a matrix?
% Sigy = [.15 0; 0 .15]; % ??
%
m = size(tsim, 2); % number of simulation runs
n = size(tfield, 2); % number of experiments

% standardize the inputs to the simulator (x and theta) to lie in [0, 1]
Rsmin = min(Rsim);
Rsrange = range(Rsim);
Rsim01 = (Rsim - Rsmin) / Rsrange;

Csmin = min(Csim);
Csrange = range(Csim);
Csim01 = (Csim - Csmin) / Csrange;

% standardize the input to the field experiment (x) the same way
Rfield01 = (Rfield - Rsmin) / Rsrange;

% standardize the simulator output to have mean zero at each height and an
% overall variance of one
tsimmean = repmat(mean(tsim,2), [1 m]);
tsimStd = tsim - tsimmean; % makes mean at each height zero
tsimsd = std(tsimStd(:));
tsimStd = tsimStd / tsimsd; % makes overall variance one (but not at each height)

% standardize the field data the same way.  To do this we'll interpolate
% the simulation mean between the heights at which the simulations were run
% in order to find the interpolated simulation mean at heights
% at which the field experiments were done.  We'll use this interpolated
% mean to standardize the field data.

% since each experiment could have a different size (different number of 
% heights where the ball was dropped), we'll loop over
% experiments and put the results in an object called tobs
for ii = 1:n
    numhts = sum(~isnan(tfield(:, ii))); % how many heights have measurements for experiment ii
    tobs(ii).tfieldmean = interp1(hsim, tsimmean(:,1), hfield(1:numhts), 'linear', 'extrap');
    tobs(ii).tfieldStd = (tfield(1:numhts, ii) - tobs(ii).tfieldmean') / tsimsd;
    tobs(ii).hfield = hfield(1:numhts); % record the relevant heights while we're at it
    tobs(ii).tfield = tfield(1:numhts, ii); % carry along the unstandardized version
end

% -- K basis --
% want to capture the variation in the height-time curves across simulation
% runs; will do this by computing the K basis (using SVD in this case)
pu = 2; % number of basis components to keep
[U, S, V] = svd(tsimStd, 0);
Ksim = U(:, 1:pu) * S(1:pu, 1:pu) ./ sqrt(m); % the pu curves capturing variation across simulation runs

% now interpolate between height grids (like we did to standardize the
% field data above) to produce a corresponding Kobs, again storing the
% results for each experiment in tobs
for ii = 1:n
    tobs(ii).Kobs = zeros(length(tobs(ii).tfieldStd), pu);
    for jj = 1:pu % compute for each basis component
        tobs(ii).Kobs(:, jj) = interp1(hsim, Ksim(:, jj), tobs(ii).hfield, 'linear', 'extrap');
    end
end


% -- D basis --
% JGM: want to model the residuals of the observed data (as modeled by the Ksim
% basis)
% JG: lay it out, and record decomposition on sim and data grids
% JG: Kernel centers and widths
% Dgrid = 0:0.1:1; % locations on which the kernels are centered
% Dwidth = 0.1; % width of each kernel
Dgrid = 0:2:max(hsim); % locations on which the kernels are centered
Dwidth = 2; % width of each kernel
% RIGHT VALUES??
pv = length(Dgrid); % number of kernels

% JG: Compute the kernel function map, for each kernel
% Designate space for the Dsim matrix, 
% one row per simulated height, one column per kernel
% (consider making the grid of heights much denser for plotting)
Dsim = zeros(length(hsim), pv); 

% designate space for the Dobs matrix for each experiment, 
% one row per experimental height, one column per kernel
for ii = 1:n
    tobs(ii).Dobs = zeros(length(tobs(ii).tfieldStd), pv); 
end

% create each kernel
for jj = 1:pv
    % first create kernel jj for each experiment
    for ii = 1:n
        % normpdf computes the value of a Gaussian with mean
        % Dgrid(jj) and variance Dwidth at each element of hfield
        tobs(ii).Dobs(:, jj) = normpdf(tobs(ii).hfield, Dgrid(jj), Dwidth);
    end
    % now create kernel jj for the simulations
    Dsim(:, jj) = normpdf(hsim, Dgrid(jj), Dwidth);
end

% JG: normalize the basis elements of D
Dmax = max(max(Dsim * Dsim'));
Dsim = Dsim / sqrt(Dmax);
for ii = 1:n
    tobs(ii).Dobs = tobs(ii).Dobs / sqrt(Dmax);
end

% now package everything up into simData and obsData structures
% -- simData --
% required fields
simData.x = [Rsim01 Csim01]; % our design (standardized)
simData.yStd = tsimStd; % output, standardized
simData.Ksim = Ksim; 

% extra fields: original data and transform stuff
simData.orig.y = tsim;
simData.orig.ymean = tsimmean;
simData.orig.ysd = tsimsd;
simData.orig.Dsim = Dsim;
simData.orig.t = hsim; 
simData.orig.xorig = [Rsim Csim]; % original scale for simulated R, C

% -- obsData --
% JG: set the x values of density
for ii = 1:n
    % required fields
    obsData(ii).x = Rfield01(ii); 
    obsData(ii).yStd = tobs(ii).tfieldStd;
    obsData(ii).Kobs = tobs(ii).Kobs;
    obsData(ii).Dobs = tobs(ii).Dobs;
    %    obsData(ii).Sigy = tobs(ii).Sigy./(tsimsd.^2);
    
    % extra fields
    obsData(ii).orig.y = tobs(ii).tfield;
    obsData(ii).orig.ymean = tobs(ii).tfieldmean;
    obsData(ii).orig.t = tobs(ii).hfield;
end

% pack up and leave
data.simData = simData;
data.obsData = obsData;


