function generateLatinTowerData()

% notes for me:
% x = R
% theta = C
% y = {h, t}, i.e., pairs of h and t that form a trace when plotted

% imagine the field experiments involve say 4 platforms --> 4 values of h.
% Then for each R, one experiment gives output of 4 h-t pairs (a curve).
% Likewise for the simulator, we have a dense grid of say 100 heights h.
% Then for each setting of {x, theta} = {R, C} we get output of 100 h-t
% pairs.

% I'll make three .m files:
% 1. one to generate the h-t pairs and write them into files.
% 2. a "runmcmc"-type file that first calls...
% 3. ...a file that reads in the data and packages it appropriately


% generate "field" data and "simulator" data, where the simulator model is 
% systematically off from reality.

% true: d2h/dt2 = g - C (dh/dt)^2 / R
% sim:  d2h/dt2 = g - C (dh/dt) / R

% inputs for field experiments: x = R
% inputs for simulator: x = R, theta = C
% We want to calibrate theta in the simulator to match the field data.

usesaved = 0; % whether to run the optimization again or just use an earlier saved run

if usesaved
    load Run25design.mat;
else
    % values for generating data
    hfield = 5:5:20; % platform heights for the field experiments
    hsim = 1.5:1.5:25; % grid of heights fed to the simulator
    hdense = [0:.01:1.99, 2:.5:25]; % a denser grid for drawing the curves
    g = 9.8; % gravity
    et = 0; %0.01; % observation error on the experimental times
    Ctrue = 0.1 / (4 * pi / 3); 
    % the coefficient of drag for a smooth sphere is 0.1, and we're
    % dividing by 4/3 pi to absorb a constant related to the volume of the
    % sphere (not including R)
    R = [.1 .2 .4]; % radii of balls to try (in meters)

    % get a Latin hypercube design of 25 points over Rsim, Csim
    % using lhsdesign(25,2).  Here's one I like on [0,1]^2:
    des = [
        0.1239    0.8024;
        0.8738    0.6473;
        0.6140    0.3337;
        0.8833    0.4783;
        0.9946    0.0548;
        0.1178    0.9382;
        0.1805    0.2411;
        0.6638    0.2861;
        0.2939    0.1208;
        0.2451    0.2397;
        0.4577    0.5696;
        0.4377    0.8874;
        0.0737    0.7384;
        0.6931    0.8683;
        0.4901    0.7070;
        0.5953    0.9828;
        0.7506    0.1009;
        0.7783    0.4225;
        0.8333    0.5318;
        0.3987    0.6312;
        0.2021    0.4990;
        0.3495    0.3680;
        0.9411    0.7935;
        0.0198    0.0218;
        0.5440    0.1925];

    % scale the first column to [0,.5] and call it Rsim
    % (this includes our field values, i.e., R \in [0,.5])
    % scale the second column to [0.05,.25] and call it Csim
    % (likewise, Ctrue \in [0.05, .25])
    Rsim = des(:,1) .* .4 + .05;
    Csim = des(:,2) .* .2 + .05;
    
%     % shit: should I divide Csim by 4/3 pi as we did for Ctrue?
%     Csim = (des(:,2) .* .2 + .05) / (4 * pi / 3);
%     % NO

    % iterate over ball sizes to get field data
    tfield = zeros(length(R), length(hfield));
    tfdense = zeros(length(R), length(hdense));
    for rr = 1:length(R)
        % generate field data
        tfield(rr,:) = invertHtrue(hfield, g, Ctrue, R(rr), et); % observed times
        tfdense(rr,:) = invertHtrue(hdense, g, Ctrue, R(rr), et); % dense grid

        ymax(rr) = max(tfdense(rr,:)); % for axis limits
    end
    
    % imagine that the biggest ball is too big to get to the highest
    % platform, so we don't observe data there
    tfield(end,end) = NaN;

    % iterate over ball sizes to get simulator data
    tsim = zeros(length(Rsim), length(hsim));
    tsdense = zeros(length(Rsim), length(hdense));
    for rr = 1:length(Rsim)
        % generate simulator data
        tsim(rr,:) = invertHsim(hsim, g, Csim(rr), Rsim(rr))';
        tsdense(rr,:) = invertHsim(hdense, g, Csim(rr), Rsim(rr));

        ymax(rr + length(R)) = max(max(tsdense(rr,:))); % for axis limits
        %     ymax(rr + length(R)) = max(max(tsim{rr})); % for axis limits
    end

    ymax = max(ymax);

    save Run25design.mat;
end

% get the three design points closest to each R
sets = [ 7 1 6; 22 20 12; 4 19 2 ];
    
colors = ['r', 'g', 'b', 'm', 'c'];

% figure for the manual
figure;
for rr = 1:length(R)

    lgnd{1} = 'Reality';
    lgnd{2} = 'Field data';
    lgnd{3} = 'Simulation runs';

    % plot field data
    subplot(2,2,rr)
    leg(1) = plot(hdense, tfdense(rr,:), 'k-'); % leg() is legend line-type entry
    hold on;
    leg(2) = plot(hfield, tfield(rr,:), 'ks','MarkerSize',10, 'MarkerFaceColor','k');

    % add curves for all 25 design points in light gray
    leg(3) = plot(hdense, tsdense(1,:),'.-', 'color', [.8 .8 .8]); % plot just one to get handle for legend
    plot(hdense, tsdense,'.-', 'color', [.8 .8 .8]); % plot all 25
    
    legend(leg, lgnd, 'Location','SouthEast');

    % overlay curves for the design points with the closest Rsim
    for pp = 1:size(sets,1)
        plot(hdense,tsdense(sets(rr,pp),:),strcat('-', colors(pp)));
        plot(hsim, tsim(sets(rr,pp),:), strcat('.', colors(pp)));
    end
    
    % add the field data again so it's on top
    plot(hdense, tfdense(rr,:), 'k-');
    plot(hfield, tfield(rr,:), 'ks','MarkerSize',10, 'MarkerFaceColor','k');


    ylim([0 ymax]);
    xlabel('Height (m)')
    ylabel ('Time (s)')
    title(['Ball radius ', num2str(R(rr))]);
    
    % add an inset plot of the design, Rsim x Csim
    if rr == 1
        h = axes('position',[0.16 0.8 0.1 0.1],'color','w');
    elseif rr == 2
        h = axes('position',[0.6 0.8 0.1 0.1],'color','w');
    else
        h = axes('position',[0.16 0.325 0.1 0.1],'color','w');
    end
        
    plot(Rsim, Csim, 'o', 'color', [.8 .8 .8], 'parent', h)
    hold on;
    plot([R(rr), R(rr)], [min(Csim), max(Csim)],'k'); % add line for field R
    % add circles for the ones were plotting
    for pp = 1:size(sets,1)
        plot(Rsim(sets(rr,pp)), Csim(sets(rr,pp)), strcat('o',colors(pp)));
    end
    xlabel('R design values');
    ylabel('C design values');

end

% write the h-t pairs into files
dirstr = './';

% sim.dat, should be length(hsim) x length(Csim)
tsim = tsim';
save([dirstr 'sim.dat'], 'tsim', '-ASCII');

% sim.height, a file with just the heights (same for all sim runs)
hsim = hsim';
save([dirstr 'sim.height'], 'hsim', '-ASCII');

% sim.design, length(Csim) x (num X's + num thetas)
designout = [Rsim Csim];
save([dirstr 'sim.design'], 'designout', '-ASCII');

% field.dat, one row per experiment (radius)
save([dirstr 'field.dat'], 'tfield', '-ASCII');

% field.height
save([dirstr 'field.height'], 'hfield', '-ASCII');

% field radii
save([dirstr 'field.radii'], 'R', '-ASCII');

