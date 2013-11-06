function [t] = invertHtrue(h, g, C, R, et)
% I have the form of h(t), but what I really want is t(h) to get data for the tower problem.
% I'll use an optimizer to do this inversion.

% inputs:
% h, starting height of the thing we're dropping (can be a vector of
% starting heights)
% g, gravity
% C, coefficient of drag (air resistance)
% R, radius of the ball we're dropping (could be a vector of radii to try)
% et, standard deviation on the observed time t given a height h

% truth:
% d^2h / dt^2 = g - (C / R) (dh / dt)^2, which I claim works out to
% h(t) = \ln(1 + exp{2t \sqrt{g C/R}) / (C/R) - t / \sqrt{(C/R) (1/g)}

% Need to solve for t given h.

% By minimizing the following function wrt t, we'll find the right time for
% each given height h and each radius R.

for ii = 1:length(h)
    for jj = 1:length(R)
        banana = @(tii) (h(ii) - dragheight(tii, g, C, R(jj))).^2;

        % Use a starting point based on no drag
        t0 = sqrt(2 * h(ii) / g);

        % Find the time.
        tii = fminsearch(banana, t0);

        % add some observation error?
        t(ii,jj) = tii + et * rand;
        
    end

end

%% given a time t, what's the height from which the ball was dropped
function h = dragheight(t, g, C, R)
%h = (g/(d^2))*(exp(-d*t)-1) + (g/d)*t;
CbyR = C / R;
%h = ( log(1 + exp(2 * t * sqrt(g * CbyR))) / CbyR ) - ( t / sqrt(CbyR / g) );
h = R/C * log( cosh( sqrt(CbyR * g) * t ) );