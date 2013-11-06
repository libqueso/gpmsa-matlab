function [t] = invertHsim(h, g, C, R)
% I have the form of h(t), but what I really want is t(h) to get data for the tower problem.
% I'll use an optimizer to do this inversion.

% inputs:
% h, starting height of the thing we're dropping (can be a vector of
% starting heights)
% g, gravity 
% C, coefficient of drag (air resistance)
% R, radius of the ball we're dropping
% et, standard deviation on the observed time t given a height h

% g, C, and R can be vectors of values to try

% simulation:
% d^2h / dt^2 = g - (C / R) (dh / dt), which I claim works out to
% h(t) = exp{-(C/R) * t} / (C/R)^2 + g * t / (C/R)

% Need to solve for t given h, and do this for each setting of g, C, and R.

% By minimizing the following function wrt t, we'll find the right time for
% each given height h.
for hh = 1:length(h)
    for cc = 1:length(C)
        for rr = 1:length(R)
            banana = @(thh) (h(hh) - dragheight(thh, g, C(cc), R(rr))).^2;

            % Use a starting point based on no drag
            t0 = sqrt(2 * h(hh) / g);%(gg));

            % Find the time.
            thh = fminsearch(banana, t0);

            % Pardon the kluge...
            if length(C) > 1 && length(R) > 1
                t(hh,cc,rr) = thh;
            elseif length(R) > 1
                t(hh,rr) = thh;
            else
                t(hh,cc) = thh;
            end
        end

    end
end

%% given a time t, what's the height from which the ball was dropped
function h = dragheight(t, g, C, R)
%h = (g/(d^2))*(exp(-d*t)-1) + (g/d)*t;

b = C / R;

% K1 and K2 are constants of integration
K1 = -log(g)/b;
K2 = -g/b^2;

h = g*t/b + exp(-b*(t+K1))/b^2 + K2;
