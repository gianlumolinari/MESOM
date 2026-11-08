function legs = patch_legs(Y, pars)
%PATCH_LEGS  Propagate the four legs of the patched trajectory.
%
%   Y = [tb1; tb2; tf1; tf2; X2]  (times in s, X2 in km and km/s)
%
%   leg 1: X1 at t = 0       forward  to t = tf1
%   leg 2: X2 at t = Tocc    backward to t = Tocc - tb1
%   leg 3: X2 at t = Tocc    forward  to t = Tocc + tf2
%   leg 4: X3 at t = T       backward to t = T - tb2
%
% Patch point A joins legs 1 and 2, patch point B joins legs 3 and 4.
% Each leg's states are returned dimensional [km, km/s], one row per step.

    persistent lastY lastLegs
    if ~isempty(lastY) && isequal(Y, lastY)
        legs = lastLegs;   % same Y as the previous call: reuse the result
        return
    end

    tb1 = Y(1); tb2 = Y(2); tf1 = Y(3); tf2 = Y(4); X2 = Y(5:10);
    TU = pars.TU;  sf = pars.nb.sfVec;
    opts = odeset('RelTol', 3e-12, 'AbsTol', 1e-14, ...
                  'Events', @(t, x) nbody_impacts(t, x, pars.nb));
    rhs = @(t, x) nbody_eom(t, x, pars.nb);

    [~, x] = ode113(rhs, [0, tf1] / TU, pars.X1 ./ sf, opts);
    legs{1} = x .* sf.';
    [~, x] = ode113(rhs, [pars.Tocc, pars.Tocc - tb1] / TU, X2 ./ sf, opts);
    legs{2} = x .* sf.';
    [~, x] = ode113(rhs, [pars.Tocc, pars.Tocc + tf2] / TU, X2 ./ sf, opts);
    legs{3} = x .* sf.';
    [~, x] = ode113(rhs, [pars.T, pars.T - tb2] / TU, pars.X3 ./ sf, opts);
    legs{4} = x .* sf.';

    lastY = Y;  lastLegs = legs;
end
