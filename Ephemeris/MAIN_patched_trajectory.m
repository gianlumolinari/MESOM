% MAIN_PATCHED_TRAJECTORY
% Moon-resonant CR3BP orbit transferred to a Sun-Earth-Moon point-mass
% ephemeris model (nbody_eom.m, SPICE ephemerides).
%
% In the ephemeris model the orbit is no longer periodic, so it is propagated
% in four legs from three CR3BP states (patch_legs.m):
%   leg 1: X1 (t = 0)    forward  to T_occ/2
%   leg 2: X2 (t = Tocc) backward to T_occ/2
%   leg 3: X2 (t = Tocc) forward  to (T + T_occ)/2
%   leg 4: X3 (t = T)    backward to (T + T_occ)/2
% The position gaps and velocity jumps at the two patch points measure how far
% the CR3BP solution is from a trajectory of the ephemeris model.
%
% Kernels: see load_kernels.m. Figures are written to ../figures/.

clear; close all;
here = fileparts(mfilename('fullpath'));
addpath(here);
load_kernels();

%% Units, epoch and CR3BP states (ECLIPJ2000, Earth-centred)
LU = 149597870.7;                  % [km]
TU = 365.25*3600*24/(2*pi);        % [s]
et0  = 1109232158.932959;          % initial epoch [s past J2000]
T    = 850464;                     % orbit period (1/3 of the synodic month) [s]
Tocc = 186399.7677438860;          % time of the occultation midpoint [s]

X1 = [198527.6679930048; 0.0000053992; 0; 0; -1.4299072557; 0];                   % t = 0
X2 = [45415.4279876202; -193534.7724471950; 0; -1.3885757048; -0.3266052740; 0];   % t = Tocc
X3 = [195688.3320424732; 33456.1154580864; 0; 0.2409696478; -1.4094567708; 0];     % t = T

pars.nb = nbody_setup(et0, '399', 'ECLIPJ2000', {'10', '399', '301'}, LU, TU);
pars.TU = TU;  pars.T = T;  pars.Tocc = Tocc;
pars.X1 = X1;  pars.X3 = X3;

%% Patched trajectory
Y = [Tocc/2; (T - Tocc)/2; Tocc/2; (T - Tocc)/2; X2];    % [tb1; tb2; tf1; tf2; X2]
legs = patch_legs(Y, pars);

gapA = norm(legs{1}(end,1:3) - legs{2}(end,1:3));   dvA = 1e3*norm(legs{1}(end,4:6) - legs{2}(end,4:6));
gapB = norm(legs{3}(end,1:3) - legs{4}(end,1:3));   dvB = 1e3*norm(legs{3}(end,4:6) - legs{4}(end,4:6));
fprintf('Patch point A: position gap %8.1f km, velocity jump %6.1f m/s\n', gapA, dvA);
fprintf('Patch point B: position gap %8.1f km, velocity jump %6.1f m/s\n', gapB, dvB);

moon    = cspice_spkpos('301', et0 + linspace(0, T, 400), 'ECLIPJ2000', 'NONE', '399');
moonOcc = cspice_spkpos('301', et0 + Tocc, 'ECLIPJ2000', 'NONE', '399');
cols = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19; 0.49 0.18 0.56];

fig1 = figure('Color', 'w', 'Position', [100 100 950 760]);
ax = axes(fig1); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on'); axis(ax, 'equal');
hm = plot(ax, moon(1,:)/1e3, moon(2,:)/1e3, ':', 'Color', [0.5 0.5 0.5]);
h = gobjects(4, 1);
for k = 1:4
    h(k) = plot(ax, legs{k}(:,1)/1e3, legs{k}(:,2)/1e3, 'Color', cols(k,:), 'LineWidth', 1.8);
end
ends = [legs{1}(end,1:2); legs{2}(end,1:2); legs{3}(end,1:2); legs{4}(end,1:2)] / 1e3;
hp = plot(ax, ends(:,1), ends(:,2), 'kd', 'MarkerFaceColor', 'y', 'MarkerSize', 8);
hx = plot(ax, X2(1)/1e3, X2(2)/1e3, 'kp', 'MarkerFaceColor', 'r', 'MarkerSize', 13);
text(ax, mean(ends(1:2,1)) + 15, mean(ends(1:2,2)), sprintf('A: %.0f km, %.0f m/s', gapA, dvA));
text(ax, mean(ends(3:4,1)) + 15, mean(ends(3:4,2)), sprintf('B: %.0f km, %.0f m/s', gapB, dvB));
hmo = plot(ax, moonOcc(1)/1e3, moonOcc(2)/1e3, 'o', 'MarkerFaceColor', [0.6 0.6 0.6], 'MarkerEdgeColor', 'k');
plot(ax, 0, 0, 'o', 'MarkerSize', 9, 'MarkerFaceColor', [0.2 0.4 0.9], 'MarkerEdgeColor', 'k');
legend(ax, [h; hp; hx; hmo; hm], {'Leg 1', 'Leg 2', 'Leg 3', 'Leg 4', 'Leg end points', ...
    'Occultation midpoint X_2', 'Moon at the occultation epoch', 'Moon orbit'}, 'Location', 'southwest');
xlabel(ax, 'x [10^3 km]  (ECLIPJ2000, Earth-centred)'); ylabel(ax, 'y [10^3 km]');
title(ax, 'CR3BP orbit propagated in the Sun–Earth–Moon ephemeris model');
exportgraphics(fig1, fullfile(here, '..', 'figures', 'ephemeris_patched.png'), 'Resolution', 150);

cspice_kclear;
