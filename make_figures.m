% MAKE_FIGURES  Figures for the README (written to figures/).
%   resonant_orbits.png   Moon-resonant periodic orbits from Periodic.csv
%   occultation_pass.png  one passage through the lunar occultation zone
% Uses cr3bp_eom.m and event_function.m from this repository.

clear; close all;
here = fileparts(mfilename('fullpath'));
outDir = fullfile(here, 'figures');
if ~isfolder(outDir), mkdir(outDir); end

% Constants (as in Occultation_events.m)
aE = 149597870.7;  aM = 384400;
TU = 365.25*3600*24/(2*pi);
nE = 1;  Tm = (3600*24*29.53)/TU;  nM = 2*pi/Tm;
mu = 3.0542e-6;
Rs = 695700/aE;  Rs2 = 1.02*Rs;  Rm = 1737.53/aE;
odeopt = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);

family = readmatrix(fullfile(here, 'Periodic.csv'));
family = family(all(~isnan(family(:,1:7)), 2), 1:7);

%% Figure 1: resonant orbits (Earth-centred, Sun-Earth rotating frame)
[~, iu] = unique(round(family(:,7), 6), 'stable');      % one orbit per resonance
fig1 = figure('Color', 'w', 'Position', [100 100 820 760]);
ax = axes(fig1); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on'); axis(ax, 'equal');
cols = lines(numel(iu));
h = gobjects(numel(iu), 1); lbl = cell(numel(iu), 1);
for k = 1:numel(iu)
    X0 = family(iu(k), 1:6); T = family(iu(k), 7);
    [~, X] = ode89(@(t, X) cr3bp_eom(t, X, mu), [0 T], X0, odeopt);
    h(k) = plot(ax, (X(:,1) - (1 - mu))*aE/1e3, X(:,2)*aE/1e3, 'Color', cols(k,:), 'LineWidth', 1.5);
    lbl{k} = sprintf('T_{orbit} / T_{Moon} = %.2f', T/Tm);
end
th = linspace(0, 2*pi, 400);
hM = plot(ax, aM*cos(th)/1e3, aM*sin(th)/1e3, 'k:', 'LineWidth', 1);
plot(ax, 0, 0, 'o', 'MarkerSize', 9, 'MarkerFaceColor', [0.2 0.4 0.9], 'MarkerEdgeColor', 'k');
text(ax, 15, -25, 'Earth');
annotation(fig1, 'textarrow', [0.25 0.15], [0.12 0.12], 'String', 'to Sun');
legend(ax, [h; hM], [lbl; {'Moon orbit'}], 'Location', 'northeastoutside');
xlabel(ax, 'x [10^3 km]  (Sun–Earth rotating frame, Earth-centred)'); ylabel(ax, 'y [10^3 km]');
title(ax, 'Moon-resonant periodic orbits, Sun–Earth CR3BP');
exportgraphics(fig1, fullfile(outDir, 'resonant_orbits.png'), 'Resolution', 150);

%% Figure 2: one occultation passage (orbit 10, Moon phase from the sweep)
X0 = family(10, 1:6);  T = family(10, 7);
thetaM0 = 2.38586508747625;                            % value found by the sweep in Occultation_events.m
ev = @(t, X) event_function(t, X, mu, aM, aE, Rs, Rs2, Rm, nE, nM, thetaM0);
[t, X, te, ~, ie] = ode45(@(t, X) cr3bp_eom(t, X, mu), [0 T], X0, odeset(odeopt, 'Events', ev, 'MaxStep', T/2000));
tIn = te(find(ie == 1, 1));  tOut = te(find(ie == 2, 1));
fprintf('Occultation: entry %.4f, exit %.4f TU, duration %.1f min\n', tIn, tOut, (tOut - tIn)*TU/60);

tMid = (tIn + tOut)/2;
[V1, V2, moon] = localZone(tMid, mu, aM, aE, Rs, Rs2, Rm, nE, nM, thetaM0);
toKm = @(p) [(p(1,:) - (1 - mu)); p(2,:)] * aE/1e3;

fig2 = figure('Color', 'w', 'Position', [100 100 900 650]);
ax = axes(fig2); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on'); axis(ax, 'equal');
p = toKm(X(:,1:3).');
hT = plot(ax, p(1,:), p(2,:), 'Color', [0.2 0.4 0.8], 'LineWidth', 1.2);
pMid = interp1(t, p.', tMid).';
hZ = plot(ax, pMid(1), pMid(2), 'o', 'MarkerSize', 8, 'MarkerFaceColor', [0.9 0.3 0.1], 'MarkerEdgeColor', 'k');
mk = toKm(moon);  v1 = toKm(V1 + [mu;0;0]);  v2 = toKm(V2 + [mu;0;0]);
hS = plot(ax, [mk(1) v1(1)], [mk(2) v1(2)], 'k--', 'LineWidth', 1);   % shadow axis from the Moon
plot(ax, [v1(1) v2(1)], [v1(2) v2(2)], 'k-', 'LineWidth', 2.5);
hMo = plot(ax, mk(1), mk(2), 'o', 'MarkerSize', 8, 'MarkerFaceColor', [0.6 0.6 0.6], 'MarkerEdgeColor', 'k');
plot(ax, 0, 0, 'o', 'MarkerSize', 9, 'MarkerFaceColor', [0.2 0.4 0.9], 'MarkerEdgeColor', 'k');
text(ax, 10, -20, 'Earth');
legend(ax, [hT hZ hS hMo], {'Resonant orbit (one period)', sprintf('Occultation (%.1f min)', (tOut - tIn)*TU/60), ...
    'Sun–Moon axis to the zone vertices V_1, V_2', 'Moon at mid-passage'}, 'Location', 'southoutside');
xlabel(ax, 'x [10^3 km]'); ylabel(ax, 'y [10^3 km]');
title(ax, sprintf('Occultation passage: %.1f min inside the zone (\\theta_{M0} = %.2f°)', ...
    (tOut - tIn)*TU/60, rad2deg(thetaM0)));
% Zoomed inset on the occultation zone
axIn = axes(fig2, 'Position', [0.20 0.30 0.28 0.30]); hold(axIn, 'on'); box(axIn, 'on'); grid(axIn, 'on');
c = (v1 + v2)/2;  w = max(3*norm(v1 - v2), 2);               % window half-width [10^3 km]
near = abs(p(1,:) - c(1)) < w & abs(p(2,:) - c(2)) < w;
plot(axIn, p(1,near), p(2,near), 'Color', [0.2 0.4 0.8], 'LineWidth', 1.2);
plot(axIn, pMid(1), pMid(2), 'o', 'MarkerSize', 8, 'MarkerFaceColor', [0.9 0.3 0.1], 'MarkerEdgeColor', 'k');
plot(axIn, [v1(1) v2(1)], [v1(2) v2(2)], 'k-', 'LineWidth', 2.5);
plot(axIn, [mk(1) v1(1)], [mk(2) v1(2)], 'k--', 'LineWidth', 1);
text(axIn, v1(1), v1(2), '  V_1', 'FontSize', 9); text(axIn, v2(1), v2(2), '  V_2', 'FontSize', 9);
axis(axIn, 'equal'); axis(axIn, [c(1) - w, c(1) + w, c(2) - w, c(2) + w]);
title(axIn, 'zoom on the zone', 'FontWeight', 'normal', 'FontSize', 9);
exportgraphics(fig2, fullfile(outDir, 'occultation_pass.png'), 'Resolution', 150);

function [V1, V2, moon] = localZone(t, mu, aM, aE, Rs, Rs2, Rm, nE, nM, thetaM0)
% Same geometry as event_function.m: Moon position and occultation-zone vertices.
    thetaE = nE*t;  thetaM = nM*t + thetaM0;
    DCM = [cos(thetaE - thetaM), sin(thetaE - thetaM), 0; -sin(thetaE - thetaM), cos(thetaE - thetaM), 0; 0 0 1];
    moon = [1 - mu; 0; 0] + DCM*[aM/aE; 0; 0];
    dSM = norm(moon);  Xh = moon/dSM;
    l1 = Rm*dSM/(Rs2 - Rm);  l2 = Rm*dSM/(Rs - Rm);
    V1 = moon + l2*Xh - [mu; 0; 0];
    V2 = moon + l1*Xh - [mu; 0; 0];
end
