clc;
clear;
close all;
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultAxesFontSize', 20);
set(0, 'DefaultLineLineWidth', 2);


clc;
clear;
close all;

set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultAxesFontSize', 20);
set(0, 'DefaultLineLineWidth', 2);

%% Extracting Events Data for each orbit

% Working directory
data_folder = pwd;

% Number of files to load
num_files = 1;

% Preallocate cell arrays to store separated data
moon_pos  = cell(num_files, 1);
t_entry   = cell(num_files, 1);
t_exit    = cell(num_files, 1);

% Loop over each .mat file
for k = 1:num_files
    filename = fullfile(data_folder, [num2str(k), '.mat']);
    S = load(filename);

    % Extract struct (assumes one variable per file)
    field_name = fieldnames(S);
    results = S.(field_name{1});

    % Extract moon positions, entry and exit times
    moon_pos{k} = [results.moon_position];
    t_entry{k}  = [results.entry_times];
    t_exit{k}   = [results.exit_times];
end
%% === Constants (Non-dimensionalized) ===
aE = 149597870.7; % Earth-Sun distance [km]
aM = 384400; % Moon-Earth distance [km]
TU = 365.25*3600*24/(2*pi); % Characteristic Time length [s]
Te = 2*pi; % Normalised Sun Earth period [ND]
Tm = (3600*24*29.53)/TU; % Normalised Earth moon period [ND]
nE = 2*pi/Te;  % Earth mean motion [rad / TU], 1 TU = 1 year
nM = 2*pi/Tm;  % Moon mean motion [rad / TU]
mu = 3.0542e-6; % Earth-Moon system mu
Rs = 695700 / aE; % Sun radius (normalized)
Rs2 = 1.02 * Rs; % Slightly bigger Sun radius
Rm = 1737.53 / aE; % Moon radius (normalized)
orbit_idx = 1;
thetaM0_ID = 3; % Select index of the chosen simulated moon true anomaly 
thetaM0 = moon_pos{orbit_idx}(thetaM0_ID); % Extract Moon initial true anomaly from corresponding event data struct

%% === Load periodic orbit family ===
file_path = fullfile(fileparts(mfilename('fullpath')), 'Initial conditions', 'Periodic.csv');
family = readmatrix(file_path);
num_orbits = size(family, 1);

%% === Plot setup ===
figure;
hold on;
grid on;
xlabel('X [ND]');
ylabel('Y [ND]');
zlabel('Z [ND]');
title('CR3BP Trajectories and Occultation Zone');
axis equal;
view(3);

%% === ODE Solver Options ===
odeopt = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);

%% === ODE Integration and Plotting for Selected Orbit ===

% Initial conditions for selected orbit
X0 = family(orbit_idx, 1:6);
period = family(orbit_idx, 7);
tspan = [0,(t_exit{orbit_idx}(thetaM0_ID)+ t_entry{orbit_idx}(thetaM0_ID))/2];
% tspan = [0,3*period];

% Integrate trajectory
[t, X] = ode89(@(t, X) cr3bp_eom(t, X, mu), tspan, X0, odeopt);

% Plot trajectory
plot3(X(:,1), X(:,2), X(:,3), 'b-', 'LineWidth', 2);

% Mark start and end points
plot3(X(1,1), X(1,2), X(1,3), 'o', 'MarkerSize', 10, ...
    'MarkerFaceColor', 'cyan', 'MarkerEdgeColor', 'black');
plot3(X(end,1), X(end,2), X(end,3), 's', 'MarkerSize', 10, ...
    'MarkerFaceColor', 'magenta', 'MarkerEdgeColor', 'black');

% Compute and plot V1 and V2 along the trajectory
V1_track = zeros(length(t), 3);
V2_track = zeros(length(t), 3);

for idx = 1:length(t)
    [V1, V2] = occultation_zone(t(idx), mu, aM, aE, Rs, Rs2, Rm, nE, nM, thetaM0);
    V1_track(idx, :) = V1';
    V2_track(idx, :) = V2';
end

% Plot V1 and V2 curves
plot3(V1_track(:,1), V1_track(:,2), V1_track(:,3), 'r:', 'LineWidth', 2);
plot3(V2_track(:,1), V2_track(:,2), V2_track(:,3), 'g:', 'LineWidth', 2);

% Mark start and end of V1
plot3(V1_track(1,1), V1_track(1,2), V1_track(1,3), '^', 'MarkerSize', 8, ...
    'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'black');
plot3(V1_track(end,1), V1_track(end,2), V1_track(end,3), 'v', 'MarkerSize', 8, ...
    'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'black');

% Mark start and end of V2
plot3(V2_track(1,1), V2_track(1,2), V2_track(1,3), '^', 'MarkerSize', 8, ...
    'MarkerFaceColor', 'green', 'MarkerEdgeColor', 'black');
plot3(V2_track(end,1), V2_track(end,2), V2_track(end,3), 'v', 'MarkerSize', 8, ...
    'MarkerFaceColor', 'green', 'MarkerEdgeColor', 'black');

% Connect final V1 and V2 with dashed line
plot3([V1_track(end,1), V2_track(end,1)], ...
      [V1_track(end,2), V2_track(end,2)], ...
      [V1_track(end,3), V2_track(end,3)], ...
      '--', 'Color', 'k', 'LineWidth', 1.5);


% Plot settings
grid on;
legend({'Trajectory', 'Start Point', 'End Point', ...
    'V1 Curve', 'V2 Curve', ...
    'V1 Start', 'V1 End', ...
    'V2 Start', 'V2 End'}, 'Location', 'bestoutside');
view(3);


% Compute Z drift 
Z_drift = (X(end,3) - X(1,3)) * aE;
disp(['Z drift: ', num2str(abs(Z_drift)), ' km']);

%Compute distance from the middle of the occultation zone at the end
P = X(end, 1:3);

V1 = V1_track(end, :);
V2 = V2_track(end, :);

V1V2 = V2 - V1;

V1P = P - V1;

% Projection of V1P onto V1V2
proj_length = dot(V1P, V1V2) / norm(V1V2);

if proj_length <= 0
    closest_point = V1;
elseif proj_length >= norm(V1V2)
    closest_point = V2;
else
    closest_point = V1 + (proj_length / norm(V1V2)) * V1V2;
end

% Distance from spacecraft to V1–V2 line
distance = norm(P - closest_point) * aE; 

% Display result
disp(['Distance from SC to V1–V2 line (final point): ', num2str(distance), ' km']);