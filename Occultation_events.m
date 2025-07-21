tic;
clc;
clear;

%This script iterates through initial values of the Moon's true anomaly
%every 0.05° and numerically integrates the trajectory with the inclusion
%of the event function. The true anomalies for which an event is triggered,
%along with the states and times at the events are then stored.

%% === Constants (ND) ===
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

%% === Load periodic orbit family ===
family = readmatrix(fullfile(fileparts(mfilename('fullpath')), 'Initial conditions', 'Periodic.csv'));
%% === Moon initial position sweep ===
thetaM0_values = [2.38289802774786, 2.38586508747625, 2.38883214720464, 2.45131493442604, 2.45428199415443];
% 0:pi/36000:2*pi; 
% Initial moon positions iterating every 0.005°
num_moon_positions = length(thetaM0_values);
% Storage for results
occultation_moon_positions = [];
event_data = struct();
% Initialize waitbar
h = waitbar(0, 'Moon position sweep in progress...');
% Update waitbar only every 100 iterations to reduce overhead
update_interval = max(1, round(num_moon_positions / 100));

for moon_idx = 1:num_moon_positions
    current_thetaM0 = thetaM0_values(moon_idx);
    
    % Update waitbar less frequently
    if mod(moon_idx, update_interval) == 0
        waitbar(moon_idx / num_moon_positions, h, ...
            sprintf(' %d of %d', moon_idx, num_moon_positions));
    end
    
    X0 = family(10, 1:6);
    period = family(10, 7);
    tspan = [0, period];
    odeopt = odeset('RelTol',1e-12, 'AbsTol',1e-14);
    
    % Combined event function for both entry and exit
    event = @(t, X) event_function(t, X, mu, aM, aE, Rs, Rs2, Rm, nE, nM, current_thetaM0);
    opts = odeset(odeopt, 'Events', event);
    
    % Integrate with event detection
    [~, ~, te, Xe, ie] = ode45(@(t, X) cr3bp_eom(t, X, mu), tspan, X0, opts);
    
    % If any events detected, store this moon position and event data
    if ~isempty(te)
        fprintf('Event detected at thetaM0 = %.14f rad (%.4f deg)\n', ...
                current_thetaM0, rad2deg(current_thetaM0));
        occultation_moon_positions(end+1) = current_thetaM0;
        event_data(length(occultation_moon_positions)).moon_position = current_thetaM0;
        event_data(length(occultation_moon_positions)).event_times = te;
        event_data(length(occultation_moon_positions)).event_states = Xe;
        event_data(length(occultation_moon_positions)).event_types = ie; % 1=entry, 2=exit
        entry_indices = find(ie == 1);
        exit_indices = find(ie == 2);
        event_data(length(occultation_moon_positions)).entry_times = te(entry_indices);
        event_data(length(occultation_moon_positions)).entry_states = Xe(entry_indices, :);
        event_data(length(occultation_moon_positions)).exit_times = te(exit_indices);
        event_data(length(occultation_moon_positions)).exit_states = Xe(exit_indices, :);
        fprintf('Duration of the event = %.14f s \n', ...
                (te(2) - te(1))*TU);
    end
end

close(h);

toc;
