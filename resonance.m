%This script plots the Jacobi constants vs the period in normalized units
%for the different periodic orbit families in the Sun Earth system.
%Different resonance ratios are plotted as vertical lines and the Orbit IDs
%of resonant orbits are displayed.

% Define folder and CSV file names
dataFolder = 'Initial conditions/';
csvFiles = {'aL1.csv', 'aL2.csv', 'aL3.csv', 'bN.csv', 'bS.csv', 'dN.csv', 'dS.csv', ...
            'dpo.csv', 'dro.csv', 'hL1.csv', 'hL2.csv', 'hL3.csv', 'lpoE.csv', 'lpoW.csv', ...
            'lyL1.csv', 'lyL2.csv', 'lyL3.csv', 'vL1.csv', 'vL2.csv', 'vL3.csv'};
% csvFiles = {'lpoW.csv'};

scatterColors = [0 1 1; 1 0 1; 0 1 0; 1 0 0; 0 0 1; 1 0.5 0; 1 1 0; 0.5 0 0.5; 0.5 1 0.5; 0 0.5 1];
scatterColors = repmat(scatterColors, ceil(length(csvFiles)/size(scatterColors,1)), 1);  

% Preallocate arrays
numFiles = length(csvFiles);
jacobiConstants = [];
periods = [];
familyLabels = [];
orbitIDs = [];

% Define Moon's orbital period and resonance ratios
moonOrbitalPeriod = 0.507987575964444; % 1TU = 1 year/2pi -> nE = 1 TU so TE = 2pi [rad/ND]] and TM [ND] = 29.53/1TU
resonanceRatios = [1, 2, 3, 4, 3/2, 0.5, 1/3, 2/3, 0.25];  
resonanceLabels = {'1:1', '2:1', '3:1', '4:1', '3:2', '1:2', '1:3', '2:3', '1:4'};

% Store resonance results
resonant_orbits = cell(numFiles, length(resonanceRatios));

% Load data from each CSV file
for i = 1:numFiles
    fullPath = fullfile(dataFolder, csvFiles{i});
    data = readmatrix(fullPath); % Read the file
    
    % Extract data columns
    fileIDs = data(:, 1);      % Orbit IDs (Column 1)
    fileJacobi = data(:, 8);   % Jacobi Constants (Column 8)
    filePeriods = data(:, 9);  % Periods (Column 9)

    % Store in main arrays
    jacobiConstants = [jacobiConstants; fileJacobi];
    periods = [periods; filePeriods];
    familyLabels = [familyLabels; repmat(i, size(fileJacobi, 1), 1)];
    orbitIDs = [orbitIDs; fileIDs];

    % Check for resonance
    for j = 1:length(resonanceRatios)
        res_period = moonOrbitalPeriod * resonanceRatios(j);
        match_idx = abs(filePeriods - res_period) < 0.001; % Find matching periods
        
        if any(match_idx)
            resonant_orbits{i, j} = fileIDs(match_idx); % Store matching orbit IDs
        end
    end
end

% Display Resonant Orbit IDs
fprintf('\nResonant Orbits:\n');
for i = 1:numFiles
    for j = 1:length(resonanceRatios)
        if ~isempty(resonant_orbits{i, j})
            fprintf('%s (%s): %s\n', csvFiles{i}, resonanceLabels{j}, num2str(resonant_orbits{i, j}'));
        end
    end
end

% Plot Jacobi Constant vs Orbital Period
figure; hold on;
scatterHandles = gscatter(periods, jacobiConstants, familyLabels, scatterColors, '.', 10);

xlabel('Orbital Period');
ylabel('Jacobi Constant');
title('Jacobi Constant vs Orbital Period');
grid on;
xlim([moonOrbitalPeriod, 4 * moonOrbitalPeriod]);

% Plot resonance lines
resonanceColors = lines(length(resonanceRatios)); 
resonanceHandles = gobjects(length(resonanceRatios), 1);
for i = 1:length(resonanceRatios)
    resonancePeriod = moonOrbitalPeriod * resonanceRatios(i);
    resonanceHandles(i) = xline(resonancePeriod, '-', resonanceLabels{i}, ...
        'Color', resonanceColors(i, :), 'LineWidth', 2, 'FontSize', 10, 'FontWeight', 'bold');
end

legend([scatterHandles; resonanceHandles], [csvFiles, resonanceLabels], 'Location', 'northeastoutside');
hold off;

