% SSTReceptorCheckSS2And10Deg.m
%
% Checks whether the @SSTReceptorHuman object returns the correct spectral
% sensitivities by comparing them against the 2 deg and 10 deg Stockman-Sharpe
% cone fundamentals.
%
% 7/26/17   ms      Written.

% Wavelength sampling
S = [380 5 81]; wls = SToWls(S);

%% 2 deg cone fundamentals
% Load from PTB
targetRaw = load('T_cones_ss2');
T_targetEnergy2 = SplineCmf(targetRaw.S_cones_ss2,targetRaw.T_cones_ss2,S,2);
T_targetQuantal2 = QuantaToEnergy(S,T_targetEnergy2')';

for i = 1:3
    T_targetQuantal2(i,:) = T_targetQuantal2(i,:)/max(T_targetQuantal2(i,:));
end

% Construct with @SSTReceptorHuman. The Stockman-Sharpe 2 deg fundamentals
% assume an observer of 32 years, a pupil diameter of 3 mm, and a field
% size of 2 deg.
observerAgeInYrs = 32;
pupilDiameterMm = 3;
fieldSizeDegrees = 2;
receptorObj = SSTReceptorHuman('S', S, 'obsAgeYrs', observerAgeInYrs, 'fieldSizeDeg', fieldSizeDegrees, 'obsPupilDiameterMm', pupilDiameterMm);

% Plot
subplot(1, 2, 1);
h1 = plot(wls, receptorObj.T.T_quantalAbsorptionsNormalized', '-k', 'LineWidth', 4);
hold on;
h2 = plot(wls, T_targetQuantal2', '-g', 'LineWidth', 1);

% Add legend and tune plot
legend([h1(1) h2(1)], '@SSTReceptorHuman', 'PTB'); legend boxoff;
pbaspect([1 1 1]); set(gca, 'TickDir', 'out');
xlabel('Wavelength [nm]'); ylabel('Normalized quantal sensitivity');
title('2 deg fundamentals');
xlim([380 780]);

%% 10 deg cone fundamentals
% Load from PTB
load T_cones_ss10
targetRaw = load('T_cones_ss10');
T_targetEnergy10 = SplineCmf(targetRaw.S_cones_ss10,targetRaw.T_cones_ss10,S,2);
T_targetQuantal10 = QuantaToEnergy(S,T_targetEnergy10')';

for i = 1:3
    T_targetQuantal10(i,:) = T_targetQuantal10(i,:)/max(T_targetQuantal10(i,:));
end

% Construct with @SSTReceptorHuman. The Stockman-Sharpe 2 deg fundamentals
% assume an observer of 32 years, a pupil diameter of 3 mm, and a field
% size of 10 deg.
observerAgeInYrs = 32;
pupilDiameterMm = 3;
fieldSizeDegrees = 10;
receptorObj = SSTReceptorHuman('S', S, 'obsAgeYrs', observerAgeInYrs, 'fieldSizeDeg', fieldSizeDegrees, 'obsPupilDiameterMm', pupilDiameterMm);

% Plot
subplot(1, 2, 2);
h1 = plot(wls, receptorObj.T.T_quantalAbsorptionsNormalized', '-k', 'LineWidth', 4);
hold on;
h2 = plot(wls, T_targetQuantal10', '-g', 'LineWidth', 1);

% Add legend and tune plot
legend([h1(1) h2(1)], '@SSTReceptorHuman', 'PTB'); legend boxoff;
pbaspect([1 1 1]); set(gca, 'TickDir', 'out');
xlabel('Wavelength [nm]'); ylabel('Normalized quantal sensitivity');
title('10 deg fundamentals');
xlim([380 780]);