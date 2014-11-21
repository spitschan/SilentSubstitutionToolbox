% SplatterToolboxDemo.m
%
% Demos the SplatterToolbox functions.
%
% See also:
%   GetHumanPhotopigmentSS CalculateSplatter PlotSplatter SaveSplatter
%   GetConeFractionBleachedFromSpectrum
%
% 1/21/14   ms      Wrote it. 11/21/14  ms      Claned up and commented.

% Clear
clear; close all;

%% First, set up some parameters.
S = [380 2 201];            % Wavelength spacing
wls = SToWls(S);
photoreceptorClasses = {'LCone' ; 'MCone' ; 'SCone' ; 'Melanopsin'}; % Photopigments we're interested in
fieldSizeDegrees = 27.5;    % Size of visual field
observerAgeInYears = 32;    % Observer age, assumed to be 32 for this modulation
pupilDiameterMm = 3;        % Pupil diameter, assumed to be 3 mm

%% Figure out where the data live.
dataDir = fullfile(fileparts(which(mfilename)), 'data');

% Make an output folder
outputDir = fullfile(fileparts(which(mfilename)), 'output');
if ~isdir(outputDir)
    mkdir(outputDir);
end

%% Load the background and modulation spectra.
% This modulation was designed to stimulate melanopsin at 45%, for a 32 year old
% observer.
bgSpd = load(fullfile(dataDir, 'spd_background.mat'));
melIsolatingSpd = load(fullfile(dataDir, 'spd_melIsolatingSpd.mat'));
% The order follows that of photoreceptorClasses
targetContrasts = {0 0 0 0.2};

%% Plot the spectra
spectraFig = figure;
plot(bgSpd.wls, bgSpd.spd, '--k'); hold on;
plot(melIsolatingSpd.wls, melIsolatingSpd.spd, '-b');
legend('Background', 'Modulation'); legend boxoff;
xlabel('Wavelength [nm]'); ylabel('Power');

%% Get the photopigment fraction bleached from the background spectrum
[fractionBleachedFromIsom, fractionBleachedFromIsomHemo] = GetConeFractionBleachedFromSpectrum(S, bgSpd.spd, ...
    fieldSizeDegrees, observerAgeInYears, pupilDiameterMm, []);

% We can now assign the fraction bleached for each photoreceptor class.
for p = 1:length(photoreceptorClasses)
    switch photoreceptorClasses{p}
        case 'LCone'
            fractionBleached(p) = fractionBleachedFromIsom(1);
        case 'MCone'
            fractionBleached(p) = fractionBleachedFromIsom(2);
        case 'SCone'
            fractionBleached(p) = fractionBleachedFromIsom(3);
        case 'LConeHemo'
            fractionBleached(p) = fractionBleachedFromIsomHemo(1);
        case 'MConeHemo'
            fractionBleached(p) = fractionBleachedFromIsomHemo(2);
        case 'SConeHemo'
            fractionBleached(p) = fractionBleachedFromIsomHemo(3);
        otherwise
            fractionBleached(p) = 0;
    end
end

%% Calculate the splatter.
% We do not pass age range, fieldSize, or pupilSize so CalculateSplatter
% will take defaults.  To explore effect of fieldSize or pupilSize, you
% need to call CalculateSplatter multiple times and compare the output.
[contrastMap, nominalLambdaMax, ageRange, lambdaMaxShiftRange] = CalculateSplatter(S, bgSpd.spd, ...
    melIsolatingSpd.spd, photoreceptorClasses, fieldSizeDegrees, [], pupilDiameterMm, [], fractionBleached);

%% Plot the splatter.
SAVEPLOTS = true;
theFig = figure;
theFig = PlotSplatter(theFig, contrastMap, photoreceptorClasses, nominalLambdaMax, ...
    observerAgeInYears, ageRange, lambdaMaxShiftRange, targetContrasts, [], 1, 1, SAVEPLOTS, '', ...
    fullfile(outputDir, 'SplatterDemo_plot'), observerAgeInYears);

%% Save some useful numbers.
SaveSplatter(outputDir, 'Demo', contrastMap, photoreceptorClasses, nominalLambdaMax, observerAgeInYears, ageRange, lambdaMaxShiftRange, targetContrasts);

%% Save the CI for all possible observers.
SaveSplatterConfidenceBounds(outputDir, 'Demo_95CI', contrastMap, photoreceptorClasses, nominalLambdaMax, ageRange, lambdaMaxShiftRange, targetContrasts, 0.9545);
SaveSplatterConfidenceBounds(outputDir, 'Demo_99CI', contrastMap, photoreceptorClasses, nominalLambdaMax, ageRange, lambdaMaxShiftRange, targetContrasts, 0.9973);