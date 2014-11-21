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

%% Load the background and modulation spectra.
% This modulation was designed to stimulate melanopsin at 45%, for a 32 year old
% observer.
bgSpd = load(fullfile(dataDir, 'spd_background.mat'));
melIsolatingSpd = load(fullfile(dataDir, 'spd_melIsolatingSpd.mat'));
% The order follows that of photoreceptorClasses
targetContrast = [0 0 0 0.45];

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
    melIsolatingSpd.spd, photoreceptorClasses, fieldSizeDegrees, [], pupilDiameterMm, [], fractionBleached)

%% Plot the splatter.  Note that some of the panels may saturate, but we
% don't really mind for demo purposes.
SAVEPLOTS = 1;
theFig = figure;
theFig = PlotSplatter(theFig, contrastMap, photoreceptorClasses, nominalLambdaMax, ...
    observerAgeInYears, ageRange, lambdaMaxShiftRange, targetContrast, [], 1, 1, SAVEPLOTS, '', ...
    'SplatterToolboxDemo_plot', observerAgeInYears);