% ContrastSplatterDemo
%
% Demos the contrast splatter calculation functions.  This runs fairly
% slowly, so it may take > 10 minutes to run.
%
% See also:
%   GetHumanPhotoreceptorSS CalculateSplatter PlotSplatter SaveSplatter
%   GetConeFractionBleachedFromSpectrum
%
% 1/21/14   ms      Wrote it.
% 11/21/14  ms      Claned up and commented.
% 12/12/14  dhb     Tuning.

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
dataDir = fullfile(fileparts(which(mfilename)), 'ContrastSplatterDemoData');

% Make an output folder
outputDir = fullfile(fileparts(which(mfilename)), 'ContrastSplatterDemoOutput');
if ~isdir(outputDir)
    mkdir(outputDir);
end

%% Load the background and modulation spectra.
% The modulation was designed to isolate melanopsin at 45%,
% for a 32 year old observer.

bgSpd = load(fullfile(dataDir,'spd_contrastsplatterdemo_bg.mat'));
melIsolatingSpd = load(fullfile(dataDir, 'spd_contrastsplatterdemo_mod.mat'));

% Specify what the isolating modulation was supposed to do on the 4
% photoreceptor classes considered above.
targetContrasts = {0 0 0 0.2};

%% Plot the bakcground and modulation spectra
spectraFig = figure;
plot(bgSpd.wls, bgSpd.spd, '--k'); hold on;
plot(melIsolatingSpd.wls, melIsolatingSpd.spd, '-b');
legend('Background', 'Modulation'); legend boxoff;
xlabel('Wavelength [nm]'); ylabel('Power');

%% Get the photopigment fraction bleached from the background spectrum
%
% This calculation does not worry about changes in fraction bleached with
% variation in photoreceptor properties, which we think is OK because
% bleaching is generally a calculation sensitive to log units of variation in
% isomerization rates.
%
% Note also that the routine called here does not currently take
% photopigment lambda-max as an argument, so that if one wanted to take
% that factor into account, we'd have to modify that routine.  It does
% take field size, observer age, and pupil diameter, so that we could vary
% those if we chose to.
[fractionBleachedFromIsom, fractionBleachedFromIsomHemo] = GetConeFractionBleachedFromSpectrum(S, bgSpd.spd, ...
    fieldSizeDegrees, observerAgeInYears, pupilDiameterMm, []);
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
% Function ContrastSplatter factorially comptues for age range and lambda-max shifts.
% Here we call it with syntax taking its age range and lambda-max shift range defaults, by
% passing the empty matrix for those two arguments.
%
% To explore effect of fieldSize, pupilSize, or fraction pigment bleached, you
% need to call CalculateSplatter multiple times and compare the output.  We
% don't demonstrate that here.
[contrastMap, nominalLambdaMax, ageRange, lambdaMaxShiftRange] = CalculateSplatter(S, bgSpd.spd, ...
    melIsolatingSpd.spd, photoreceptorClasses, fieldSizeDegrees, [], pupilDiameterMm, [], fractionBleached);

%% Plot the splatter map.
SAVEPLOTS = true;
theFig = figure;
theFig = PlotSplatter(theFig, contrastMap, photoreceptorClasses, nominalLambdaMax, ...
    observerAgeInYears, ageRange, lambdaMaxShiftRange, targetContrasts, [], 1, 1, SAVEPLOTS, '', ...
    fullfile(outputDir, 'SplatterDemo_plot'), observerAgeInYears);

%% Save some useful numbers.
SaveSplatter(outputDir, 'Demo', contrastMap, photoreceptorClasses, nominalLambdaMax, observerAgeInYears, ageRange, lambdaMaxShiftRange, targetContrasts);

%% Save information such as min/max contrast within the confidence ellipses.
SaveSplatterConfidenceBounds(outputDir, 'Demo_95CI', contrastMap, photoreceptorClasses, nominalLambdaMax, ageRange, lambdaMaxShiftRange, targetContrasts, 0.9545);
SaveSplatterConfidenceBounds(outputDir, 'Demo_99CI', contrastMap, photoreceptorClasses, nominalLambdaMax, ageRange, lambdaMaxShiftRange, targetContrasts, 0.9973);