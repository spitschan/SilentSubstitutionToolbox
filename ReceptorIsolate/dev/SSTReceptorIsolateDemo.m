%% Do a little bit of house keeping
clear all; close all; clc

%% Set the wavelength spacing
S = [380 2 201];
wls = SToWls(S);

%% Construct the receptor object
receptorObj = SSTReceptorHuman('S', S, 'verbosity', 'high', 'obsAgeYrs', 32);


%%
% Construct LEDs
NLEDs = 8;
peakWls = [450 472.5 502.5 530 590 615 632.5 660];
fwhm = 12*ones(1, NLEDs);
maxPower = ones(1, NLEDs);
for i = 1:length(fwhm)
    % Figure out the standard deviation.
    standardDeviation(i) = FWHMToStd(fwhm(i));
end


% Make the spectrum.
for i = 1:length(fwhm)
    spd(:, i) = normpdf(wls, peakWls(i), fwhm(i));
    spd(:, i) = spd(:, i)./max(spd(:, i))*maxPower(i);
end

B_primary = spd;

% Set background to the monitor midpoint, and use the ambient
% spectrum from the calibration file.
backgroundPrimary = 0.5*ones(NLEDs, 1);
ambientSpd = wls; ambientSpd(:) = 0;

% Don't pin any primaries.  Do enforce a constraint that we don't
% go right to the edge of the gamut.  The head room parameter is
% defined in the [0-1] device primary space.  Using a little head
% room keeps us a bit away from the hard edge of the device.
whichPrimariesToPin = [];
primaryHeadRoom = 0.02;

% No smoothness constraint envforced here.  It really wouldn't make
% to much sense for a three-primary monitor, as the smoothness of a
% monitor spectrum is pretty much determined by the spectral shape
% of its primarites.
maxPowerDiff = 10000;

%%
% Isolate the receptors by calling the wrapper
% initialPrimary = backgroundPrimary;
% [modulationPrimary backgroundPrimary] = ReceptorIsolateOptimBackgroundMulti(receptorObj.T.T_energy, whichReceptorsToTarget, ...
%     whichReceptorsToIgnore,whichReceptorsToMinimize,B_primary,backgroundPrimary,...
%     initialPrimary,whichPrimariesToPin,primaryHeadRoom,maxPowerDiff,...
%     desiredContrasts,ambientSpd,directionsYoked,directionsYokedAbs,pegBackground);

%backgroundPrimary = modulationPrimary{1};
whichDirection = 'SDirected';
whichReceptorsToTarget = [4];
whichReceptorsToIgnore = [5];
whichReceptorsToMinimize = [];
desiredContrasts = [];
pegBackground = false;

modulationPrimary = ReceptorIsolate(receptorObj.T.T_energy,whichReceptorsToTarget, whichReceptorsToIgnore, whichReceptorsToMinimize, ...
   B_primary, backgroundPrimary, backgroundPrimary, whichPrimariesToPin,...
   primaryHeadRoom, 0.1, desiredContrasts, ambientSpd);
%
% Compute the contrasts that we got.
fprintf('\n');
for ii = 1:size(modulationPrimary, 2)
    backgroundReceptors = receptorObj.T.T_energy*(B_primary*backgroundPrimary + ambientSpd);
    modulationReceptors = receptorObj.T.T_energy*B_primary*(modulationPrimary - backgroundPrimary);
    contrastReceptors = modulationReceptors ./ backgroundReceptors;
    for j = 1:size(receptorObj.T.T_energy,1)
        fprintf('\t%s: contrast = %0.4f\n',receptorObj.labels{j},contrastReceptors(j));
    end
end

modPrimaryPos = modulationPrimary;
modPrimaryNeg = backgroundPrimary - (modulationPrimary - backgroundPrimary);

%%
%
gcFig1G = figure;
plot(wls, B_primary*modPrimaryNeg, 'LineWidth', 2, 'Color', [97 103 100]/255); hold on;
%plot(wls, B_primary*backgroundPrimary, '-', 'LineWidth', 2, 'Color', [108 94 219]/255);
plot(wls, B_primary*modPrimaryPos, 'LineWidth', 2, 'Color', [38 47 219]/255);
xlabel('Wavelength [nm]'); ylabel('Relative power');
set(gca, 'XTick', [400 500 600 700]);
set(gca, 'YTick', [0 1]);
set(gca, 'TickDir', 'out');
ylim([-0.01 1.01]);
xlim([400 750]);
box off;
pbaspect([1 0.5 1]);

% Save the figure
set(gcFig1G, 'PaperPosition', [0 0 4 3]);
set(gcFig1G, 'PaperSize', [4 3]);
cd(fullfile(fileparts(mfilename('fullpath'))));
saveas(gcFig1G, 'Fig1G.pdf', 'pdf');
close(gcFig1G);



%%
maxPowerDiff = 10000;
whichDirection = 'MelDirected';
whichReceptorsToTarget = {[4 5]};
whichReceptorsToIgnore = {[]};
whichReceptorsToMinimize = {[]};
desiredContrasts = {[]};
directionsYoked = [0];
directionsYokedAbs = [0];
pegBackground = false;

% Isolate the receptors by calling the wrapper
initialPrimary = backgroundPrimary;
[modulationPrimary backgroundPrimary] = ReceptorIsolateOptimBackgroundMulti(receptorObj.T.T_energy, whichReceptorsToTarget, ...
    whichReceptorsToIgnore,whichReceptorsToMinimize,B_primary,backgroundPrimary,...
    initialPrimary,whichPrimariesToPin,primaryHeadRoom,maxPowerDiff,...
    desiredContrasts,ambientSpd,directionsYoked,directionsYokedAbs,pegBackground);

backgroundPrimary = modulationPrimary{1};
whichDirection = 'MelDirected';
whichReceptorsToTarget = [4 5];
whichReceptorsToIgnore = [];
whichReceptorsToMinimize = [];
desiredContrasts = [0.5 0.3];
pegBackground = false;

modulationPrimary = ReceptorIsolate(receptorObj.T.T_energy,whichReceptorsToTarget, whichReceptorsToIgnore, whichReceptorsToMinimize, ...
   B_primary, backgroundPrimary, modulationPrimary{2}, whichPrimariesToPin,...
   primaryHeadRoom, maxPowerDiff, desiredContrasts, ambientSpd);

% Compute the contrasts that we got.
for ii = 1:size(modulationPrimary, 2)
    backgroundReceptors = receptorObj.T.T_energy*(B_primary*backgroundPrimary + ambientSpd);
    modulationReceptors = receptorObj.T.T_energy*B_primary*(modulationPrimary - backgroundPrimary);
    contrastReceptors = modulationReceptors ./ backgroundReceptors;
    for j = 1:size(receptorObj.T.T_energy,1)
        fprintf('\t%s: contrast = %0.4f\n',receptorObj.labels{j},contrastReceptors(j));
    end
end

modPrimaryPos = modulationPrimary;
modPrimaryNeg = backgroundPrimary - (modulationPrimary - backgroundPrimary);
