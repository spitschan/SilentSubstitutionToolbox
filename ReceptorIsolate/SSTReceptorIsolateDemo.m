% SSTReceptor based demo of ReceptorIsolate
%
% Demonstration of the ReceptorIsolate engine for finding modulations that
% isolate various receptors; the receptor fundamentals in this
% demonstration are generated using the SSTReceptor class(es).
%
% NOTES:
%  NOTE - JV: the SSTReceptorHuman class does not (currently - 10/30/17) correct
%  for photopigment bleaching, unlike the legacy GetHumanPhotoreceptorSS
%  code that the legacy ReceptorIsolateDemo used.
%
%  NOTE - JV: the SSTReceptorHuman class (currently - 10/30/17) uses fixed
%  standard parameters for vessel oxygenation fraction (.85), and vessel
%  overall thickness (5um).
%
% KNOWN BUGS:
%  Seems to only work for S = [380 2 201]

%% Setup
% Clear workspace and command window, close plots
clearvars; close all; clc;
CheckSSTInstalled();

%% Get receptor fundamentals using SSTReceptor object
% Select receptor model Description of models can be added to this cell
% array; all will be presented as options.
availableModels = {'Human cones, penumbral cones, rods, and melanopsin',...
    };

% Print each entry
fprintf('<strong>Available receptor models:</strong>\n');
for i = 1:numel(availableModels)
    fprintf('\t[%i] %s\n',i,availableModels{i})
end

% Prompt, and check input
receptorModel = GetWithDefault('Enter model number',1);
assert(receptorModel > 0 && receptorModel <= numel(availableModels),'SST:InputError','Unknown receptor model entered');

% Prompt for values for relevant parameters
age = GetWithDefault('Observer age (years)?', 32);
fieldSize = GetWithDefault('Field size (degrees)?', 27.5);
pupilDiameter = GetWithDefault('Pupil diameter (mm)?', 4.7);
%vesselOxyFraction = GetWithDefault('Oxygenation fraction for vessel hemoglobin?', 0.85);
%vesselOverallThicknessUm = GetWithDefault('Vessel thickness?', 5);
doPenumbralCones = logical(GetWithDefault('Include penumbral cones [1 = yes, 0 = no]?',1));

% Get receptor fundamentals using SSTReceptorHuman object
receptors = SSTReceptorHuman('verbosity','high',...
    'obsAgeInYrs',age,...
    'obsPupilDiameter',pupilDiameter,...
    'fieldSizeDeg',fieldSize,...
    'doPenumbralConesTrueFalse',doPenumbralCones);

%% Plot receptor fundamentals
f1 = figure(); clf;
subplot(2,3,1); plot_quantalIsomerizations = receptors.plotSpectralSensitivities('ax',gca,'whichFormat','T_quantalIsomerizations','saveFig',false);
subplot(2,3,2); plot_quantalAbsorptions = receptors.plotSpectralSensitivities('ax',gca,'whichFormat','T_quantalAbsorptions','saveFig',false);
subplot(2,3,3); plot_quantalAbsorptionsNormalized = receptors.plotSpectralSensitivities('ax',gca,'whichFormat','T_quantalAbsorptionsNormalized','saveFig',false);
subplot(2,2,3); plot_energy = receptors.plotSpectralSensitivities('ax',gca,'whichFormat','T_energy','saveFig',false);
subplot(2,2,4); plot_energyNormalized = receptors.plotSpectralSensitivities('ax',gca,'whichFormat','T_energyNormalized','saveFig',false);


%% Define devices
% devices are defined as structs, with the following fields:
%   .name = shorthand name for device
%
%   .description = more detailed (human-readable) description of the device
%   for the users convenience
%
%   .cal = calibration struct for actual devices (e.g., OneLight, monitors)
%   , empty structs for hypothetical devices.
%
%   .S = wavelength sampling, in standard PTB format ( [start delta end] )
%
%   .B_primary = device primary basis vectors, in [0-1] range.
%
%   .ambientSpd = spectral power distribution of ambient light
%
%   .primaryHeadRoom = a constraint that we don't go right to the edge of
%   the gamut.  The head room parameter is defined in the [0-1] device
%   primary space.  Using a little head room keeps us a bit away from the
%   hard edge of the device.
%
%   .maxPowerDiff = a smoothness constraint.

% OneLight
% OneLight device in Brainard/Aguirre lab, with the corresponding
% calibration data. A small headroom is implemented. Smoothness constraint
% value is a magic constant defined by hand, that determins the maximum
% absolute value of the change in spectral power between two adjacent
% wavelength samples. Thus it's appropriate value depends on the overall
% power of the viewed light as well as on the wavelength sampling step.
devices(1).name = 'OneLight';
devices(1).description = 'Brainard/Aguirre lab OneLight device, and corresponding calibration data';
devices(1).cal = LoadCalFile('OneLightDemoCal.mat',[],fullfile(fileparts(mfilename('fullpath')),'ReceptorIsolateDemoData'));
devices(1).S = devices(1).cal.describe.S;
devices(1).B_primary = devices(1).cal.computed.pr650M;
devices(1).ambientSpd = devices(1).cal.computed.pr650MeanDark;
devices(1).primaryHeadRoom = .02;
devices(1).maxPowerDiff = 10^-1.5;

% Ideal spectrum producing device
% Hypothetical ideal spectrum producing device, with delta functions of
% unit power at each wavelength. No headroom or smoothness constraint
% implemented (as the device is ideal).
devices(2).name = 'Spectral';
devices(2).description = 'Hypothetical ideal spectrum producing device (delta function primaries at each wavelength with unit power)';
devices(2).cal = struct(); % hypothetical devices don't have a calibration file
devices(2).S = [380 2 201];
devices(2).B_primary = eye(devices(2).S(3));
devices(2).ambientSpd = zeros(size(devices(2).B_primary,2),1);
devices(2).primaryHeadRoom = 0;
devices(2).maxPowerDiff = 10000;

% Some typical monitor 
% Typical monitor with calibration data (including ambientSpd) supplied by
% the PTB. Headroom implemented. No smoothness constraint implemented, as
% the smoothness of a monitor spectrum is pretty much determined by the
% shape of the primaries.
devices(3).name = 'Monitor';
devices(3).description = 'Some typical monitor (supplied by PTB)';
devices(3).cal = LoadCalFile('PTB3TestCal');
devices(3).S = [380 2 201];
devices(3).B_primary = SplineSpd(devices(3).cal.S_device,devices(3).cal.P_device,devices(3).S);
devices(3).ambientSpd = SplineSpd(devices(3).cal.S_ambient,devices(3).cal.P_ambient,devices(3).S);
devices(3).primaryHeadRoom = 0.02;
devices(3).maxPowerDiff = 10000;

% Five LED primaries
% No smoothness constraint implement, as the smoothness of a monitor
% spectrum is pretty much determined by the shape of the primaries.
devices(4).name = 'FiveLED';
devices(4).description = 'Hypothetical device with five LED primaries';
devices(4).cal = struct(); % hypothetical devices don't have a calibration file
devices(4).S = [380 2 201];
devices(4).B_primary = LEDPrimaries(devices(4).S,[456 488 540 592 632],[10 10 10 17 17]/2,[.57 .125 .156 .27 .75]);
devices(4).ambientSpd = zeros(devices(4).S(3),1);
devices(4).primaryHeadRoom = 0;
devices(4).maxPowerDiff = 10000;

%% Select device
% Print each entry
fprintf('\n<strong>Available devices:</strong>\n');
for i = 1:numel(devices)
    fprintf('<strong>\t[%i] %s</strong>: %s\n',i,devices(i).name,devices(i).description)
end

% Prompt, and check input
input = GetWithDefault('Enter device number',1);
assert(input > 0 && input <= numel(devices),...
    'SST:InputError','Unknown device entered: %i. Range [1-%i] expected',input,numel(devices));
device = devices(input);
fprintf('<strong>Device [%i] selected:</strong> %s (%s)\n',input,device.name,device.description);

%% Set background spectrum
% Half-on in whatever primary space we are working in
backgroundPrimary = repmat(.5,size(device.B_primary,2),1);

%% Correct for bleaching FIXME
% correctBleaching = logical( GetWithDefault('Correct for photopigment bleaching [1 = yes, 0 = no]?',1) );

%% Prompt for direction of modulation
receptorSelection = [receptors.labels',cell(size(receptors.labels'))];
fprintf('\n<strong>Target [1], Silence [2], or Ignore [3] receptor:</strong>\n')
for i = 1:size(receptorSelection,1)
    receptorSelection{i,2} = GetWithDefault(sprintf('\t%s',receptorSelection{i,1}),3);
    assert(receptorSelection{i,2} > 0 && receptorSelection{i,2} <= 3,...
        'SST:InputError','Receptors should either be targeted [1], silenced [2], or Ignored [3]');
end
targetReceptors = find([receptorSelection{:,2}] == 1);
silenceReceptors = find([receptorSelection{:,2}] == 2);
ignoredReceptors = find([receptorSelection{:,2}] == 3);

fprintf('<strong>Seeking direction that isolates</strong>')
fprintf(' %ss,',receptors.labels{targetReceptors});
fprintf('\n<strong>while silencing</strong>')
fprintf(' %ss,',receptors.labels{silenceReceptors});
fprintf('...');

%% Maximize contrast vs. desired contrast TODO
desiredContrast = [];

%% Create direction spectrum (by calling ReceptorIsolate)
directionPrimary = ReceptorIsolate(receptors.T.T_energyNormalized,targetReceptors,ignoredReceptors,[],...
    device.B_primary, backgroundPrimary, backgroundPrimary, [],...
    device.primaryHeadRoom, device.maxPowerDiff, desiredContrast, device.ambientSpd);
fprintf('done.\n');

%% Calculate nominal contrasts
backgroundReceptor = receptors.T.T_quantalIsomerizations * (device.B_primary * backgroundPrimary + device.ambientSpd);
directionReceptor = receptors.T.T_quantalIsomerizations * (device.B_primary * directionPrimary + device.ambientSpd);
contrasts = (directionReceptor - backgroundReceptor) ./ backgroundReceptor;
fprintf('\n<strong>Nominal contrasts:</strong>\n')
for i = 1:numel(receptors.labels)
    fprintf('<strong>\t%s:</strong> %.2f%%\n',receptors.labels{i},contrasts(i)*100);
end

%% Generate some plots
% Plot modulation spectra
plt_modulationSpectra = figure; hold on;
plot(SToWls(device.S),device.B_primary*directionPrimary,'r','LineWidth',2);
plot(SToWls(device.S),device.B_primary*backgroundPrimary,'k','LineWidth',2);
title('Modulation spectra');
legend({'Direction','Background'});
xlim([380 780]);
xlabel('Wavelength');
ylabel('Power');
pbaspect([1 1 1]);

% Plot primaries
plt_modulationPrimaries = figure; hold on;
stem(directionPrimary,'r','LineWidth',2);
stem(backgroundPrimary,'k','LineWidth',2);
title('Modulation primaries settings');
legend({'Direction','Background'});
xlim([1 length(backgroundPrimary)]);
ylim([0 1]);
xlabel('Primary number (nominal)');
ylabel('Setting');