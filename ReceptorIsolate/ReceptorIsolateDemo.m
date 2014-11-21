% ReceptorIsolateDemo
%
% Let's find a modulation that isolates various
% photopigments, for various device models.
%
% 3/30/12  dhb      Wrote it.
% 2/3/13   ms       Updated to match changed ReceptorIsolate syntax.
% 4/19/13  dhb, ms  Got this working again.

%% Clear and close
clear; close all;

%% Run us in our home directory
cd(fileparts(mfilename('fullpath')));

%% Which to do?
fprintf('Available receptor models:\n');
fprintf('\t[1]  Melanopsin (human cones and melanopsin)\n');
fprintf('\t[2]  Dog\n');
whichModelNumber = GetWithDefault('Enter model',1);
switch (whichModelNumber)
    case 1
        whichModel = 'Melanopsin';
    case 2
        whichModel = 'Dog';
    otherwise
        error('Unknown primaries entered');         
end

%% Prompt for device to compute with respect to
fprintf('\nAvailable devices:\n');
fprintf('\t[1]  Spectral - ideal spectral producing device\n');
fprintf('\t[2]  OneLight - calibration data from our OneLight\n');
fprintf('\t[3]  Monitor - some typical monitor\n');
whichPrimaryNumber = GetWithDefault('Enter device',1);
switch (whichPrimaryNumber)
    case 1
        whichPrimaries = 'Spectral';
    case 2
        whichPrimaries = 'OneLight';
    case 3
        whichPrimaries = 'Monitor';
    otherwise
        error('Unknown primaries entered');         
end

%% Define primaries and conditions on them
switch (whichPrimaries)
    case 'Spectral'
        % Shows computations in the spectral domain, for
        % a fictional device that can produce delta function
        % primaries at each wavelength with unit power.
        S = WlsToS((400:4:700)');
        B_primary = eye(S(3));
        backgroundPrimary = 0.5*ones(size(B_primary,2),1);   
        
        % Pin first and last primaries at their background value.
        whichPrimariesToPin = [1 size(B_primary,1)];
        
        % Leave no headroom for this ideal device
        primaryHeadRoom = 0;
        
        % Set ambient to zero for this ideal device
        ambientSpd = zeros(size(B_primary,2),1); 
        
        % No smoothness
        maxPowerDiff = Inf;

     case 'OneLight'
        % Get some OneLight primary basis
        cal = LoadCalFile('OLEyeTrackerLongCable');
        S = cal.describe.S;
        B_primary = cal.computed.pr650M;
        ambientSpd = cal.computed.pr650MeanDark;
        
        % Half on in OneLight primary space
        backgroundPrimary = 0.5*ones(size(B_primary,2),1);   
        
        % Don't pin
        whichPrimariesToPin = [];
        primaryHeadRoom = 0.02;
        
        % No smoothness
        maxPowerDiff = Inf;
        
    case 'Monitor'
        S = WlsToS((400:4:700)');
        cal = LoadCalFile('DogScreen1NoLights');
        B_primary = SplineSpd(cal.S_device,cal.P_device,S);
        backgroundPrimary = [0.5 0.5 0.5]';
        ambientSpd = SplineSpd(cal.S_ambient,cal.P_ambient,S);

        whichPrimariesToPin = [];
        primaryHeadRoom = 0;
                
        % No smoothness
        maxPowerDiff = Inf;
end

%% Get sensitivities and set other relvant parameters
switch (whichModel)
    case 'Melanopsin';
        
        % Make sensitivities for L, M, S, Mel
        %
        % ComputeCIEConeFundamentals has subtle behavior
        % with respect to size of lambdaMax passed to it,
        % see the help text.  For that reason, we
        % need to call twice with and piece together the
        % answer.
        fieldSizeDeg = 10;
        whichNomogram = 'StockmanSharpe';
        lambdaMax1 = [558.9, 530.3, 420.7];
        lambdaMax2 = [558.9, 530.3, 480];
        T_quanta1 = ComputeCIEConeFundamentals(S,fieldSizeDeg,32,3,lambdaMax1,whichNomogram);
        T_energy1 = EnergyToQuanta(S,T_quanta1')';
        T_quanta2 = ComputeCIEConeFundamentals(S,fieldSizeDeg,32,3,lambdaMax2,whichNomogram);
        T_energy2 = EnergyToQuanta(S,T_quanta2')';
        T_receptors = [T_energy1(1:3,:) ; T_energy2(3,:)];
        receptorNames = {'L cones', 'M cones', 'S cones', 'Melanopsin'};
        
        % Which receptor to ignore?
        whichReceptorsToIgnore = [];
        
         % Desired contrast in isolated direction
        desiredContrast = 0.25;
    case 'Dog'; 
        % Dog cone and rod receptors
        load T_dogrec
        T_receptors = SplineCmf(S_dogrec,T_dogrec,S);
        receptorNames = {'L cones', 'S cones' 'Rods' };
        
        % Which receptor to ignore?
        whichReceptorsToIgnore = [];
        
        % Desired contrast in isolated direction
        desiredContrast = [];
end

%% Normalize receptors
for i = 1:size(T_receptors,1)
    T_receptors(i,:) = T_receptors(i,:)/max(T_receptors(i,:));
end

%% Prompt for which to isolate
fprintf('\n%s, %s, isolation options:\n',whichModel,whichPrimaries);
for i = 1:length(receptorNames)
    fprintf('\t[%d] %s\n',i,receptorNames{i});
end
whichReceptorsToIsolate = GetWithDefault('Enter which receptor to isolate',1);

%% Background spd.  Make sure is within primaries.
% Need to make sure we start optimization at background,
% or else the constraints don't work so well.
backgroundReceptors = T_receptors*(B_primary*backgroundPrimary + ambientSpd);

%% Properly expand desiredContrast if necessary
if (~isempty(desiredContrast))
    desiredContrast = desiredContrast*ones(size(whichReceptorsToIsolate));
end

%% Set up the background as the initial guess to start the optimization from.
% There are cases in which the initial guess will not be the background.
initialPrimary = backgroundPrimary;

%% Do it, pruit
isolatingPrimary = ReceptorIsolate(T_receptors,whichReceptorsToIsolate,whichReceptorsToIgnore,B_primary,backgroundPrimary,initialPrimary,whichPrimariesToPin,primaryHeadRoom,maxPowerDiff,desiredContrast,ambientSpd);

%% Report responses
diffReceptors = T_receptors*B_primary*(isolatingPrimary - backgroundPrimary);
contrastReceptors = diffReceptors ./ backgroundReceptors;
for i = 1:size(T_receptors,1)
    fprintf('%s: contrast = %0.3f\n',receptorNames{i},contrastReceptors(i));
end

%% Plot
plotDir = 'ReceptorIsolateDemoPlots';
if ~isdir(plotDir);
   mkdir(plotDir); 
end
curDir = pwd;
cd(plotDir);

% Sensitivities
theFig1 = figure; clf; hold on
plot(SToWls(S),T_receptors);
savefig(sprintf('%s_%s_%s_Sensitivities.pdf',whichModel,whichPrimaries,receptorNames{whichReceptorsToIsolate}),theFig1,'pdf');

% Modulation spectra
theFig2 = figure; hold on
plot(SToWls(S),B_primary*isolatingPrimary,'r','LineWidth',2);
plot(SToWls(S),B_primary*backgroundPrimary,'k','LineWidth',2);
title(sprintf('%s, %s, Isolating: %s; isolated contrast: %0.1f',whichModel,whichPrimaries, receptorNames{whichReceptorsToIsolate},contrastReceptors(whichReceptorsToIsolate)));
xlim([380 780]);
xlabel('Wavelength');
ylabel('Power');
pbaspect([1 1 1]);
savefig(sprintf('%s_%s_%s_Modulation.pdf',whichModel,whichPrimaries,receptorNames{whichReceptorsToIsolate}),theFig2,'pdf');

% Primaries
theFig3 = figure; hold on
plot(isolatingPrimary,'r','LineWidth',2);
plot(backgroundPrimary,'k','LineWidth',2);
title(sprintf('%s, %s, Isolating: %s; primary settings',whichModel,whichPrimaries, receptorNames{whichReceptorsToIsolate}));
xlim([0 length(backgroundPrimary)]);
ylim([0 1]);
xlabel('Primary');
ylabel('Setting');
savefig(sprintf('%s_%s_%s_Primaries.pdf',whichModel,whichPrimaries,receptorNames{whichReceptorsToIsolate}),'pdf');

cd(curDir);
