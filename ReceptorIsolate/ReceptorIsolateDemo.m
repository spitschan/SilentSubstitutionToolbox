% ReceptorIsolateDemo
%
% Find modulations that isolates various photopigments, for various device
% models. The engine that drives this demo is ReceptorIsolate.
%
% This demonstrates various cases, including when you've just got a
% standard three primary monitor and when you have a device that can
% produce arbitrary (within gamut limits) spectra across the visible range
% of wavelengths.
% 
% It also demonstrates our code for producing individually customized (for
% age, field size, pupil diameter, and background light level) estimates of
% human cone spectral sensitivities.
%
% 3/30/12  dhb      Wrote it.
% 2/3/13   ms       Updated to match changed ReceptorIsolate syntax.
% 4/19/13  dhb, ms  Got this working again.
% 12/11/14 dhb      Clean up for release.

%% Clear and close
clear; close all;

%% Run in the directory that contains this file.
cd(fileparts(mfilename('fullpath')));

%% User choice of which of several cases to do
fprintf('Available receptor models:\n');
fprintf('\t[1]  Human cones and melanopsin\n');
fprintf('\t[2]  Dog\n');
whichModelNumber = GetWithDefault('Enter model',1);
switch (whichModelNumber)
    case 1
        whichModel = 'HumanPhotopigments';
    case 2
        whichModel = 'Dog';
    otherwise
        error('Unknown receptor case entered');
end

%% Prompt for device to compute with respect to
fprintf('\nAvailable devices:\n');
fprintf('\t[1]  Spectral - ideal spectral producing device\n');
fprintf('\t[2]  OneLight - calibration data from Brainard/Aguirre lab OneLight device\n');
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
        % Shows computations in the spectral domain, for a fictional device
        % that can produce delta function primaries at each wavelength with
        % unit power.
        S = WlsToS((400:4:700)');
        B_primary = eye(S(3));
        backgroundPrimary = 0.5*ones(size(B_primary,2),1);
        
        % Pin first and last primaries at their background value during the
        % optimization
        whichPrimariesToPin = [1 size(B_primary,1)];
        
        % Leave no headroom for this ideal device
        primaryHeadRoom = 0;
        
        % Set ambient light to zero for this ideal device
        ambientSpd = zeros(size(B_primary,2),1);
        
        % No smoothness constraint enforced here. Use a big number (Inf
        % does not work).
        maxPowerDiff = 10000;    
    case 'OneLight'
        % Get a OneLight calibration file, stored here for demo purposes.
        % Extract the descrption of spectral primaries, which is what we
        % need for this demo.  Gamma correction of primary values to
        % settings would need to be handled to get the device to actually
        % produce the spectrum, but here we are just finding the desired
        % modulation.
        calPath = fullfile(fileparts(mfilename('fullpath')), 'cals', []);
        cal = LoadCalFile('OneLightDemoCal.mat',[],calPath);
        S = cal.describe.S;
        B_primary = cal.computed.pr650M;
        ambientSpd = cal.computed.pr650MeanDark;
        
        % Half on in OneLight primary space
        backgroundPrimary = 0.5*ones(size(B_primary,2),1);
        
        % Don't pin any primaries.  Do enforce a constraint that we don't
        % go right to the edge of the gamut.  The head room parameter is
        % defined in the [0-1] device primary space.  Using a little head
        % room keeps us a bit away from the hard edge of the device.
        whichPrimariesToPin = [];
        primaryHeadRoom = 0.02;
        
        % Set smoothness constraint value.  This is a magic constant
        % defined by hand, that determins the maximum absolute value of the
        % change in spectral power between two adjacent wavelength samples.
        % Thus it's appropriate value depends on the overall power of the
        % viewed light as well as on the wavelength sampling step.
        maxPowerDiff = 10^-1.5; 
    case 'Monitor'
        % Typical monitor calibration file.  Use the one supplied with PTB.
        S = WlsToS((400:4:700)');
        cal = LoadCalFile('PTB3TestCal');
        B_primary = SplineSpd(cal.S_device,cal.P_device,S);
        
        % Set background to the monitor midpoint, and use the ambient
        % spectrum from the calibration file.
        backgroundPrimary = [0.5 0.5 0.5]';
        ambientSpd = SplineSpd(cal.S_ambient,cal.P_ambient,S);
        
         % Don't pin any primaries.  Do enforce a constraint that we don't
        % go right to the edge of the gamut.  The head room parameter is
        % defined in the [0-1] device primary space.  Using a little head
        % room keeps us a bit away from the hard edge of the device.
        whichPrimariesToPin = [];
        primaryHeadRoom = 0;
        
        % No smoothness constraint envforced here.  It really wouldn't make
        % to much sense for a three-primary monitor, as the smoothness of a
        % monitor spectrum is pretty much determined by the spectral shape
        % of its primarites. 
        maxPowerDiff = 10000;
end

%% Get sensitivities and set other relvant parameters
switch (whichModel)
    case 'HumanPhotopigments';
        % The routines that do these computations are in the
        % ContrastSplatter directory of the SilentSubstitutionToolbox. They
        % provide pre-defined receptor types and compute spectral
        % sensitivities using the routines provided in the Psychtoolbox.
        % The routines here, however, also allow computation of fraction
        % cone bleached, which may be used to adjust pigment peak optical
        % density.  They can also compute photopigment variants corrected
        % for filtering by blood vessels.
        
        % Prompt user for key parameters that affect the spectral
        % sensitivities.
        observerAgeInYears = GetWithDefault('\tObserver age in years?', 32);
        fieldSizeDegrees = GetWithDefault('\tField size in degrees?', 27.5);
        pupilDiameterMm = GetWithDefault('\tPupil diameter?', 4.7);
                correctBleaching = GetWithDefault('\tCorrect for photopigment bleaching [1 = yes, 0 = no]?', 1);
                
        % Define photoreceptor classes that we'll consider.
        photoreceptorClasses = {'LCone', 'MCone', 'SCone', 'Melanopsin', 'Rods', 'LConeHemo', 'MConeHemo', 'SConeHemo'};
        
        % Correct for pigment bleaching if desired.  This is done
        % separately for open-field and penumbral cones.  The bleaching
        % correction routine only knows about L, M, and S cones.
        if correctBleaching
            % The hard work is done by routine
            % GetConeFractionBleachedFromSpectrum, which is in the
            % ContrastSplatter directory. We set the verbose flag here.
            verbose = true;
            [fractionBleachedFromIsom, fractionBleachedFromIsomHemo] = ...
                GetConeFractionBleachedFromSpectrum(S, B_primary*backgroundPrimary + ambientSpd, fieldSizeDegrees, observerAgeInYears, pupilDiameterMm, [], verbose);
            
            % Assign the fraction bleached for each photoreceptor class.
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
        else
            fractionBleached = [];
        end
        
        % Make sensitivities.  The wrapper routine is
        % GetHumanPhotopigmentSS, which is in the ContrastSplatter
        % directory.  Each row of the matrix T_receptors provides the
        % spectral sensitivity of the photoreceptor class in the
        % corresponding entry of the cell array photoreceptorClasses.
        %
        % The last two arguments are the oxygenation fraction and the
        % vessel thickness. We set them to be empty here, prompting the
        % user to enter these values later.
        oxygenationFraction = [];
        vesselThickness = [];
        T_receptors = GetHumanPhotopigmentSS(S, photoreceptorClasses, fieldSizeDegrees, observerAgeInYears, pupilDiameterMm, [], fractionBleached, oxygenationFraction, vesselThickness);
        
        %% Let user choose a photoreceptor class to target
        fprintf('Available photoreceptor classes to target:\n');
        fprintf('\t[1]  Melanopsin, silience open-field cones; ignore rods and penumbral cones\n');
        fprintf('\t[2]  Melanopsin, silence open-field and penumbral cones; ignore rods)\n');
        fprintf('\t[3]  S cones, silence open-field cones, melanopsin, and prenumbral L and M cones; ignore rods and penumbral S cones\n');
        fprintf('\t[4]  Penumbral L and M cones, silence open-field cones, melanopsin, and prenumbral S cones; ignore rods\n');

        whichDirectionNumber = GetWithDefault('Enter direction',1);
        
        % Depending on which direction is chosen, specify the indices
        % into the rows of T_receptors to define the various classes.
        %  Targeted photoreceptors - contrast maximized or driven to specified target contrast
        %  Ignored photoreceptors - ignored in calculation
        %  Minmized photoreceptors - legacy variable no longer used
        %  All remaining photoreceptors are silenced.
        switch (whichDirectionNumber)
            case 1
                whichDirection = 'MelanopsinDirectedLegacy';
                whichReceptorsToTarget = [4];
                whichReceptorsToIgnore = [5 6 7 8];
                whichReceptorsToMinimize = [];
            case 2
                whichDirection = 'MelanopsinDirected';
                whichReceptorsToTarget = [4];
                whichReceptorsToIgnore = [5];
                whichReceptorsToMinimize = [];
            case 3
                whichDirection = 'SDirected';
                whichReceptorsToTarget = [3];
                whichReceptorsToIgnore = [5 8];
                whichReceptorsToMinimize = [];
            case 4
                whichDirection = 'PenumbralLM';
                whichReceptorsToTarget = [6 7];
                whichReceptorsToIgnore = [5];
                whichReceptorsToMinimize = [];
            otherwise
                error('Unknown direction entered');
        end
        
        
    case 'Dog';
        % Load in PTB's view of the photoreceptor spectral sensitivities in
        % canine.
        load T_dogrec
        T_receptors = SplineCmf(S_dogrec,T_dogrec,S);
        photoreceptorClasses = {'L cones', 'S cones' 'Rods' };
        
        % Which receptor to ignore?
        whichReceptorsToIgnore = [];
        
        % Desired contrast in isolated direction
        desiredContrast = [];
        
        %% Which to do?
        fprintf('Available directions:\n');
        fprintf('\t[1]  Dog L cones\n');
        fprintf('\t[2]  Dog S cones\n');
        fprintf('\t[3]  Dog Rods\n');
        whichDirectionNumber = GetWithDefault('Enter direction',1);
        switch (whichDirectionNumber)
            case 1
                whichDirection = 'DogL';
                whichReceptorsToTarget = [1];
                whichReceptorsToIgnore = [];
                whichReceptorsToMinimize = [];
            case 2
                whichDirection = 'DogS';
                whichReceptorsToTarget = [2];
                whichReceptorsToIgnore = [];
                whichReceptorsToMinimize = [];
            case 3
                whichDirection = 'DogRods';
                whichReceptorsToTarget = [3];
                whichReceptorsToIgnore = [];
                whichReceptorsToMinimize = [];
            otherwise
                error('Unknown direction entered');
        end
        
end

% User chooses whether to maximize contrast in targeted receptor classes or
% or get it as close to a specified value as possible.
maximizeTargetContrast = GetWithDefault('\tMaximize contrast? [1 = yes, 0 = no]', 1);
if maximizeTargetContrast
    desiredContrast = [];
else
    desiredContrast = GetWithDefault('\tDesired contrast?', 0.45)*ones(size(whichReceptorsToTarget));
end

% Nice message for user
fprintf('\nGenerating stimuli which isolate receptor classes');
for i = 1:length(whichReceptorsToTarget)
    fprintf('\n  - %s', photoreceptorClasses{whichReceptorsToTarget(i)});
end
fprintf('\nGenerating stimuli which ignore receptor classes');
if (~length(whichReceptorsToIgnore) == 0)
    for i = 1:length(whichReceptorsToIgnore)
        fprintf('\n  - %s', photoreceptorClasses{whichReceptorsToIgnore(i)});
    end
else
    fprintf('\n  - None');
end
fprintf('\nThe remaining classes will be silenced\n');

%% Call the optimization routine.
%
% Careful examaination of the arguments will reveal that the initialGuess
% for the primaries is set to the background value for the primaries, so
% that the constraints are all met at the start of the search.  The
% optimization routine is much happier when it is done this way -- bad
% things happen if you start with a guess that violates constraints.
modulationPrimary = ReceptorIsolate(T_receptors,whichReceptorsToTarget, whichReceptorsToIgnore, whichReceptorsToMinimize, ...
    B_primary, backgroundPrimary, backgroundPrimary, whichPrimariesToPin,...
    primaryHeadRoom, maxPowerDiff, desiredContrast, ambientSpd);

%% Compute the contrasts that we got.
backgroundReceptors = T_receptors*(B_primary*backgroundPrimary + ambientSpd);
modulationReceptors = T_receptors*B_primary*(modulationPrimary - backgroundPrimary);
contrastReceptors = modulationReceptors ./ backgroundReceptors;
for j = 1:size(T_receptors,1)
    fprintf('\t%s: contrast = %0.4f\n',photoreceptorClasses{j},contrastReceptors(j));
end

%% Plots
plotDir = 'ReceptorIsolateDemoPlots';
if ~isdir(plotDir);
    mkdir(plotDir);
end
curDir = pwd;
cd(plotDir);

% Photoreceptor sensitivities
theFig1 = figure; clf; hold on
plot(SToWls(S),T_receptors,'LineWidth',2);
xlabel('Wavelength (nm)')
ylabel('Sensitivity');
title('Normalized photoreceptor sensitivities');
saveas(theFig1,sprintf('%s_%s_%s_Sensitivities.pdf',whichModel,whichPrimaries,whichDirectionNumber),'pdf');

% Modulation spectra
theFig2 = figure; hold on
plot(SToWls(S),B_primary*modulationPrimary,'r','LineWidth',2);
plot(SToWls(S),B_primary*backgroundPrimary,'k','LineWidth',2);
title('Modulation spectra');
xlim([380 780]);
xlabel('Wavelength');
ylabel('Power');
pbaspect([1 1 1]);
saveas(theFig2,sprintf('%s_%s_%s_Modulation.pdf',whichModel,whichPrimaries,whichDirectionNumber),'pdf');

% Primaries
theFig3 = figure; hold on
plot(modulationPrimary,'r','LineWidth',2);
plot(backgroundPrimary,'k','LineWidth',2);
title('Primary settings');
xlim([0 length(backgroundPrimary)]);
ylim([0 1]);
xlabel('Primary Number (nominal)');
ylabel('Setting');
saveas(theFig3,sprintf('%s_%s_%s_Primaries.pdf',whichModel,whichPrimaries,whichDirectionNumber),'pdf');

%% Return to the directory from whence we started
cd(curDir);