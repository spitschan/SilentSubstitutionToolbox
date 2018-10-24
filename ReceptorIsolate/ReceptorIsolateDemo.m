function status = ReceptorIsolateDemo(varargin)
% Demonstrate calculation of receptor isolating directions
%
% Syntax:
%    ReceptorIsolateDemo
%
% Description:
%    Find modulations that isolates various photopigments, for various device
%    models. The engine that drives this demo is called ReceptorIsolate.
%
%    This demonstrates various cases, including when you've just got a
%    standard three primary monitor and when you have a device that can
%    produce arbitrary (within gamut limits) spectra across the visible range
%    of wavelengths.
%
%    It also demonstrates our code for producing individually customized (for
%    age, field size, pupil diameter, and background light level) estimates of
%    human cone spectral sensitivities.
%
% Inputs:
%    None
%
% Outputs:
%    status - Returns 1 in interactive mode, or if validations pass. Returns
%             0 if a validation fails.
%
% Optional key/value pairs:
%    'validate'                - String (default 'none'). Can pass strings to run
%                                preset validations. When string is none, runs in
%                                interactive mode and prompts for values.
%                                 - 'none' Interactive mode, no validation checks.
%                                 - 'basichuman' Basic human validation
%    'whichDirectionNumber'    - Number (default 1). Which preset direction
%                                to validate, when 'validate' is not 'none'
%    'validationFractionTolerance' - Number (default 0.001). Fractional
%                                tolerance for validation check, when validating.
%
% See also:
%

% History:
%   3/30/12  dhb      Wrote it.
%   2/3/13   ms       Updated to match changed ReceptorIsolate syntax.
%   4/19/13  dhb, ms  Got this working again.
%   12/11/14 dhb      Clean up for release.
%   1/6/17   ms       Updated name convention for photopigments.
%   08/10/18 dhb      Add L-M direction demo.  This works for a specified
%                     contrast. I think we have code in
%                     ReceptorIsolateOptimBackgroundMulti that allows
%                     maximizing along a direction. For another day to
%                     remember how that works.
%   10/17/18 dhb      Add validation mode.  

%% Clear and close
close all;

%% Parse input
p = inputParser;
p.addParameter('validate','none',@ischar);
p.addParameter('whichDirectionNumber',1,@isnumeric);
p.addParameter('validationFractionTolerance',0.001,@isnumeric)
p.parse(varargin{:});

%% Set status to OK
status = 1;

%% Run in the directory that contains this file.
curDir = pwd;
cd(fileparts(mfilename('fullpath')));

%% User choice of which of several cases to do
switch (p.Results.validate)
    case 'basichuman'
        whichModel = 'HumanPhotopigments';
        whichPrimaries = 'Spectral';
        basicHuman.maximizeTargetContrast = 1;
        basicHuman.Results.desiredContrast = [];
        
    case {'none'}
        % Interactive mode
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
        
        % Prompt for device to compute with respect to
        fprintf('\nAvailable devices:\n');
        fprintf('\t[1]  Spectral - ideal spectral producing device\n');
        fprintf('\t[2]  OneLight - calibration data from Brainard/Aguirre lab OneLight device\n');
        fprintf('\t[3]  Monitor - some typical monitor\n');
        fprintf('\t[4]  FiveLEDPrimaries - 5 LEDs\n');
        whichPrimaryNumber = GetWithDefault('Enter device',1);
        switch (whichPrimaryNumber)
            case 1
                whichPrimaries = 'Spectral';
            case 2
                whichPrimaries = 'OneLight';
            case 3
                whichPrimaries = 'Monitor';
            case 4
                whichPrimaries = 'FiveLEDPrimaries';
            otherwise
                error('Unknown primaries entered');
        end
    otherwise
        error('Unknown validation type specified');
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
        calPath = fullfile(fileparts(mfilename('fullpath')), 'ReceptorIsolateDemoData', []);
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
    case 'FiveLEDPrimaries'
        % Construct LEDs
        peakWls = [456 488 540 592 632];
        fwhm = [10 10 10 17 17]/2;
        maxPower = [.57 .125 .156 .27 0.75];
        for i = 1:length(fwhm)
            % Figure out the standard deviation.
            standardDeviation(i) = FWHMToStd(fwhm(i));
        end
        
        S = [380 2 201];
        wls = SToWls(S);
        
        % Make the spectrum.
        for i = 1:length(fwhm)
            spd(:, i) = normpdf(wls, peakWls(i), fwhm(i));
            spd(:, i) = spd(:, i)./max(spd(:, i))*maxPower(i);
        end
        
        B_primary = spd;
        
        % Set background to the monitor midpoint, and use the ambient
        % spectrum from the calibration file.
        backgroundPrimary = [0.5 0.5 0.5 0.5 0.5]';
        ambientSpd = zeros(201, 1);
        
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
    case 'HumanPhotopigments'
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
        %
        % Note that we don't typically vary or pass the blood vessel
        % parameters but rather simply accept the defaults used by
        % GetHumanPhotoreceptorSS.  It's mainly for fun that we show
        % how to do this here.
        switch (p.Results.validate)
            case 'basichuman'
                observerAgeInYears = 32;
                fieldSizeDegrees = 27.5;
                pupilDiameterMm = 4.7;
                vesselOxyFraction = 0.85;
                vesselOverallThicknessUm = 5;
                correctBleaching = 0;
            case {'none'}
                % Interactive mode
                observerAgeInYears = GetWithDefault('\tObserver age in years?', 32);
                fieldSizeDegrees = GetWithDefault('\tField size in degrees?', 27.5);
                pupilDiameterMm = GetWithDefault('\tPupil diameter?', 4.7);
                vesselOxyFraction = GetWithDefault(['\tOxygenation fraction for vessel hemoglobin [typical 0.85]?'], 0.85);
                vesselOverallThicknessUm = GetWithDefault(['\tVessel thickness [typical 5 um]?'], 5);
                correctBleaching = GetWithDefault(['\tCorrect for cone photopigment bleaching [1 = yes, 0 = no]?'],0);
            otherwise
                error('Unknown validation type specified');
        end

        % Define photoreceptor classes that we'll consider.
        % ReceptorIsolate has a few more built-ins than these.
        photoreceptorClasses = {'LConeTabulatedAbsorbance', 'MConeTabulatedAbsorbance', 'SConeTabulatedAbsorbance', 'Melanopsin', 'Rods', ...
            'LConeTabulatedAbsorbancePenumbral', 'MConeTabulatedAbsorbancePenumbral', 'SConeTabulatedAbsorbancePenumbral'};
        
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
        % GetHumanPhotoreceptorSS, which is in the ContrastSplatter
        % directory.  Each row of the matrix T_receptors provides the
        % spectral sensitivity of the photoreceptor class in the
        % corresponding entry of the cell array photoreceptorClasses.
        %
        % The last two arguments are the oxygenation fraction and the
        % vessel thickness. We set them to be empty here, prompting the
        % user to enter these values later.
        oxygenationFraction = [];
        vesselThickness = [];
        T_receptors = GetHumanPhotoreceptorSS(S, photoreceptorClasses, fieldSizeDegrees, observerAgeInYears, pupilDiameterMm, [], fractionBleached, oxygenationFraction, vesselThickness);
        
        %% Let user choose a photoreceptor class to target
        LMinusMTargetContrast = 0.06;
        switch (p.Results.validate)
            case 'basichuman'
                whichDirectionNumber = p.Results.whichDirectionNumber;
            case {'none'}
                % Interactive
                fprintf('Available photoreceptor classes to target:\n');
                fprintf('\t[1]  Melanopsin, silence open-field cones; ignore rods and penumbral cones\n');
                fprintf('\t[2]  Melanopsin, silence open-field and penumbral cones; ignore rods)\n');
                fprintf('\t[3]  S cones, silence open-field L and M cones, melanopsin, and prenumbral L and M cones; ignore rods and penumbral S cones\n');
                fprintf('\t[4]  S cones, silence open field L and M cones, ingore all others\n');
                fprintf('\t[5]  Penumbral L and M cones, silence open-field cones, melanopsin, and prenumbral S cones; ignore rods\n');
                fprintf('\t[6]  Rods\n');
                fprintf('\t[7]  L minus M at %0.2f contrast\n',LMinusMTargetContrast);
                whichDirectionNumber = GetWithDefault('Enter direction',1);
            otherwise
                error('Unknown validation type specified');
        end
        
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
                whichDirection = 'SDirected';
                whichReceptorsToTarget = [3];
                whichReceptorsToIgnore = [4 5 6 7 8];
                whichReceptorsToMinimize = [];
            case 5
                whichDirection = 'PenumbralLM';
                whichReceptorsToTarget = [6 7];
                whichReceptorsToIgnore = [5];
                whichReceptorsToMinimize = [];
            case 6
                whichDirection = 'Rods';
                whichReceptorsToTarget = [5];
                whichReceptorsToIgnore = [6 7 8];
                whichReceptorsToMinimize = [];
            case 7
                % Here we want to obtain desired contrasts on two
                % photoreceptor classes.  We don't have code in place to
                % maximize in a specified direction, but we can explicitly
                % specify the desired contrasts on the targeted directions.
                % One can then maximize by hand, finding out how far one
                % can go before it is out of gamut.  A little loop and
                % check could automate this process.
                whichDirection = 'LMinusM';
                whichReceptorsToTarget = [1 2];
                whichReceptorsToIgnore = [4 5 6 7 8];
                whichReceptorsToMinimize = [];
                desiredContrast = [LMinusMTargetContrast -LMinusMTargetContrast];
            otherwise
                error('Unknown direction entered');
        end
        
        
    case 'Dog'
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
%
% If we target, here we specify the same contrast for all targeted classes.
% This is not necessary, they can differ.  It just makes the demo code a
% bit simpler to yoke them since we only have to prompt for one number.
if (~exist('desiredContrast','var') | isempty(desiredContrast))
    switch (p.Results.validate)
        case 'basichuman'
            maximizeTargetContrast = basicHuman.maximizeTargetContrast;
            if maximizeTargetContrast
                desiredContrast = [];
            else
                desiredContrast = basicHuman.desiredContrast;
            end
            
        case {'none'}
            % Interactive mode
            maximizeTargetContrast = GetWithDefault('\tMaximize contrast? [1 = yes, 0 = no]', 1); 
            if maximizeTargetContrast
                desiredContrast = [];
            else
                desiredContrast = GetWithDefault('\tDesired contrast (applies to all targeted classes)?', 0.45)*ones(size(whichReceptorsToTarget));
            end
        otherwise
            error('Unknown validation type specified');
    end   
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

% Positive modulation of receptors
modulationReceptors = T_receptors*B_primary*(modulationPrimary - backgroundPrimary);
contrastReceptors = modulationReceptors ./ backgroundReceptors;
fprintf('Positive modulation contrasts\n');
for j = 1:size(T_receptors,1)
    fprintf('\t%s: contrast = %0.4f\n',photoreceptorClasses{j},contrastReceptors(j));
end

% Negative modulation of receptors
modulationReceptors = T_receptors*B_primary*(-(modulationPrimary - backgroundPrimary));
contrastReceptors = modulationReceptors ./ backgroundReceptors;
fprintf('Negative modulation contrasts\n');
for j = 1:size(T_receptors,1)
    fprintf('\t%s: contrast = %0.4f\n',photoreceptorClasses{j},contrastReceptors(j));
end

%% Plots and/or validation check
switch (p.Results.validate)
    case 'none'
        % Make lots in interactive mode
        plotDir = 'ReceptorIsolateDemoOutput';
        if ~isdir(plotDir)
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
        saveas(theFig1,sprintf('%s_%s_%s_Sensitivities.pdf',whichModel,whichPrimaries,whichDirection),'pdf');
        
        % Modulation spectra
        theFig2 = figure; hold on
        plot(SToWls(S),B_primary*modulationPrimary,'r','LineWidth',2);
        plot(SToWls(S),B_primary*backgroundPrimary,'k','LineWidth',2);
        title('Modulation spectra');
        xlim([380 780]);
        xlabel('Wavelength');
        ylabel('Power');
        pbaspect([1 1 1]);
        saveas(theFig2,sprintf('%s_%s_%s_Modulation.pdf',whichModel,whichPrimaries,whichDirection),'pdf');
        
        % Primaries
        theFig3 = figure; hold on
        plot(modulationPrimary,'r','LineWidth',2);
        plot(backgroundPrimary+(-(modulationPrimary-backgroundPrimary)),'g','LineWidth',2);
        plot(backgroundPrimary,'k','LineWidth',2);
        title('Primary settings');
        xlim([0 length(backgroundPrimary)]);
        ylim([0 1]);
        xlabel('Primary Number (nominal)');
        ylabel('Setting');
        legend({'Positive', 'Negative', 'Background'},'Location','NorthEastOutside');
        saveas(theFig3,sprintf('%s_%s_%s_Primaries.pdf',whichModel,whichPrimaries,whichDirection),'pdf');
    case 'basichuman'
        % Check that some numbers computed now match what they were when we
        % set this up.
        switch p.Results.whichDirectionNumber
            case 1
                receptorsCheck = sum(T_receptors(:));
                if (abs(max(receptorsCheck - 189.39))/189.39 > p.Results.validationFractionTolerance)
                    status = 0;
                end
                modulationPrimaryCheck = sum(modulationPrimary(:));
                if (abs(max(modulationPrimaryCheck - 47.553))/47.553 > p.Results.validationFractionTolerance)
                    status = 0;
                end
        end
    otherwise
        error('Unknown value for validate type passed.')
        
end

%% Return to the directory from whence we started
cd(curDir);