function [T_energyNormalized,T_quantalIsomerizations,nominalLambdaMax] = GetHumanPhotoreceptorSS(S, photoreceptorClasses, fieldSizeDegrees, ageInYears,...
    pupilDiameterMm, lambdaMaxShift, fractionPigmentBleached, vesselOxyFraction, vesselOverallThicknessUm)
% [T_energyNormalized,T_quantalIsomerizations,nominalLambdaMax] = GetHumanPhotoreceptorSS(S, photoreceptorClasses, fieldSizeDegrees, ageInYears, pupilDiameterMm,
%   lambdaMaxShift, fractionPigmentBleached, vesselOxyFraction, vesselOverallThicknessUm)
%
% Produces photopigment sensitivities that we often need, and allowing
% variation in age and lambda-max.  T_energyNormalized are the sensitivities
% in energy units, normalized to max of 1.  T_quantalIosmerizations are the
% probability of an isomerization for quantal unit input.
%
% If empty variables are passed for any of the following variables, defaults will be assumed.
%
% Input:
%   S (1x3)                         - Wavelength spacing.
%                                     Default: [380 2 201]
%   photoreceptorClasses (cell)     - Cell with names of photoreceptor classes.
%                                     Supported options: 'LCone', 'MCone', 'SCone', 'Melanopsin', 'Rods', ...
%                                                        'LCone10DegTabulatedSS', 'MCone10DegTabulatedSS', 'SCone10DegTabulatedSS', ...
%                                                        'MelanopsinLegacy', 'RodsLegacy', 'CIE1924VLambda' ...
%                                                        'LConeHemo', 'MConeHemo', 'SConeHemo'
%                                     Default: {'LCone' ; 'MCone' ; 'SCone'}
%   fieldSizeDegrees (1x1)          - Field size in degrees.
%                                     Default: 10
%   ageInYears (1x1)                - Observer age in years.
%                                     Default: 32
%   pupilDiameterMm (1x1)           - Pupil diameter in mm.
%                                     Default: 3
%   lambdaMaxShift                  - Shift of lambda-max.
%                                     Can be scalar, in which case same value is used for all classes.  Or a
%                                     vector of same length as photoreceptorClasses.
%                                     Default: 0
%   fractionPigmentBleached         - Fraction of pigment bleached.
%                                     Vector of same size as number of photoreceptor classes.
%                                     Default: 0 for each entry.
%   vesselOxyFraction               - Fraction of oxygenated blood assumed.
%                                     Default: 0.85
%   vesselOverallThicknessUm        - Thickness of vessel in um
%                                     Default: 5
% Output:
%   T_energyNormalized              - Spectral sensitivities in energy
%                                     units (normalized to max.).
%   T_quantaIsomerization           - Spectral sensitivities in quanta
%                                     units.  These may be used to compute
%                                     isomerizations from retinal illuminance
%   nominalLambdaMax                - Peak of photopigment absorption spectrum spectral sensitivities (no density taken into account).
%
% NOTES:
%   A) The Hemo variants take into account an estimate of the absorption spectrum of hemoglobin as seen through retinal
%   blood vessels.
%   B) Not all variants have a meaningful T_quantalIsomerizations variable returned.  When we don't have that sensitivity
%   easily, a vector of NaN's of the right size is returned instead.
%   C) The 'Legacy' versions are estimates that we used in some of our early work and that we now don't think are as good as the current ones.
%   They are still here because sometimes we need to figure out what we did in the past.
%   D) Passed fractionBleached numbers don't affect rod, melanopsin, or tabulated cone calculations.
%   We throw an error if either of these ever receives a non-zero value.
%   E) If you pass lambdaMaxShift as a scalar, the same value of lambdaMaxShift is applied to all classes computed within a single call to this function.
%   This seems unlikely to be a good behavior, but I am keeping it for backwards compatibility.  It is OK if the scalar contains If you want to shift each
%   separately, use a loop around this function and pass a different value for each call.  We could generalize this code to take a vector of
%   lambda max shifts in, and that might be smart to do.
%
% 1/21/14   ms    Wrote it based on old code.
% 5/24/14   dhb   Remove vestigal references to a returned labels variable.
% 5/26/14   dhb   Fix bug: 'Melanopsin-2' was being computed with a shift of -1.
%           dhb   Simplify return interface.  Add many comments.
%           dhb   Return isomerization sensitivities for hemoglobin variants.
% 11/21/14  ms    Cleaned up and commented
% 12/12/14  dhb   More cleaning and comments.
%                 This now takes hemoglobin parameters.
%                 Manuel did this but didn't leave a comment.  I pushed the
%                 responsibility for prompting up to the caller.
% 10/22/15  dhb   Added note about how the lambdaMaxShift parameter works.

%% Set defaults

% Wavelength sampling
if (nargin < 1 | isempty(S))
    S = [380 2 201];
end

% Photoreceptor classes to generate
if (nargin < 2 | isempty(photoreceptorClasses))
    photoreceptorClasses = {'LCone' ; 'MCone' ; 'SCone'};
end

% Field size
if (nargin < 3 | isempty(fieldSizeDegrees))
    fieldSizeDegrees = 10;
end

% Observer age
% If the passed observer age is <20 or >80, we assume that the observer is
% 20 or 80 respectively, which are the maximum ages given by the CIE standard.
if (nargin < 4 | isempty(ageInYears))
    ageInYears = 32;
end
if ageInYears < 20
    ageInYears = 20;
    %fprintf('Observer age truncated at 20\n');
end
if ageInYears > 80
    ageInYears = 80;
    %fprintf('Observer age truncated at 80\n');
end

% Pupil diameter
if (nargin < 5 | isempty(pupilDiameterMm))
    pupilDiameterMm = 3;
end

% Shift of pigment lambda max from nominal value, and some
% sanity checks on what we get.
if (nargin < 6 | isempty(lambdaMaxShift))
    lambdaMaxShift = 0;
end
if (length(lambdaMaxShift) == 1 & length(photoreceptorClasses) ~= 1 & lambdaMaxShift ~= 0)
    error('A scalar but non-zero lambdaMaxShift was passed.  This will not lead to good things.  Fix it.')
end
if (length(lambdaMaxShift) ~= 1 & length(lambdaMaxShift) ~= length(photoreceptorClasses))
    error('lambdaMaxShift passed as a vector with length not equally to number of photoreceptor classes.');
end

% Fraction pigment bleached
if (nargin < 7 | isempty(fractionPigmentBleached))
    if length(photoreceptorClasses) > 1
        fractionPigmentBleached = zeros(length(photoreceptorClasses),1);
    else
        fractionPigmentBleached = 0;
    end
end

% Vessel oxygenation
if (nargin < 8 | isempty(vesselOxyFraction))
    vesselOxyFraction = 0.85;
end

% Vessel thickness
if (nargin < 9 | isempty(vesselOverallThicknessUm))
    vesselOverallThicknessUm = 5;
end

%% Assign empty vectors
T_quanta = [];
T_energyNormalized = [];
T_quantalIsomerizations = [];
nominalLambdaMax = [];

%% Fussing for bleaching calcs
%
% The fractionPigmentBleached vectors come in the same dimensions as
% photoreceptors. However, ComputeCIEConeFundamentals expects LMS triplets.
% So, we sort out the hemo vs. non-hemo fractions. We do that because we do
% not know the order of photopigment classes passed into this function. In
% case we only have one cone type passed, which sometimes happens, we
% still extract the triplet of fractions of pigment bleached, under the
% assumption that in the input vector, the order is LMS. This is a bit
% kludge-y, but works.
%
% I don't think this is quite right for cases where multiple "R" variants
% of a class get passed (e.g. both LCone and Lcone2Deg are passed in one call,
% becuase later values in the passed vector clobber earlier ones in the vectors that get set here.
% This needs some thought.
if length(photoreceptorClasses) > 1
    for i = 1:length(photoreceptorClasses)
        switch photoreceptorClasses{i}
            case {'LCone'}
                fractionConeBleachedFromIsom(1) = fractionPigmentBleached(i);
            case {'MCone'}
                fractionConeBleachedFromIsom(2) = fractionPigmentBleached(i);
            case {'SCone'}
                fractionConeBleachedFromIsom(3) = fractionPigmentBleached(i);
            case {'LCone2Deg'}
                fractionConeBleachedFromIsom(1) = fractionPigmentBleached(i);
            case {'MCone2Deg'}
                fractionConeBleachedFromIsom(2) = fractionPigmentBleached(i);
            case {'SCone2Deg'}
                fractionConeBleachedFromIsom(3) = fractionPigmentBleached(i);
            case {'LConeHemo'}
                fractionConeBleachedFromIsomHemo(1) = fractionPigmentBleached(i);
            case {'MConeHemo'}
                fractionConeBleachedFromIsomHemo(2) = fractionPigmentBleached(i);
            case {'SConeHemo'}
                fractionConeBleachedFromIsomHemo(3) = fractionPigmentBleached(i);
            case {'LConeHemoLegacy'}
                fractionConeBleachedFromIsomHemo(1) = fractionPigmentBleached(i);
            case {'MConeHemoLegacy'}
                fractionConeBleachedFromIsomHemo(2) = fractionPigmentBleached(i);
            case {'SConeHemoLegacy'}
                fractionConeBleachedFromIsomHemo(3) = fractionPigmentBleached(i);
            case {'Melanopsin', 'MelanopsinLegacy'}
                if (fractionPigmentBleached(i) ~= 0)
                    error('Non-zero fractionPigmentBleached passed for photoreceptor class that does not support this');
                end
            case {'Rods', 'RodsLegacy', 'CIE1924VLambda'}
                if (fractionPigmentBleached(i) ~= 0)
                    error('Non-zero fractionPigmentBleached passed for photoreceptor class that does not support this');
                end
            case {'LCone10DegTabulatedSS', 'MCone10DegTabulatedSS', 'SCone10DegTabulatedSS'}
                if (fractionPigmentBleached(i) ~= 0)
                    error('Non-zero fractionPigmentBleached passed for photoreceptor class that does not support this');
                end
            otherwise
                error('Unknown photoreceptor class passed');
        end
    end
    
% If only one cone class is passed, which can happen in splatter
% calculations, we set the fraction pigment bleached for the pigments that
% are not passed to be 0. This is because PTB machinery expects triplets.
elseif length(photoreceptorClasses) == 1
    switch photoreceptorClasses{1}
        case 'LCone'
            fractionConeBleachedFromIsom(1) = fractionPigmentBleached;
            fractionConeBleachedFromIsom(2) = 0;
            fractionConeBleachedFromIsom(3) = 0;
        case 'MCone'
            fractionConeBleachedFromIsom(1) = 0;
            fractionConeBleachedFromIsom(2) = fractionPigmentBleached;
            fractionConeBleachedFromIsom(3) = 0;
        case 'SCone'
            fractionConeBleachedFromIsom(1) = 0;
            fractionConeBleachedFromIsom(2) = 0;
            fractionConeBleachedFromIsom(3) = fractionPigmentBleached;
        case 'LCone2Deg'
            fractionConeBleachedFromIsom(1) = fractionPigmentBleached;
            fractionConeBleachedFromIsom(2) = 0;
            fractionConeBleachedFromIsom(3) = 0;
        case 'MCone2Deg'
            fractionConeBleachedFromIsom(1) = 0;
            fractionConeBleachedFromIsom(2) = fractionPigmentBleached;
            fractionConeBleachedFromIsom(3) = 0;
        case 'SCone2Deg'
            fractionConeBleachedFromIsom(1) = 0;
            fractionConeBleachedFromIsom(2) = 0;
            fractionConeBleachedFromIsom(3) = fractionPigmentBleached;
        case 'LConeHemo'
            fractionConeBleachedFromIsomHemo(1) = fractionPigmentBleached;
            fractionConeBleachedFromIsomHemo(2) = 0;
            fractionConeBleachedFromIsomHemo(3) = 0;
        case 'MConeHemo'
            fractionConeBleachedFromIsomHemo(1) = 0;
            fractionConeBleachedFromIsomHemo(2) = fractionPigmentBleached;
            fractionConeBleachedFromIsomHemo(3) = 0;
        case 'SConeHemo'
            fractionConeBleachedFromIsomHemo(1) = 0;
            fractionConeBleachedFromIsomHemo(2) = 0;
            fractionConeBleachedFromIsomHemo(3) = fractionPigmentBleached;
        case 'LConeHemoLegacy'
            fractionConeBleachedFromIsomHemo(1) = fractionPigmentBleached;
            fractionConeBleachedFromIsomHemo(2) = 0;
            fractionConeBleachedFromIsomHemo(3) = 0;
        case 'MConeHemoLegacy'
            fractionConeBleachedFromIsomHemo(1) = 0;
            fractionConeBleachedFromIsomHemo(2) = fractionPigmentBleached;
            fractionConeBleachedFromIsomHemo(3) = 0;
        case 'SConeHemoLegacy'
            fractionConeBleachedFromIsomHemo(1) = 0;
            fractionConeBleachedFromIsomHemo(2) = 0;
            fractionConeBleachedFromIsomHemo(3) = fractionPigmentBleached;
    end
end

% Transpose if we can. We do this because the PTB machinery expects this.
if exist('fractionConeBleachedFromIsom', 'var')
    fractionConeBleachedFromIsom = fractionConeBleachedFromIsom';
end

if exist('fractionConeBleachedFromIsomHemo', 'var')
    fractionConeBleachedFromIsomHemo = fractionConeBleachedFromIsomHemo';
end

%% Iterate over the photoreceptor classes that have been passed.
for i = 1:length(photoreceptorClasses)
    theClass = photoreceptorClasses{i};
    
    % Get lambdaMaxShift to use for this class.
    if (length(lambdaMaxShift) == 1)
        lambdaMaxShiftUse = lambdaMaxShift;
    elseif (length(lambdaMaxShift) == length(photoreceptorClasses))
        lambdaMaxShiftUse = lambdaMaxShift(theClass);
    else
        error('Input lambdaMaxShift does not have an allowable dimension.');
    end
    
    switch theClass
        case 'LCone'
            whichNomogram = 'StockmanSharpe';
            lambdaMax = [558.9 530.3 420.7]';
            
            %% Construct cones, pull out L cone
            [T_quantalNormalized1,~,T_quantalIsomerizations1] = ComputeCIEConeFundamentals(S,fieldSizeDegrees,ageInYears,pupilDiameterMm,lambdaMax+lambdaMaxShiftUse,whichNomogram,[],[],[],fractionConeBleachedFromIsom);
            T_energy1 = EnergyToQuanta(S,T_quantalNormalized1')';
            
            % Add to the receptor vector
            T_energyNormalized = [T_energyNormalized ; T_energy1(1,:)];
            T_quantalIsomerizations = [T_quantalIsomerizations ; T_quantalIsomerizations1(1,:)];
            nominalLambdaMax = [nominalLambdaMax lambdaMax(1)+lambdaMaxShiftUse];
        case 'MCone'
            whichNomogram = 'StockmanSharpe';
            lambdaMax = [558.9 530.3 420.7]';
            
            %% Construct cones, pull out M cone
            [T_quantalNormalized1,~,T_quantalIsomerizations1] = ComputeCIEConeFundamentals(S,fieldSizeDegrees,ageInYears,pupilDiameterMm,lambdaMax+lambdaMaxShiftUse,whichNomogram,[],[],[],fractionConeBleachedFromIsom);
            T_energy1 = EnergyToQuanta(S,T_quantalNormalized1')';
            
            % Add to the receptor vector
            T_energyNormalized = [T_energyNormalized ; T_energy1(2,:)];
            T_quantalIsomerizations = [T_quantalIsomerizations ; T_quantalIsomerizations1(2,:)];
            nominalLambdaMax = [nominalLambdaMax lambdaMax(2)+lambdaMaxShiftUse];
        case 'SCone'
            whichNomogram = 'StockmanSharpe';
            lambdaMax = [558.9 530.3 420.7];
            
            %% Construct cones, pull out S cone
            [T_quantalNormalized1,~,T_quantalIsomerizations1] = ComputeCIEConeFundamentals(S,fieldSizeDegrees,ageInYears,pupilDiameterMm,lambdaMax+lambdaMaxShiftUse,whichNomogram,[],[],[],fractionConeBleachedFromIsom);
            T_energy1 = EnergyToQuanta(S,T_quantalNormalized1')';
            
            % Add to the receptor vector
            T_energyNormalized = [T_energyNormalized ; T_energy1(3,:)];
            T_quantalIsomerizations = [T_quantalIsomerizations ; T_quantalIsomerizations1(3,:)];
            nominalLambdaMax = [nominalLambdaMax lambdaMax(3)+lambdaMaxShiftUse];
        case 'LCone2Deg'
            whichNomogram = 'StockmanSharpe';
            lambdaMax = [558.9 530.3 420.7]';
            
            %% Construct cones, pull out L cone
            [T_quantalNormalized1,~,T_quantalIsomerizations1] = ComputeCIEConeFundamentals(S,2,ageInYears,pupilDiameterMm,lambdaMax+lambdaMaxShiftUse,whichNomogram,[],[],[],fractionConeBleachedFromIsom);
            T_energy1 = EnergyToQuanta(S,T_quantalNormalized1')';
            
            % Add to the receptor vector
            T_energyNormalized = [T_energyNormalized ; T_energy1(1,:)];
            T_quantalIsomerizations = [T_quantalIsomerizations ; T_quantalIsomerizations1(1,:)];
            nominalLambdaMax = [nominalLambdaMax lambdaMax(1)+lambdaMaxShiftUse];
        case 'MCone2Deg'
            whichNomogram = 'StockmanSharpe';
            lambdaMax = [558.9 530.3 420.7]';
            
            %% Construct cones, pull out M cone
            [T_quantalNormalized1,~,T_quantalIsomerizations1] = ComputeCIEConeFundamentals(S,2,ageInYears,pupilDiameterMm,lambdaMax+lambdaMaxShiftUse,whichNomogram,[],[],[],fractionConeBleachedFromIsom);
            T_energy1 = EnergyToQuanta(S,T_quantalNormalized1')';
            
            % Add to the receptor vector
            T_energyNormalized = [T_energyNormalized ; T_energy1(2,:)];
            T_quantalIsomerizations = [T_quantalIsomerizations ; T_quantalIsomerizations1(2,:)];
            nominalLambdaMax = [nominalLambdaMax lambdaMax(2)+lambdaMaxShiftUse];
        case 'SCone2Deg'
            whichNomogram = 'StockmanSharpe';
            lambdaMax = [558.9 530.3 420.7];
            
            %% Construct cones, pull out S cone
            [T_quantalNormalized1,~,T_quantalIsomerizations1] = ComputeCIEConeFundamentals(S,2,ageInYears,pupilDiameterMm,lambdaMax+lambdaMaxShiftUse,whichNomogram,[],[],[],fractionConeBleachedFromIsom);
            T_energy1 = EnergyToQuanta(S,T_quantalNormalized1')';
            
            % Add to the receptor vector
            T_energyNormalized = [T_energyNormalized ; T_energy1(3,:)];
            T_quantalIsomerizations = [T_quantalIsomerizations ; T_quantalIsomerizations1(3,:)];
            nominalLambdaMax = [nominalLambdaMax lambdaMax(3)+lambdaMaxShiftUse];
        case {'Melanopsin'}
            % Melanopsin
            photoreceptors = DefaultPhotoreceptors('LivingHumanMelanopsin');
            photoreceptors.nomogram.S = S;
            nominalLambdaMaxTmp = photoreceptors.nomogram.lambdaMax+lambdaMaxShiftUse;
            photoreceptors.nomogram.lambdaMax = photoreceptors.nomogram.lambdaMax+lambdaMaxShiftUse;
            photoreceptors.fieldSizeDegrees = fieldSizeDegrees;
            photoreceptors.ageInYears = ageInYears;
            photoreceptors.pupilDiameter.value = pupilDiameterMm;
            photoreceptors = FillInPhotoreceptors(photoreceptors);
            
            % Add to the receptor vector
            T_energyNormalized = [T_energyNormalized ; photoreceptors.energyFundamentals];
            T_quantalIsomerizations = [T_quantalIsomerizations ; photoreceptors.isomerizationAbsorptance];
            nominalLambdaMax = [nominalLambdaMax nominalLambdaMaxTmp];
        case 'Rods'
            % Rods
            photoreceptors = DefaultPhotoreceptors('LivingHumanRod');
            photoreceptors.nomogram.S = S;
            nominalLambdaMaxTmp = photoreceptors.nomogram.lambdaMax+lambdaMaxShiftUse;
            photoreceptors.nomogram.lambdaMax = photoreceptors.nomogram.lambdaMax+lambdaMaxShiftUse;
            photoreceptors.fieldSizeDegrees = fieldSizeDegrees;
            photoreceptors.ageInYears = ageInYears;
            photoreceptors.pupilDiameter.value = pupilDiameterMm;
            photoreceptors = FillInPhotoreceptors(photoreceptors);
            
            % Add to the receptor vector
            T_energyNormalized = [T_energyNormalized ; photoreceptors.energyFundamentals];
            T_quantalIsomerizations = [T_quantalIsomerizations ; photoreceptors.isomerizationAbsorptance];
            nominalLambdaMax = [nominalLambdaMax nominalLambdaMaxTmp];
         
        %%%%% NEED TO FIX NOMINAL LAMBDA MAX HERE AND BELOW.  %%%%%
        case 'LCone10DegTabulatedSS'
            % Load in the tabulated 10-deg S-S fundamentals
            targetRaw = load('T_cones_ss10');
            T_energy_tmp = SplineCmf(targetRaw.S_cones_ss10,targetRaw.T_cones_ss10(1,:),S,2);
            T_energyNormalized = [T_energyNormalized ; T_energy_tmp];
            T_quanta = [T_quanta ; QuantaToEnergy(S,T_energy_tmp')'];
            T_quantalIsomerizations = [T_quantalIsomerizations ; NaN*ones(size(T_quanta))];
        case 'MCone10DegTabulatedSS'
            % Load in the tabulated 10-deg S-S fundamentals
            targetRaw = load('T_cones_ss10');
            T_energy_tmp = SplineCmf(targetRaw.S_cones_ss10,targetRaw.T_cones_ss10(2,:),S,2);
            T_energyNormalized = [T_energyNormalized ; T_energy_tmp];
            T_quanta = [T_quanta ; QuantaToEnergy(S,T_energy_tmp')'];
            T_quantalIsomerizations = [T_quantalIsomerizations ; NaN*ones(size(T_quanta))];
        case 'SCone10DegTabulatedSS'
            % Load in the tabulated 10-deg S-S fundamentals
            targetRaw = load('T_cones_ss10');
            T_energy_tmp = SplineCmf(targetRaw.S_cones_ss10,targetRaw.T_cones_ss10(3,:),S,2);
            T_energyNormalized = [T_energyNormalized ; T_energy_tmp];
            T_quantalIsomerizations = [T_quantalIsomerizations ; NaN*ones(size(T_quanta))];
        case 'MelanopsinLegacy'
            % Construct the melanopsin receptor
            whichNomogram = 'StockmanSharpe';
            lambdaMax = [558.9, 530.3, 480+lambdaMaxShiftUse];
            
            % Make a call to ComputeCIEConeFundamentals() which makes appropriate calls
            T_quanta_tmp = ComputeCIEConeFundamentals(S,10,ageInYears,3,lambdaMax,whichNomogram);
            T_energyNormalized = [T_energyNormalized ; EnergyToQuanta(S,T_quanta_tmp(3,:)')'];
            T_quantalIsomerizations = [T_quantalIsomerizations ; NaN*ones(size(T_quanta))];
            nominalLambdaMax = [nominalLambdaMax 480];
        case 'RodsLegacy'
            whichNomogram = 'StockmanSharpe';
            lambdaMax = [558.9, 530.3, 480];
            ageInYears = 32;
            pupilSize = 3;
            fieldSize = 10;
            lambdaMaxRods = 500;
            
            DORODS = true;
            T_quanta_tmp = ComputeCIEConeFundamentals(S,fieldSize,ageInYears,pupilSize,lambdaMaxRods,whichNomogram,[],DORODS);
            T_energyNormalized = [T_energyNormalized ; EnergyToQuanta(S,T_quanta_tmp')'];
            T_quantalIsomerizations = [T_quantalIsomerizations ; NaN*ones(size(T_quanta))];
        case 'CIE1924VLambda'
            
            % Load in the CIE 1959 scotopic luminosity function
            targetRaw = load('T_rods');
            T_energyNormalized = [T_energyNormalized ; SplineCmf(targetRaw.S_rods,targetRaw.T_rods,S,2)];
            T_quantalIsomerizations = [T_quantalIsomerizations ; NaN*ones(size(T_quanta))];
        case 'LConeHemo'
            whichNomogram = 'StockmanSharpe';
            lambdaMax = [558.9 530.3 420.7]';
            
            %% Construct cones, pull out L cone
            [T_quantalNormalized1,~,T_quantalIsomerizations1] = ComputeCIEConeFundamentals(S,fieldSizeDegrees,ageInYears,pupilDiameterMm,lambdaMax+lambdaMaxShiftUse,whichNomogram,[],[],[],fractionConeBleachedFromIsomHemo);
            T_energy1 = EnergyToQuanta(S,T_quantalNormalized1')';
            
            % Multiply with blood transmittance, which is parametrized by
            % the oxygenation fraction and vessel thickness.
            source = 'Prahl';
            trans_Hemoglobin = GetHemoglobinTransmittance(S,vesselOxyFraction,vesselOverallThicknessUm,source);
            
            % Add to the receptor vector
            T_energyNormalized = [T_energyNormalized ; T_energy1(1,:) .* trans_Hemoglobin'];
            T_quantalIsomerizations = [T_quantalIsomerizations ; T_quantalIsomerizations1(1,:) .* trans_Hemoglobin'];
            nominalLambdaMax = [nominalLambdaMax lambdaMax(1)];
        case 'MConeHemo'
            whichNomogram = 'StockmanSharpe';
            lambdaMax = [558.9 530.3 420.7]';
            
            %% Construct cones, pull out M cone
            [T_quantalNormalized1,~,T_quantalIsomerizations1] = ComputeCIEConeFundamentals(S,fieldSizeDegrees,ageInYears,pupilDiameterMm,lambdaMax+lambdaMaxShiftUse,whichNomogram,[],[],[],fractionConeBleachedFromIsomHemo);
            T_energy1 = EnergyToQuanta(S,T_quantalNormalized1')';
            
            % Multiply with blood transmittance
            source = 'Prahl';
            trans_Hemoglobin = GetHemoglobinTransmittance(S,vesselOxyFraction,vesselOverallThicknessUm,source);
            
            % Add to the receptor vector
            T_energyNormalized = [T_energyNormalized ; T_energy1(2,:) .* trans_Hemoglobin'];
            T_quantalIsomerizations = [T_quantalIsomerizations ; T_quantalIsomerizations1(2,:) .* trans_Hemoglobin'];
            nominalLambdaMax = [nominalLambdaMax lambdaMax(2)];
        case 'SConeHemo'
            whichNomogram = 'StockmanSharpe';
            lambdaMax = [558.9 530.3 420.7];
            
            %% Construct cones, pull out S cone
            [T_quantalNormalized1,~,T_quantalIsomerizations1] = ComputeCIEConeFundamentals(S,fieldSizeDegrees,ageInYears,pupilDiameterMm,lambdaMax+lambdaMaxShiftUse,whichNomogram,[],[],[],fractionConeBleachedFromIsomHemo);
            T_energy1 = EnergyToQuanta(S,T_quantalNormalized1')';
            
            % Multiply with blood transmittance
            source = 'Prahl';
            trans_Hemoglobin = GetHemoglobinTransmittance(S,vesselOxyFraction,vesselOverallThicknessUm,source);
            
            % Add to the receptor vector
            T_energyNormalized = [T_energyNormalized ; T_energy1(3,:).* trans_Hemoglobin'];
            T_quantalIsomerizations = [T_quantalIsomerizations ; T_quantalIsomerizations1(3,:) .* trans_Hemoglobin'];
            nominalLambdaMax = [nominalLambdaMax lambdaMax(3)];
        case 'LConeHemoLegacy'
            whichNomogram = 'StockmanSharpe';
            lambdaMax = [558.9 530.3 420.7]';
            
            %% Construct cones, pull out L cone
            [T_quantalNormalized1,~,T_quantalIsomerizations1] = ComputeCIEConeFundamentals(S,fieldSizeDegrees,ageInYears,pupilDiameterMm,lambdaMax+lambdaMaxShiftUse,whichNomogram,[],[],[],fractionConeBleachedFromIsomHemo);
            T_energy1 = EnergyToQuanta(S,T_quantalNormalized1')';
            
            % Multiply with blood transmittance, which is parametrized by
            % the oxygenation fraction and vessel thickness.
            source = 'Prahl';
            trans_Hemoglobin = GetHemoglobinTransmittance(S,vesselOxyFraction,vesselOverallThicknessUm*2.303,source);
            
            % Add to the receptor vector
            T_energyNormalized = [T_energyNormalized ; T_energy1(1,:) .* trans_Hemoglobin'];
            T_quantalIsomerizations = [T_quantalIsomerizations ; T_quantalIsomerizations1(1,:) .* trans_Hemoglobin'];
            nominalLambdaMax = [nominalLambdaMax lambdaMax(1)];
        case 'MConeHemoLegacy'
            whichNomogram = 'StockmanSharpe';
            lambdaMax = [558.9 530.3 420.7]';
            
            %% Construct cones, pull out M cone
            [T_quantalNormalized1,~,T_quantalIsomerizations1] = ComputeCIEConeFundamentals(S,fieldSizeDegrees,ageInYears,pupilDiameterMm,lambdaMax+lambdaMaxShiftUse,whichNomogram,[],[],[],fractionConeBleachedFromIsomHemo);
            T_energy1 = EnergyToQuanta(S,T_quantalNormalized1')';
            
            % Multiply with blood transmittance
            source = 'Prahl';
            trans_Hemoglobin = GetHemoglobinTransmittance(S,vesselOxyFraction,vesselOverallThicknessUm*2.303,source);
            
            % Add to the receptor vector
            T_energyNormalized = [T_energyNormalized ; T_energy1(2,:) .* trans_Hemoglobin'];
            T_quantalIsomerizations = [T_quantalIsomerizations ; T_quantalIsomerizations1(2,:) .* trans_Hemoglobin'];
            nominalLambdaMax = [nominalLambdaMax lambdaMax(2)];
        case 'SConeHemoLegacy'
            whichNomogram = 'StockmanSharpe';
            lambdaMax = [558.9 530.3 420.7];
            
            %% Construct cones, pull out S cone
            [T_quantalNormalized1,~,T_quantalIsomerizations1] = ComputeCIEConeFundamentals(S,fieldSizeDegrees,ageInYears,pupilDiameterMm,lambdaMax+lambdaMaxShiftUse,whichNomogram,[],[],[],fractionConeBleachedFromIsomHemo);
            T_energy1 = EnergyToQuanta(S,T_quantalNormalized1')';
            
            % Multiply with blood transmittance
            source = 'Prahl';
            trans_Hemoglobin = GetHemoglobinTransmittance(S,vesselOxyFraction,vesselOverallThicknessUm*2.303,source);
            
            % Add to the receptor vector
            T_energyNormalized = [T_energyNormalized ; T_energy1(3,:).* trans_Hemoglobin'];
            T_quantalIsomerizations = [T_quantalIsomerizations ; T_quantalIsomerizations1(3,:) .* trans_Hemoglobin'];
            nominalLambdaMax = [nominalLambdaMax lambdaMax(3)];
    end
end

%% Normalize energy sensitivities.
%
% They might already be normalized in most cases, but this makes sure.
for i = 1:size(T_energyNormalized)
    T_energyNormalized(i,:) = T_energyNormalized(i,:)/max(T_energyNormalized(i,:));
end

