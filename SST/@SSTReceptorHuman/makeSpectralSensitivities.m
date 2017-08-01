function obj = makeSpectralSensitivities(obj)
% obj = makeSpectralSensitivities(obj)
%
% This method of @SSTReceptorHuman creates the point estimate of the cone
% fundamentals using machinery from Psychtoolbox-3. By default, both the
% LMS cone fundamentals, and the melanopsin and rod spectral sensitivities
% are returned.
%
% The outputs are returned to the field "T" of the receptor object, in the
% following formats:
%     T.T_quantalIsomerizations - Quantal isomerizations
%     T.T_quantalAbsorptions - Quantal absorptions
%     T.T_quantalAbsorptionsNormalized - Normalized quantal absoprtions
%     T.T_energy - Energy fundamentals
%     T.T_energyNormalized - Normalized energy fundamentals%
%
% 7/25/17   ms  Commented.

%% LMS cones
[T_quantalAbsorptionsNormalizedLMS,T_quantalAbsorptionsLMS,T_quantalIsomerizationsLMS] = ComputeCIEConeFundamentals(obj.S,...
    obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm);

%% Melanopsin
[T_quantalAbsorptionsNormalizedMel,T_quantalAbsorptionsMel,T_quantalIsomerizationsMel] = ComputeCIEMelFundamental(obj.S,...
    obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,[]);

%% Rod
[T_quantalAbsorptionsNormalizedRod,T_quantalAbsorptionsRod,T_quantalIsomerizationsRod] = ComputeCIERodFundamental(obj.S,...
    obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,[]);

%% Assemble the sensitivities
% Normalized quantal sensitivities
T_quantalAbsorptionsNormalized = [T_quantalAbsorptionsNormalizedLMS ; T_quantalAbsorptionsNormalizedMel ; T_quantalAbsorptionsNormalizedRod];

% Quantal isomerizations
T_quantalIsomerizations = [T_quantalIsomerizationsLMS ; T_quantalIsomerizationsMel ; T_quantalIsomerizationsRod];

% Quantal absorption
T_quantalAbsorptions = [T_quantalAbsorptionsLMS ; T_quantalAbsorptionsMel ; T_quantalAbsorptionsRod];

% Convert to energy fundamentals
T_energy = EnergyToQuanta(obj.S,T_quantalAbsorptionsNormalized')';

% And normalize the energy fundamentals
T_energyNormalized = bsxfun(@rdivide,T_energy,max(T_energy, [], 2));

% Assign the fields in the receptor object
obj.T.T_quantalIsomerizations = T_quantalIsomerizations;
obj.T.T_quantalAbsorptions = T_quantalAbsorptions;
obj.T.T_quantalAbsorptionsNormalized = T_quantalAbsorptionsNormalized;
obj.T.T_energy = T_energy;
obj.T.T_energyNormalized = T_energyNormalized;
obj.labels = {'LCone' 'MCone' 'SCone' 'Mel' 'Rod'};