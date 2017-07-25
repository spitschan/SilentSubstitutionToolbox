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

% LMS cones
[T_quantalAbsorptionsNormalized,T_quantalAbsorptions,T_quantalIsomerizations] = ComputeCIEConeFundamentals(obj.S,...
    obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,[],[],[], ...
    false,[],[],[]);

% Melanopsin
photoreceptorMel = DefaultPhotoreceptors('LivingHumanMelanopsin');
photoreceptorMel.nomogram.S = obj.S;
photoreceptorMel.fieldSizeDegrees = obj.fieldSizeDeg;
photoreceptorMel.ageInYears = obj.obsAgeInYrs;
photoreceptorMel.pupilDiameter.value = obj.obsPupilDiameterMm;
photoreceptorMel = FillInPhotoreceptors(photoreceptorMel);

% Rod
photoreceptorRod = DefaultPhotoreceptors('LivingHumanRod');
photoreceptorRod.nomogram.S = obj.S;
photoreceptorRod.fieldSizeDegrees = obj.fieldSizeDeg;
photoreceptorRod.ageInYears = obj.obsAgeInYrs;
photoreceptorRod.pupilDiameter.value = obj.obsPupilDiameterMm;
photoreceptorRod = FillInPhotoreceptors(photoreceptorRod);

% Normalized quantal sensitivities
T_quantalAbsorptionsNormalized = [T_quantalAbsorptionsNormalized ; photoreceptorMel.quantalFundamentals ; photoreceptorRod.quantalFundamentals];

% Quantal isomerizations
T_quantalIsomerizations = [T_quantalIsomerizations ; photoreceptorMel.isomerizationAbsorptance ; photoreceptorRod.isomerizationAbsorptance];

% Quantal absorption
T_quantalAbsorptions = [T_quantalAbsorptions ; photoreceptorMel.effectiveAbsorptance ; photoreceptorRod.effectiveAbsorptance];

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
obj.labels = {'LCone' 'MCone' 'SCone' 'Melanopsin' 'Rod'};