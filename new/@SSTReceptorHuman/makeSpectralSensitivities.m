function obj = makeSpectralSensitivities(obj)
% First three are cones
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

% Convert to energy
T_energy = EnergyToQuanta(obj.S,T_quantalAbsorptionsNormalized')';
T_energyNormalized = bsxfun(@rdivide,T_energy,max(T_energy, [], 2));

% Assign the fields
obj.T.T_quantalAbsorptionsNormalized = T_quantalAbsorptionsNormalized;
obj.T.T_quantalAbsorptions = T_quantalAbsorptions;
obj.T.T_quantalIsomerizations = T_quantalIsomerizations;
obj.T.T_energy = T_energy;
obj.T.T_energyNormalized = T_energyNormalized;