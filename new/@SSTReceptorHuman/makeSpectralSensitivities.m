function obj = makeSpectralSensitivities(obj)
% First three are cones
[T_quantalAbsorptionsNormalized,T_quantalAbsorptions,T_quantalIsomerizations] = ComputeCIEConeFundamentals(obj.S,...
    obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,[],[],[], ...
    false,[],[],[]);
T_energy = EnergyToQuanta(obj.S,T_quantalAbsorptionsNormalized')';
T_energyNormalized = bsxfun(@rdivide,T_energy,max(T_energy, [], 2));

% obj = makeSpectralSensitivities(obj)
obj.T.T_quantalAbsorptionsNormalized = T_quantalAbsorptionsNormalized;
obj.T.T_quantalAbsorptions = T_quantalAbsorptions;
obj.T.T_quantalIsomerizations = T_quantalIsomerizations;
obj.T.T_energy = T_energy;
obj.T.T_energyNormalized = T_energyNormalized;