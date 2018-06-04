function makeSpectralSensitivities(obj)
% makeSpectralSensitivities
%
% Usage:
%     receptorObj.makeSpectralSensitivities;
%     makeSpectralSensitivities(receptorObj);
%
% Description:
%     This method of @SSTReceptorHuman creates the point estimate of the cone
%     fundamentals using machinery from Psychtoolbox-3. By default, the
%     LMS cone fundamentals, and the melanopsin and rod spectral sensitivities
%     are returned in the rows of the created T matrix, in that order.
%
%     The field obj.doPenumbralConesTrueFalse determines whether penumbral L*M*S*
%     fundamentals are computed and added as the last rows of the T matrix.  This
%     can be controlled by a key/value pair on the creation of the object.
%    
%     The outputs are returned to the field "T" of the receptor object, in
%     the following formats/units:
%       T.T_quantalIsomerizations - Quantal isomerizations
%       T.T_quantalAbsorptions - Quantal absorptions
%       T.T_quantalAbsorptionsNormalized - Normalized quantal absoprtions
%       T.T_energy - Energy fundamentals
%       T.T_energyNormalized - Normalized energy fundamentals
%       T.T_absorbance - Pigment absorbances
%       T.trans_lens - Lens transmittance
%       T.trans_mac - Macular transmittance
%
%     This routine gets called when the receptor object gets created by
%     @SSTReceptorHuman.
%
% Input:
%     obj - The receptorObj (e.g. from @SSTReceptor or @SSTReceptorHuman)
%
% Output:
%     None.
%
% Optional key/value pairs:
%     None.
%
% See also:
%     @SSTReceptorHuman, makeSpectralSensitivities,
%     makeSpectralSensitivitiesParametricVariation

% 7/25/17   ms  Commented.
% 9/7/17    ms  Updated header comments.

%% Generate the spectral sensitivities
%
% LMS cones
[T_quantalAbsorptionsNormalizedLMS,T_quantalAbsorptionsLMS,T_quantalIsomerizationsLMS,adjIndDiffParamsLMS] = ComputeCIEConeFundamentals(obj.S,...
    obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm);

% Melanopsin
[T_quantalAbsorptionsNormalizedMel,T_quantalAbsorptionsMel,T_quantalIsomerizationsMel,adjIndDiffParamsMel] = ComputeCIEMelFundamental(obj.S,...
    obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,[]);

% Rod
[T_quantalAbsorptionsNormalizedRod,T_quantalAbsorptionsRod,T_quantalIsomerizationsRod,adjIndDiffParamsRod] = ComputeCIERodFundamental(obj.S,...
    obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,[]);

% L*M*S* (penumbral) cones, if required
if obj.doPenumbralConesTrueFalse
    % We assume standard parameters here.
    source = 'Prahl';
    vesselOxyFraction = 0.85;
    vesselOverallThicknessUm = 5;
    trans_Hemoglobin = GetHemoglobinTransmittance(obj.S, vesselOxyFraction, vesselOverallThicknessUm, source);
    
    % Expand for the three cones
    trans_Hemoglobin = repmat(trans_Hemoglobin, 1, size(T_quantalAbsorptionsNormalizedLMS, 1));
    
    T_quantalAbsorptionsNormalizedLMSPenumbral = T_quantalAbsorptionsNormalizedLMS .* trans_Hemoglobin';
    T_quantalAbsorptionsNormalizedLMSPenumbral = bsxfun(@rdivide,T_quantalAbsorptionsNormalizedLMSPenumbral,max(T_quantalAbsorptionsNormalizedLMSPenumbral, [], 2));
    T_quantalAbsorptionsLMSPenumbral = T_quantalAbsorptionsLMS .* trans_Hemoglobin';
    T_quantalIsomerizationsLMSPenumbral = T_quantalIsomerizationsLMS .* trans_Hemoglobin';
end

%% Assemble the sensitivities
%
% Normalized quantal sensitivities
T_quantalAbsorptionsNormalized = [T_quantalAbsorptionsNormalizedLMS ; T_quantalAbsorptionsNormalizedMel ; T_quantalAbsorptionsNormalizedRod];

% Quantal isomerizations
T_quantalIsomerizations = [T_quantalIsomerizationsLMS ; T_quantalIsomerizationsMel ; T_quantalIsomerizationsRod];

% Quantal absorption
T_quantalAbsorptions = [T_quantalAbsorptionsLMS ; T_quantalAbsorptionsMel ; T_quantalAbsorptionsRod];

%% Add the penumbral cones if required
if obj.doPenumbralConesTrueFalse
    T_quantalAbsorptionsNormalized = [T_quantalAbsorptionsNormalized ; T_quantalAbsorptionsNormalizedLMSPenumbral];
    T_quantalAbsorptions = [T_quantalAbsorptions ; T_quantalAbsorptionsLMSPenumbral];
    T_quantalIsomerizations= [T_quantalIsomerizations ; T_quantalIsomerizationsLMSPenumbral];
end

%% Convert to energy fundamentals
T_energy = EnergyToQuanta(obj.S,T_quantalIsomerizations')';

% And normalize the energy fundamentals
T_energyNormalized = bsxfun(@rdivide,T_energy,max(T_energy, [], 2));

%% Assign the fields in the receptor object
obj.T.T_quantalIsomerizations = T_quantalIsomerizations;
obj.T.T_quantalAbsorptions = T_quantalAbsorptions;
obj.T.T_quantalAbsorptionsNormalized = T_quantalAbsorptionsNormalized;
obj.T.T_energy = T_energy;
obj.T.T_energyNormalized = T_energyNormalized;
obj.T.T_absorbance = [adjIndDiffParamsLMS.absorbance ; adjIndDiffParamsMel ; adjIndDiffParamsRod];
obj.T.trans_lens = adjIndDiffParamsLMS.lens;
obj.T.trans_mac = adjIndDiffParamsLMS.mac;

if obj.doPenumbralConesTrueFalse == 0
    obj.labels = {'LCone' 'MCone' 'SCone' 'Mel' 'Rod'};
elseif obj.doPenumbralConesTrueFalse == 1
    obj.labels = {'LCone' 'MCone' 'SCone' 'Mel' 'Rod' 'LConePenumbral' 'MConePenumbral' 'SConePenumbral'};
end