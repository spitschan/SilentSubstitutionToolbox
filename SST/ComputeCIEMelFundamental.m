function [T_quantalAbsorptionsNormalizedMel,T_quantalAbsorptionsMel,T_quantalIsomerizationsMel,adjIndDiffParams,params,staticParams] = ComputeCIEMelFundamental(S,fieldSizeDegrees,ageInYears,pupilDiameterMm,indDiffParams)
% [T_quantalAbsorptionsNormalizedMel,T_quantalAbsorptionsMel,T_quantalIsomerizationsMel,params,staticParams] = ComputeCIEMelFundamental(S,fieldSizeDegrees,ageInYears,pupilDiameterMm,indDiffParams)
%
% This function computes the melanopsin fundamental. It uses the same logic as
% ComputeCIEConeFundamentals.
%
% 8/1/17    ms      Written.

%% Override default values so that FillInPhotoreceptors does
% our work for us.  The CIE standard uses field size,
% age, and pupil diameter to computer other values.
% to compute other quantities.
photoreceptors = DefaultPhotoreceptors('LivingHumanMelanopsin');
photoreceptors.nomogram.S = S;
photoreceptors.fieldSizeDegrees = fieldSizeDegrees;
photoreceptors.ageInYears = ageInYears;
photoreceptors.pupilDiameter.value = pupilDiameterMm;
photoreceptors = FillInPhotoreceptors(photoreceptors);

%% Set up for call into the low level routine that computes the CIE fundamentals.
staticParams.S = photoreceptors.nomogram.S;
staticParams.fieldSizeDegrees = photoreceptors.fieldSizeDegrees;
staticParams.ageInYears = photoreceptors.ageInYears;
staticParams.pupilDiameterMM = photoreceptors.pupilDiameter.value;
staticParams.lensTransmittance = photoreceptors.lensDensity.transmittance;
staticParams.macularTransmittance = photoreceptors.macularPigmentDensity.transmittance;
staticParams.quantalEfficiency = photoreceptors.quantalEfficiency.value;
staticParams.whichNomogram = photoreceptors.nomogram.source;
params.axialDensity = photoreceptors.axialDensity.value;
params.lambdaMax = photoreceptors.nomogram.lambdaMax;
params.indDiffParams = indDiffParams;
params.DORODS = true;

%% Drop into more general routine to compute
%
% See comment in ComputeRawConeFundamentals about the fact that
% we ought to unify this routine and what FillInPhotoreceptors does.
[T_quantalAbsorptionsNormalizedMel,T_quantalAbsorptionsMel,T_quantalIsomerizationsMel,adjIndDiffParams] = ComputeRawConeFundamentals(params,staticParams);