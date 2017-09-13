function [T_quantalAbsorptionsNormalizedRod,T_quantalAbsorptionsRod,T_quantalIsomerizationsRod,adjIndDiffParams,params,staticParams] = ComputeCIERodFundamental(S,fieldSizeDegrees,ageInYears,pupilDiameterMm,indDiffParams)
% [T_quantalAbsorptionsNormalizedRod,T_quantalAbsorptionsRod,T_quantalIsomerizationsRod,params,staticParams] = ComputeCIERodFundamental(S,fieldSizeDegrees,ageInYears,pupilDiameterMm,indDiffParams)
%
% Usage:
%     ComputeCIERodFundamental(S,fieldSizeDegrees,ageInYears,pupilDiameterMm,indDiffParams)
%
% Description:
%     This function computes the Rodanopsin fundamental. It uses the same
%     logic as ComputeCIEConeFundamentals.
%
% Input:
%     S - Wavelength specification in standard PTB 'S' format
%     fieldSizeDegrees - Field size in degrees
%     ageInYears - Assumed observer age in years
%     pupilDiameterMm - Assumed diameter in mm of the observer's pupil
%     indDiffParams - struct containing the default individual difference
%                     parameters in its field. These are:
%                           indDiffParams.dlens - lens density
%                           indDiffParams.dmac - macular pigment density
%                           indDiffParams.dphotopigment - 1x3 vector of
%                                                         photopigment
%                                                         densities for LMS
%                                                         cones
%                           indDiffParams.lambdaMaxShift - 1x3 vector of
%                                                          lambda max
%                                                          shifts for LMS
%                                                          cones
%                           indDiffParams.shiftType - shift type for
%                                                     lambda-max adjustment
%     
%
% Output:
%     T_quantalAbsorptionsNormalizedRod - Rodanopsin normalized quantal
%                                         absorptions in the same wavelength
%                                         spacing as in S
%     T_quantalAbsorptionsRod - Rodanopsin quantal absorptions in the same
%                               wavelength spacing as in S
%     T_quantalIsomerizationsRod - Rodanopsin quantal isomerizations in the
%                                  same wavelength spacing as in S
%     adjIndDiffParams - Individual difference parameters
%     params, staticParams - Params structs returned from the DefaultPhotoreceptors and
%                            FillInPhotoreceptors
%
% See also:
%     ComputeCIEMelFundamental, ComputeCIEConeFundamentals, ComputeRawConeFundamentals

%
% Optional key/value pairs:
%     None.

% 8/1/17    ms  Written.
% 9/8/17    ms  Added header comments.

%% Override default values so that FillInPhotoreceptors does
% our work for us.  The CIE standard uses field size,
% age, and pupil diameter to computer other values.
% to compute other quantities.
photoreceptors = DefaultPhotoreceptors('LivingHumanRod');
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

%% Call ComputeRawConeFundamentals
[T_quantalAbsorptionsNormalizedRod,T_quantalAbsorptionsRod,T_quantalIsomerizationsRod,adjIndDiffParams] = ComputeRawConeFundamentals(params,staticParams);