function [fractionBleachedFromIsom, fractionBleachedFromIsomHemo] = GetConeFractionBleachedFromSpectrum(S, spd, fieldSizeDegrees, observerAgeInYears, pupilDiameterMm, ...
    desiredPhotopicLuminanceCdM2, verbose)
% [fractionBleachedFromIsom, fractionBleachedFromIsomHemo] = GetConeFractionBleachedFromSpectrum(S, spd, fieldSizeDegrees, observerAgeInYears, pupilDiameterMm, ...
%   [desiredPhotopicLuminanceCdM2], [verbose])
%
% Returns the fraction bleached for LMS cones and LMS penumbral cones for a
% given spectrum, age, field size and pupil diameter.
%
% Input:
%   S                   - 'S' representations of wls
%   spd                 - Spectral power distribution in Watts per m2 per
%                         steradian. Needs to match what is specified in S
%   fieldSizeDegrees    - Size of the visual field
%   observerAgeInYears  - Observer age
%   pupilDiameterMm     - Assumed pupil diameter
%   desiredPhotopicLuminanceCdM2 - Adjust background luminance by scaling
%                         Optional, and ignored if not passed.
%   verbose             - Prints out the fraction bleached.  Default false.
%
% Output:
%   fractionBleachedFromIsom - Fraction bleached for open-field cones
%   fractionBleachedFromIsomHemo - Fraction bleached for penumbral cones
%
% The purpose of desired PhotopicLuminancCdM2 is to handle small deviations
% between a calibration and the current luminance, without redoing thw
% whole calibration.  The bleaching fractions are not highly sensitive to
% small changes in luminance, so this is just fine.
%
% Note:
%   This routine just does the computations for standard versions of hte
%   cones.  It doesn't take a shifted lambdaMax into account. We could pass
%   one more argument and do that, but at present we don't.
%
% See also:
%   ComputePhotopigmentBleaching
%
% References:
%   Rushton WA & Henry GH (1968) Bleaching and regeneration of cone
%   pigments in man. Vision Res 8(6):617-631.
%
%   Kaiser PK & Boynton RM (1996) Human Color Vision (Optical Society of
%   America, Washington, D.C.) 2nd Ed.
%
% 8/1/14    ms      Wrote it, based on code by DHB
% 11/21/14  ms      Cleaned up and commented.
% 12/12/14  dhb     Make verbose arg (added yesterday by ms) optional.
%           dhb     desiredPhotopicLuminanceCdM2 was ignored completely.
%                   Now used, but ignored if not passed.

%% Set default args
if (nargin < 7 | isempty(verbose))
    verbose = false;
end

%% Do a few conversions
backgroundSpd = spd;
radianceWattsPerM2Sr = backgroundSpd;
radianceWattsPerM2Sr(radianceWattsPerM2Sr < 0) = 0;
radianceWattsPerCm2Sr = (10.^-4)*radianceWattsPerM2Sr;
radianceQuantaPerCm2SrSec = EnergyToQuanta(S,radianceWattsPerCm2Sr);

%% Get the fraction bleached for each cone type. See
% OLGetBGConeIsomerizations for reference.

%% Load CIE functions.
load T_xyz1931
T_xyz = SplineCmf(S_xyz1931,683*T_xyz1931,S);
photopicLuminanceCdM2 = T_xyz(2,:)*radianceWattsPerM2Sr;
chromaticityXY = T_xyz(1:2,:)*radianceWattsPerM2Sr/sum(T_xyz*radianceWattsPerM2Sr);

%% Optional adjust of background luminance by scaling.
% Handles small shifts from original calibration, just by scaling.  This is close enough for purposes
% of computing fraction of pigment bleached.
if (nargin < 6 | isempty(desiredPhotopicLuminanceCdM2))
    desiredPhotopicLuminanceCdM2 = photopicLuminanceCdM2;
end
scaleFactor = desiredPhotopicLuminanceCdM2/photopicLuminanceCdM2;
radianceWattsPerM2Sr = scaleFactor*radianceWattsPerM2Sr;
radianceWattsPerCm2Sr = scaleFactor*radianceWattsPerCm2Sr;
radianceQuantaPerCm2SrSec = scaleFactor*radianceQuantaPerCm2SrSec;
photopicLuminanceCdM2 = scaleFactor*photopicLuminanceCdM2;

%% Get cone spectral sensitivities to use to compute isomerization rates
[T_cones, T_quantalIsom]  = GetHumanPhotoreceptorSS(S, {'LCone' 'MCone' 'SCone'}, fieldSizeDegrees, observerAgeInYears, pupilDiameterMm, [], [], [], []);
[T_conesHemo, T_quantalIsomHemo]  = GetHumanPhotoreceptorSS(S, {'LConeHemo' 'MConeHemo' 'SConeHemo'}, fieldSizeDegrees, observerAgeInYears, pupilDiameterMm, [], [], [], []);

%% Compute irradiance, trolands, etc.
pupilAreaMm2 = pi*((pupilDiameterMm/2)^2);
eyeLengthMm = 17;
degPerMm = RetinalMMToDegrees(1,eyeLengthMm);
irradianceWattsPerUm2 = RadianceToRetIrradiance(radianceWattsPerM2Sr,S,pupilAreaMm2,eyeLengthMm);
irradianceScotTrolands = RetIrradianceToTrolands(irradianceWattsPerUm2, S, 'Scotopic', [], num2str(eyeLengthMm));
irradiancePhotTrolands = RetIrradianceToTrolands(irradianceWattsPerUm2, S, 'Photopic', [], num2str(eyeLengthMm));
irradianceQuantaPerUm2Sec = EnergyToQuanta(S,irradianceWattsPerUm2);
irradianceWattsPerCm2 = (10.^8)*irradianceWattsPerUm2;
irradianceQuantaPerCm2Sec = (10.^8)*irradianceQuantaPerUm2Sec;
irradianceQuantaPerDeg2Sec = (degPerMm^2)*(10.^-2)*irradianceQuantaPerCm2Sec;

%% This is just to get cone inner segment diameter
photoreceptors = DefaultPhotoreceptors('CIE10Deg');
photoreceptors = FillInPhotoreceptors(photoreceptors);

%% Get isomerizations
theLMSIsomerizations = PhotonAbsorptionRate(irradianceQuantaPerUm2Sec,S, ...
    T_quantalIsom,S,photoreceptors.ISdiameter.value);
theLMSIsomerizationsHemo = PhotonAbsorptionRate(irradianceQuantaPerUm2Sec,S, ...
    T_quantalIsomHemo,S,photoreceptors.ISdiameter.value);

%% Get fraction bleached
fractionBleachedFromTrolands = ComputePhotopigmentBleaching(irradiancePhotTrolands,'cones','trolands','Boynton');
fractionBleachedFromIsom = zeros(3,1);
fractionBleachedFromIsomHemo = zeros(3,1);
for i = 1:3
    fractionBleachedFromIsom(i) = ComputePhotopigmentBleaching(theLMSIsomerizations(i),'cones','isomerizations','Boynton');
    fractionBleachedFromIsomHemo(i) = ComputePhotopigmentBleaching(theLMSIsomerizationsHemo(i),'cones','isomerizations','Boynton');
end

%% Optional printout
if verbose
    fprintf('    * Stimulus luminance used %0.1f candelas/m2\n',photopicLuminanceCdM2);
    fprintf('    * Fraction bleached computed from trolands (applies to L and M cones): %0.2f\n',fractionBleachedFromTrolands);
    fprintf('    * Fraction bleached from isomerization rates: L, %0.2f; M, %0.2f; S, %0.2f\n', ...
        fractionBleachedFromIsom(1),fractionBleachedFromIsom(2),fractionBleachedFromIsom(3));
    fprintf('    * Fraction bleached from isomerization rates: LHemo, %0.2f; MHemo, %0.2f; SHemo, %0.2f\n', ...
        fractionBleachedFromIsomHemo(1),fractionBleachedFromIsomHemo(2),fractionBleachedFromIsomHemo(3));
end