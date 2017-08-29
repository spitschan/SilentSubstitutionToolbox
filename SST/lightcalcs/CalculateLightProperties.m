function out = CalculateLightProperties(S, spd, pupilDiameterMm)
% out = CalculateLightProperties(S, spd, pupilDiameterMm)
%
% Calculates various quantities based on a radiance spectrum.
%
% 8/29/17   ms      Written.

%% Start with the spd
out.spd.radianceWattsPerM2Sr = spd;
out.spd.radianceWattsPerM2Sr(out.spd.radianceWattsPerM2Sr < 0) = 0;
out.spd.radianceWattsPerCm2Sr = (10.^-4)*out.spd.radianceWattsPerM2Sr;
out.spd.radianceQuantaPerCm2SrSec = EnergyToQuanta(S,out.spd.radianceWattsPerCm2Sr);

%% Load CIE functions.
load T_xyz1931
T_xyz = SplineCmf(S_xyz1931,683*T_xyz1931,S);
out.photopicLuminanceCdM2 = T_xyz(2,:)*out.spd.radianceWattsPerM2Sr;
out.chromaticityXY = T_xyz(1:2,:)*out.spd.radianceWattsPerM2Sr/sum(T_xyz*out.spd.radianceWattsPerM2Sr);

%% Load cone spectral sensitivities
load T_cones_ss2
T_cones = SplineCmf(S_cones_ss2,T_cones_ss2,S);

%% Load in a directly measured sunlight through window
% and off piece of white paper towel on floor for comparison.
% Surely that is safe to look at.
load spd_phillybright
spd_phillybright = SplineSpd(S_phillybright,spd_phillybright,S);
photopicLuminancePhillyBrightCdM2 = T_xyz(2,:)*spd_phillybright;
OLSLratio = out.spd.radianceWattsPerM2Sr./spd_phillybright;

%% Compute irradiance, trolands, etc.
if ~exist('pupilDiameterMm', 'var') | isempty(pupilDiameterMm)
    pupilDiameterMm = 4.7;
    pupilDiameterMm = GetWithDefault('Enter observer pupil diameter in mm',pupilDiameterMm);
end
pupilAreaMm2 = pi*((pupilDiameterMm/2)^2);
eyeLengthMm = 17;
degPerMm = RetinalMMToDegrees(1,eyeLengthMm);
out.spd.irradianceWattsPerUm2 = RadianceToRetIrradiance(out.spd.radianceWattsPerM2Sr,S,pupilAreaMm2,eyeLengthMm);
out.irradianceScotTrolands = RetIrradianceToTrolands(out.spd.irradianceWattsPerUm2, S, 'Scotopic', [], num2str(eyeLengthMm));
out.irradiancePhotTrolands = RetIrradianceToTrolands(out.spd.irradianceWattsPerUm2, S, 'Photopic', [], num2str(eyeLengthMm));
out.spd.irradianceQuantaPerUm2Sec = EnergyToQuanta(S,out.spd.irradianceWattsPerUm2);
out.spd.irradianceWattsPerCm2 = (10.^8)*out.spd.irradianceWattsPerUm2;
out.spd.irradianceQuantaPerCm2Sec = (10.^8)*out.spd.irradianceQuantaPerUm2Sec;
out.spd.irradianceQuantaPerDeg2Sec = (degPerMm^2)*(10.^-2)*out.spd.irradianceQuantaPerCm2Sec;

%% Pupil adjustment factor for Ansi MPE
mpepupilDiameterMm = 3;
%mpepupilDiameterMm  = GetWithDefault('Enter ANSI 2007 MPE caclulations assumed pupil diameter in mm',mpepupilDiameterMm );
pupilAdjustFactor = (pupilDiameterMm/mpepupilDiameterMm).^2;

%% Get trolands another way.  For scotopic trolands, this just uses scotopic vlambda (in PTB as T_rods)
% and the magic factor of 1700 scotopic lumens per Watt from Wyszecki & Stiles (2cd edition),
% p. 257.  (This is the analog of 683 photopic lumens per Watt.  Then apply formula from
% page 103 of same book.
%
% Same idea for photopic trolands, although we already have luminance in cd/m2 from above so
% we can short cut a little.
%
% The agreement is good to integer scotopic trolands and I'm will to write off the rest
% as round off error.
load T_rods
T_scotopicVlambda = SplineCmf(S_rods,T_rods,S);
out.irradianceScotTrolands_check = pupilAreaMm2*1700*(T_scotopicVlambda*out.spd.radianceWattsPerM2Sr);
out.irradiancePhotTrolands_check = pupilAreaMm2*out.photopicLuminanceCdM2;

%% Get cone coordinates from radiance, and also adjust by pupil area.
% Useful for comparing to light levels produced by monochromatic lights
% in other papers
%theLMS = T_cones*out.spd.radianceWattsPerM2Sr;
%theLMSTimesPupilArea = pupilAreaMm2*theLMS;

%% Compute irradiance arriving at cornea
%
% According to OSA Handbook of Optics, 2cd Edition, Chaper 24 (vol 2), pp. 24.13-24.15, the
% conversion is (assuming some approximations), irradiance = radiance*stimulusArea/distance^2.
% This is implemented in RadianceAndDistanceAreaToCornIrradiance
stimulusRadiusMm = 6;
stimulusDistanceMm = 25;
stimulusRadiusM = stimulusRadiusMm/1000;
stimulusAreaM2 = pi*(stimulusRadiusM^2);
stimulusDistanceM = stimulusDistanceMm/1000;
stimulusRadiusDeg = rad2deg(stimulusRadiusMm/stimulusDistanceMm);
stimulusAreaDegrees2 = pi*(stimulusRadiusDeg^2);
out.spd.cornealIrradianceWattsPerM2 = RadianceAndDistanceAreaToCornIrradiance(out.spd.radianceWattsPerM2Sr,stimulusDistanceM,stimulusAreaM2);
out.spd.cornealIrradianceWattsPerCm2 = (10.^-4)*out.spd.cornealIrradianceWattsPerM2;
out.spd.cornealIrradianceQuantaPerCm2Sec = EnergyToQuanta(S,out.spd.cornealIrradianceWattsPerCm2);

%% Calculate derivatives of the spectra
out.sumRadianceWattsPerM2Sr = sum(out.spd.radianceWattsPerM2Sr);
out.log10SumRadianceWattsPerM2Sr = log10(sum(out.spd.radianceWattsPerM2Sr));

out.sumRadianceWattsPerCm2Sr = sum(out.spd.radianceWattsPerCm2Sr);
out.log10SumRadianceWattsPerCm2Sr = log10(sum(out.spd.radianceWattsPerCm2Sr));

out.sumIrradianceWattsPerCm2 = sum(out.spd.irradianceWattsPerCm2);
out.log10SumIrradianceWattsPerCm2 = log10(sum(out.spd.irradianceWattsPerCm2));

out.sumIrradianceQuantaPerCm2Sec = sum(out.spd.irradianceQuantaPerCm2Sec);
out.log10SumIrradianceQuantaPerCm2Sec = log10(sum(out.spd.irradianceQuantaPerCm2Sec));

out.sumIrradianceQuantaPerDeg2Sec = sum(out.spd.irradianceQuantaPerDeg2Sec);
out.log10SumIrradianceQuantaPerDeg2Sec = log10(sum(out.spd.irradianceQuantaPerDeg2Sec));

out.sumCornealIrradianceWattsPerCm2 = sum(out.spd.cornealIrradianceWattsPerCm2);
out.log10SumCornealIrradianceWattsPerCm2 = log10(sum(out.spd.cornealIrradianceWattsPerCm2));

out.sumCornealIrradianceQuantaPerCm2Sec = sum(out.spd.cornealIrradianceQuantaPerCm2Sec);
out.log10SumCornealIrradianceQuantaPerCm2Sec = log10(sum(out.spd.cornealIrradianceQuantaPerCm2Sec));

%% Melanopic analysis
% We have retinal and corneal spectral irradiances above and calculate a
% melanopsin-weighted melanopic irradiance from the quantities irradianceQuantaPerCm2Sec and cornealIrradianceWattsPerCm2
melanopsinAssumedFieldSizeDeg = 10;
melanopsonAssumedAgeYears = 32;

% Get the spectral sensitivities from the 
receptorObj = SSTReceptorHuman('verbosity', 'low', 'obsAgeInYrs', melanopsonAssumedAgeYears, 'fieldSizeDeg', melanopsinAssumedFieldSizeDeg);
T_melanopsinQuantal = receptorObj.T.T_quantalAbsorptionsNormalized(4, :)

% Retinal irradiance
out.sumMelIrradianceQuantaPerCm2Sec = T_melanopsinQuantal*out.spd.irradianceQuantaPerCm2Sec;
out.log10SumMelIrradianceQuantaPerCm2Sec = log10(out.sumMelIrradianceQuantaPerCm2Sec);

% Corneal irradiance
out.sumMelCornealIrradianceQuantaPerCm2Sec = T_melanopsinQuantal*out.spd.cornealIrradianceQuantaPerCm2Sec;
out.log10SumMelCornealIrradianceQuantaPerCm2Sec = log10(out.sumMelCornealIrradianceQuantaPerCm2Sec);