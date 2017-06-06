function [photopicLuminanceCdM2, irradianceScotTrolands] = GetLuminanceAndTrolandsFromSpd(S, radianceWattsPerM2Sr, pupilDiameterMm, print);
% [photopicLuminanceCdM2, irradianceScotTrolands] = GetLuminanceAndTrolandsFromSpd(S, radianceWattsPerM2Sr, pupilDiameterMm, print);
%
% This function prints out light levels.
%
% 7/10/16       ms      Wrote it.

% Calculate the photopic luminance
load T_xyz1931
T_xyz = SplineCmf(S_xyz1931,683*T_xyz1931,S);
photopicLuminanceCdM2 = T_xyz(2,:)*radianceWattsPerM2Sr;
chromaticityXY = T_xyz(1:2,:)*radianceWattsPerM2Sr/sum(T_xyz*radianceWattsPerM2Sr);

% Calculate scotopic and photopic trolands
pupilAreaMm2 = pi*((pupilDiameterMm/2)^2);
eyeLengthMm = 17;
degPerMm = RetinalMMToDegrees(1,eyeLengthMm);
irradianceWattsPerUm2 = RadianceToRetIrradiance(radianceWattsPerM2Sr,S,pupilAreaMm2,eyeLengthMm);
irradiancePhotTrolands = RetIrradianceToTrolands(irradianceWattsPerUm2, S, 'Photopic', [], num2str(eyeLengthMm));
irradianceScotTrolands = RetIrradianceToTrolands(irradianceWattsPerUm2, S, 'Scotopic', [], num2str(eyeLengthMm));

if print
    fprintf('  * Luminance <strong>%0.6f</strong> cd/m2\n',photopicLuminanceCdM2);
    fprintf('  *   x=<strong>%0.2f</strong>, y=<strong>%0.2f</strong>\n',chromaticityXY(1), chromaticityXY(2));
    fprintf('  * Retinal irradiance [%g mm pupil] <strong>%0.1f</strong> ph td (<strong>%0.1f</strong> log ph td)\n',pupilDiameterMm,irradiancePhotTrolands,log10(irradiancePhotTrolands));
    fprintf('  * Retinal irradiance [%g mm pupil] <strong>%0.1f</strong> sc td (<strong>%0.1f</strong> log sc td)\n',pupilDiameterMm,irradianceScotTrolands, log10(irradianceScotTrolands)); 
end