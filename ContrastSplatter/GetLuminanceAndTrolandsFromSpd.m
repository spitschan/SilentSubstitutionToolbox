function [photopicLuminanceCdM2, irradianceScotTrolands] = GetLuminanceAndTrolandsFromSpd(S, radianceWattsPerM2Sr, pupilDiameterMm);

load T_xyz1931
T_xyz = SplineCmf(S_xyz1931,683*T_xyz1931,S);
photopicLuminanceCdM2 = T_xyz(2,:)*radianceWattsPerM2Sr;

pupilAreaMm2 = pi*((pupilDiameterMm/2)^2);
eyeLengthMm = 17;
degPerMm = RetinalMMToDegrees(1,eyeLengthMm);
irradianceWattsPerUm2 = RadianceToRetIrradiance(radianceWattsPerM2Sr,S,pupilAreaMm2,eyeLengthMm);
irradianceScotTrolands = RetIrradianceToTrolands(irradianceWattsPerUm2, S, 'Scotopic', [], num2str(eyeLengthMm));