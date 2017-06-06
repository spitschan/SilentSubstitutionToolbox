% SSTReceptorCheckSS10Deg
load T_cones_ss10
T_cones_ss10 = SplineCmf(S_cones_ss10, T_cones_ss10, wls);

observerAgeInYrs = 32;
fractionBleached = [0 0 0];
pupilDiameterMm = 3;
fieldSizeDegrees = 10;
receptorObj = SSTReceptorHuman('obsAgeYrs', observerAgeInYrs, 'fieldSizeDeg', fieldSizeDegrees, 'obsPupilDiameterMm', pupilDiameterMm);

plot(receptorObj.T.T_energyNormalized');
hold on;
plot(T_cones_ss10', '-r');