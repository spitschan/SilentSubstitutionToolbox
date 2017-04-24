% Photopigment density adjustment
PhotopigmentAxialDensity({'LCone' 'MCone' 'SCone'},'Human','CIE',27.5);
PhotopigmentAxialDensity({'LCone' 'MCone' 'SCone'},'Human','CIE',64);

% Macular density adjustment
a = MacularTransmittance([380 2 201],'Human','CIE',27.5);
b = MacularTransmittance([380 2 201],'Human','CIE',64);

S = [380 2 201]; wls = SToWls(S);
species = 'Human';
source = 'CIE';
pupilDiameterMM = 3;

theAges = [20:10:60];
for ageInYears = theAges
    [lensTransmit, lensDensity] = LensTransmittance(S,species,source,ageInYears,pupilDiameterMM);
    subplot(1, 2, 1); plot(wls, lensTransmit); hold on; xlabel('Wavelength [nm]'); ylabel('Lens tranmittance');
    subplot(1, 2, 2); plot(wls, lensDensity); hold on; xlabel('Wavelength [nm]'); ylabel('Lens density');
end