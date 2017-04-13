% Photopigment density adjustment
PhotopigmentAxialDensity({'LCone' 'MCone' 'SCone'},'Human','CIE',27.5);
PhotopigmentAxialDensity({'LCone' 'MCone' 'SCone'},'Human','CIE',64);

% Macular density adjustment
a = MacularTransmittance([380 2 201],'Human','CIE',27.5);
b = MacularTransmittance([380 2 201],'Human','CIE',64);