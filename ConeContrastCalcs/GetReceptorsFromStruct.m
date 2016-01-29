function T_receptors = GetReceptorsFromStruct(S, photoreceptorClasses, theConeStruct)
% T_receptors = GetReceptorsFromStruct(S, photoreceptorClasses, theConeStruct)
%
% Convert the cone struct to a set of cone spectral sensitivities.

theConeStruct.lambdaMaxShift = [theConeStruct.Lshift theConeStruct.Mshift theConeStruct.Sshift theConeStruct.Melshift];
theConeStruct.fractionBleached = [0 0 0 0];
theConeStruct.oxygenationFraction = [];
theConeStruct.vesselThickness = [];
T_receptors = GetHumanPhotoreceptorSS(S, photoreceptorClasses, ...
    theConeStruct.fieldSizeDegrees, theConeStruct.observerAgeInYears, theConeStruct.pupilDiameterMm, theConeStruct.lambdaMaxShift, ...
    theConeStruct.fractionBleached, theConeStruct.oxygenationFraction, theConeStruct.vesselThickness);

switch (theConeStruct.type)
    case 'SmithPokorny10'
        [T_cones_sp10,S_cones_sp10] = GetSmithPokornyCones10(S,true);
        T_receptors(1:3,:) = T_cones_sp10;
    case 'StockmanSharpe10'
        load T_cones_ss10
        T_receptors(1:3,:) = SplineCmf(S_cones_ss10,T_cones_ss10,S);
end
