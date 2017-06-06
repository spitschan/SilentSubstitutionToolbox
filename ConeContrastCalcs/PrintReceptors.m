function PrintReceptors(theConeStruct)
% PrintReceptors(theConeStruct)
%
% Function to print information about the cone spectral sensitivities.
fprintf('\tField size degrees: %0.1f\n', theConeStruct.fieldSizeDegrees);
fprintf('\tObserver age years: %0.1f\n', theConeStruct.observerAgeInYears);
fprintf('\tPupil diameter mm: %0.1f\n', theConeStruct.pupilDiameterMm);
fprintf('\tL cone lambda max shift: %0.1f\n', theConeStruct.Lshift);
fprintf('\tM cone lambda max shift: %0.1f\n', theConeStruct.Mshift);
