%% Nominal ones
pupilDiameterMm = 6;
fieldSizeDeg = 10;
obsAgeInYrs = 50;
% [1] Spectral sensitivities throguh SSTReceptorHuman
receptorObj = SSTReceptorHuman('verbosity', 'high', 'fieldSizeDeg', fieldSizeDeg, ...
    'obsPupilDiameterMm', pupilDiameterMm, 'obsAgeInYrs', obsAgeInYrs);

% [2] Spectral sensitivities throguh GetHumanPhotoreceptorSS
[T_energyNormalized,T_quantalIsomerizations] = GetHumanPhotoreceptorSS([380 2 201], ...
    {'LConeTabulatedAbsorbance' 'MConeTabulatedAbsorbance' 'SConeTabulatedAbsorbance' 'Melanopsin' 'Rods'}, ...
    fieldSizeDeg, obsAgeInYrs, pupilDiameterMm);
%%
% Check for agreement
fprintf('* Checking for agreement:\n');
for ii = 1:size(T_energyNormalized, 1)
    difft(ii) = sum(abs(receptorObj.T.T_energyNormalized(ii, :) - T_energyNormalized(ii, :)));
    fprintf('    - <strong>%s</strong>:\t%d\n', receptorObj.labels{ii}, difft(ii));
end