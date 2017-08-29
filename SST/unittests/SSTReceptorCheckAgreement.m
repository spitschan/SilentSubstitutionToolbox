%% SSTReceptorCheckAgreement
%
% This program checks for agreement between the old
% 'GetHumanPhotoreceptorSS' method and the new 'SSTReceptorHuman' object by
% serially iterating through various parameters sets (varying pupil
% diameter, field size and observer age, which are the front-end
% parameters).
%
% 8/29/17   ms      Written.

% Iterate over
for pupilDiameterMm = [2 4 6 8] % pupil diameter
    for fieldSizeDeg = [2 10 20 40 60] % field size
        for obsAgeInYrs = [20 30 40 50 60 70 80] % observer age
            % [1] Spectral sensitivities throguh SSTReceptorHuman
            receptorObj = SSTReceptorHuman('verbosity', 'low', 'fieldSizeDeg', fieldSizeDeg, ...
                'obsPupilDiameterMm', pupilDiameterMm, 'obsAgeInYrs', obsAgeInYrs, 'doPenumbralConesTrueFalse', true);
            
            % [2] Spectral sensitivities throguh GetHumanPhotoreceptorSS
            [T_energyNormalized,T_quantalIsomerizations] = GetHumanPhotoreceptorSS([380 2 201], ...
                {'LConeTabulatedAbsorbance' 'MConeTabulatedAbsorbance' 'SConeTabulatedAbsorbance' ...
                'Melanopsin' 'Rods' ...
                'LConeTabulatedAbsorbancePenumbral' 'MConeTabulatedAbsorbancePenumbral' 'SConeTabulatedAbsorbancePenumbral'}, ...
                fieldSizeDeg, obsAgeInYrs, pupilDiameterMm);
            %%
            % Check for agreement
            for ii = 1:size(T_energyNormalized, 1)
                difft0(ii) = sum(abs(receptorObj.T.T_energyNormalized(ii, :) - T_energyNormalized(ii, :)));
                difft1(ii) = sum(abs(receptorObj.T.T_quantalIsomerizations(ii, :) - T_quantalIsomerizations(ii, :)));
                if any(difft0 > 1e-7) || any(difft1 > 1e-7)
                    error('X Not consistent');
                else
                    fprintf('X Pupil diam [mm] <strong>%d</strong>, field size [deg] <strong>%d</strong>, obs. age [yrs] <strong>%d</strong>: <strong>OK</strong>\n', pupilDiameterMm, fieldSizeDeg, obsAgeInYrs);
                end
            end
            clear difft0 difft1;
        end
    end
end