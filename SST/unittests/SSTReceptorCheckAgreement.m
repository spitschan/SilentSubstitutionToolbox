for pupilDiameterMm = [2:1:9];
    for fieldSizeDeg = [1:1:60]
        for obsAgeInYrs = [20:80];
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
            %fprintf('* Checking for agreement:\n');
            for ii = 1:size(T_energyNormalized, 1)
                difft0(ii) = sum(abs(receptorObj.T.T_energyNormalized(ii, :) - T_energyNormalized(ii, :)));
                difft1(ii) = sum(abs(receptorObj.T.T_quantalIsomerizations(ii, :) - T_quantalIsomerizations(ii, :)));
                if any(difft0 > 1e-7) || any(difft1 > 1e-7)
                    error('Not consistent');
                else
                    %fprintf('    - <strong>%s</strong>:\t%d [energy], %d [isomer.]\n', receptorObj.labels{ii}, difft0(ii), difft1(ii));
                end
            end
            clear difft0 difft1;
        end
    end
end