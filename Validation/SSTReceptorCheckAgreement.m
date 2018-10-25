function status = SSTReceptorCheckAgreement
%
% Description:
%   This program checks for agreement between the old
%   'GetHumanPhotoreceptorSS' method and the new 'SSTReceptorHuman' object by
%   serially iterating through various parameters sets (varying pupil
%   diameter, field size and observer age, which are the front-end
%   parameters).
%
% Inputs:
%    None.
%
% Outputs:
%    status      -  True if no error, false otherwise.

% History:
%   8/29/17   ms      Written.
%   9/8/17    dhb     Comment out pause.  This seems useful for debugging but once it's working
%                     we just want to run through the whole thing and see all the OKs scroll by.
%   10/25/18  dhb     Make a function, return status rather than throw an
%                     error.  Makes it callable from a higher level
%                     validation function.

% Define wavelength spacing
S = [380 2 201];
wls = SToWls(S);

% Iterate over ...
status = true;
for pupilDiameterMm = [2 4 6 8] % ... pupil diameter
    for fieldSizeDeg = [2 10 20 40 60] % ... field size
        for obsAgeInYrs = [20 30 40 50 60 70 80] % ... observer age
            % [1] Spectral sensitivities through SSTReceptorHuman
            receptorObj = SSTReceptorHuman('S', S, 'verbosity', 'low', 'fieldSizeDeg', fieldSizeDeg, ...
                'obsPupilDiameterMm', pupilDiameterMm, 'obsAgeInYrs', obsAgeInYrs, 'doPenumbralConesTrueFalse', true);
            
            % [2] Spectral sensitivities throguh GetHumanPhotoreceptorSS
            [T_energyNormalized,T_quantalIsomerizations] = GetHumanPhotoreceptorSS(S, ...
                {'LConeTabulatedAbsorbance' 'MConeTabulatedAbsorbance' 'SConeTabulatedAbsorbance' ...
                'Melanopsin' 'Rods' ...
                'LConeTabulatedAbsorbancePenumbral' 'MConeTabulatedAbsorbancePenumbral' 'SConeTabulatedAbsorbancePenumbral'}, ...
                fieldSizeDeg, obsAgeInYrs, pupilDiameterMm);
            
            % Check for agreement by calculating the sum of the absolute
            % difference between the receptor sensitivities produced by
            % SSTReceptorHuman and those porudced by
            % GetHumanPhotoreceptorSS. These are in agreement, and would throw an error otherwise.
            for ii = 1:size(T_energyNormalized, 1)
                difft0(ii) = sum(abs(receptorObj.T.T_energyNormalized(ii, :) - T_energyNormalized(ii, :)));
                difft1(ii) = sum(abs(receptorObj.T.T_quantalIsomerizations(ii, :) - T_quantalIsomerizations(ii, :)));
            end
            if any(difft0 > 1e-7) || any(difft1 > 1e-7)
                status = false;
            else
                fprintf('X Pupil diam [mm] <strong>%d</strong>, field size [deg] <strong>%d</strong>, obs. age [yrs] <strong>%d</strong>: <strong>OK</strong>\n', pupilDiameterMm, fieldSizeDeg, obsAgeInYrs);
            end
            
            % Make some plots
            subplot(1, 2, 1);
            plot(wls, receptorObj.T.T_energyNormalized, '-r', 'LineWidth', 2); hold on;
            plot(wls, T_energyNormalized, '--k', 'LineWidth', 2);
            xlabel('Wavelength [nm]'); ylabel('Sensitivity');
            
            subplot(1, 2, 2);
            plot(wls, receptorObj.T.T_energyNormalized-T_energyNormalized, '-r'); hold on;
            xlabel('Wavelength [nm]'); ylabel('\Delta Sensitivity');
            %pause;
            
            % Clear the drawing
            subplot(1, 2, 1); hold off;
            subplot(1, 2, 2); hold off;
            
            clear difft0 difft1;
        end
    end
end