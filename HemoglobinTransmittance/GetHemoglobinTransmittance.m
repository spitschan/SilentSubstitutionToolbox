function [transmittance,absorptance,S] = GetHemoglobinTransmittance(S,oxyFraction,overallThicknessUm,source)
% [transmittance,absorbtance] = GetHemoglobinTransmittance(S,oxyFraction,overallThicknessUm,source)
% 
% Return hemoglobin transmittance = 10.^(absorptance), given a fraction of oxy/deoxy hemoglobin
% and an an overall thickness of the blood layer.
%
% S may be passed as empty, in which case the S at which the underlying
% data are tabulated is used and returned.
%
% Source:
%   'Prahl' - Based on data posted by Prahl and corresponding constants. 
%             See  http://omlc.org/spectra/hemoglobin/.  [Default]
%
% See HemoglobinTransmittanceDemo for a layout of the calculations that are encapsulated
% here as well as an example of this in use.
%
% See also HemoglobinTransmittanceDemo, GetHemoglobin.
%
% 6/4/14  dhb  Wrote it, based on what ms and I worked out in HemoglobinTransmittanceDemo.
% 1/15/14 ms   Fixed error of 2.303 constant, which assumed our absorption
%              coefficient was in natural log. ln 10 = 2.303.

%% Set source
if (nargin < 4 | isempty(source))
    source = 'Prahl';
end

switch (source)
    case 'Prahl'
        % Obtain the Prahl hemoglobin estimates.  See
        % http://omlc.org/spectra/hemoglobin/.
        %
        % These are tabulated in the function GetHemoglobinPrahl, and are returned as
        % molar extinction coefficients (cm-1/M), whatever an 'M' is.
        [S_Prahl, wls_Prahl, oxyMolarExtinction_Prahl, deoxyMolarExtinction_Prahl] = GetHemoglobin('Prahl');
        
        % Prahl writes:
        %     To convert this data, denoted by e for extinction coeffecient) to absorption coefficient, denoted by ua, in (cm-1),
        %     multiply e by the molar concentration x (in grams/liter) and
        %     the pathlength l (in cm).  This whole thing then gets divided  by the constant 64,500 (g/mole) which we think
        %     describes hemoglobin.
        %     Prahl also has an extra factor of 2.303 in his example calculation.
        %         ua = (2.303) * e * x (g/liter) * l (cm) / 64,500 (g/mole)
        %     x for whole blood is x = 150 g/liter.
        %
        % We note that the 2.303 constant is in his formulat to convert the
        % absorption coefficients from decimal log units to natural log
        % units.  We want to stay in decimal units, so we omit that factor
        % in our calculations.  This factor seems to float in and out of
        % calculations we have seen in the literature because which type of
        % units (decimal or natural log) are being used seems to vary from
        % one lab to the next without it being explicitly stated.
        absorptivityPerCm_oxy_Prahl = (oxyMolarExtinction_Prahl * 150)/64500;
        absorptivityPerCm_deoxy_Prahl = (deoxyMolarExtinction_Prahl * 150)/64500;
        
        % Calculate absorptance for oxy and deoxy.  We split the total thickness into effective
        % layers, one for oxy and one for deoxy.
        UmToCm = 1e-4;
        absorptance_oxy_Prahl = absorptivityPerCm_oxy_Prahl*overallThicknessUm*UmToCm;
        absorptance_deoxy_Prahl = absorptivityPerCm_deoxy_Prahl*overallThicknessUm*UmToCm;
        absorptance = oxyFraction*absorptance_oxy_Prahl + (1-oxyFraction)*absorptance_deoxy_Prahl;
          
        % Spline to specified S, if S passed
        if (~isempty(S))
            absorptance = SplineSrf(S_Prahl,absorptance,S);
        else
            S = S_Prahl;
        end
        
        % Get transmittance.
        transmittance = 10.^(-absorptance);
    otherwise
        error('Unknown source specified');
end
