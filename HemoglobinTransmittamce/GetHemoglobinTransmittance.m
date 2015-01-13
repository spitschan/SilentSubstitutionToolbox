function [transmittance,absorptance] = GetHemoglobinTransmittance(S,oxyFraction,overallThicknessUm,source)
% [transmittance,absorbtance] = GetHemoglobinTransmittance(S,oxyFraction,overallThicknessUm,source)
% 
% Return hemoglobin transmittance = 10.^(absorptance), given a fraction of oxy/deoxy hemoglobin
% and an an overall thickness of the blood layer.
%
% Source:
%   'Prahl' - Based on data posted by Prahl and corresponding constants. [Default]
%
% See HemoglobinTransmissivityDemo for a layout of the calculations that are encapsulated
% here as well as an example of this in use.
%
% See also HemoglobinTransmissivityDemo, GetHemoglobin.
%
% 6/4/14  dhb  Wrote it, based on what ms and I worked out in HemoglobinTransmissivityDemo.

%% Set source
if (nargin < 4 | isempty(source))
    source = 'Prahl';
end

switch (source)
    case 'Prahl'
        % Obtain the Prahl hemoglobin estimates.
        %
        % These are tabulated in the function GetHemoglobinPrahl, and are returned as
        % molar extinction coefficients (cm-1/M), whatever an 'M' is.
        [S_Prahl, wls_Prahl, oxyMolarExtinction_Prahl, deoxyMolarExtinction_Prahl] = GetHemoglobin('Prahl');
        
        % Prahl writes:
        %     To convert this data to absorption coefficient in (cm-1), multiply by the molar concentration and 2.303,
        %         ua = (2.303) e (x g/liter)/(64,500 g Hb/mole)
        %     where x is the number of grams per liter. A typical value of x for whole blood is x=150 g Hb/liter.
        absorptivityPerCm_oxy_Prahl = (2.303 .* oxyMolarExtinction_Prahl * 150)/64500;
        absorptivityPerCm_deoxy_Prahl = (2.303 .* deoxyMolarExtinction_Prahl * 150)/64500;
        
        % Calculate absorptance for oxy and deoxy.  We split the total thickness into effective
        % layers, one for oxy and one for deoxy.
        UmToCm = 1e-4;
        absorptance_oxy_Prahl = absorptivityPerCm_oxy_Prahl*overallThicknessUm*UmToCm;
        absorptance_deoxy_Prahl = absorptivityPerCm_deoxy_Prahl*overallThicknessUm*UmToCm;
        absorptance = oxyFraction*absorptance_oxy_Prahl + (1-oxyFraction)*absorptance_deoxy_Prahl;
          
        % Spline to specified S
        absorptance = SplineSrf(S_Prahl,absorptance,S);
        
        % Get transmittance.
        transmittance = 10.^(-absorptance);
    otherwise
        error('Unknown source specified');
end