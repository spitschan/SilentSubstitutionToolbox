function [contrastMap, nominalLambdaMax, ageRange, lambdaMaxShiftRange] = CalculateSplatter(S, backgroundSpd, modulationSpd, photoreceptorClasses, fieldSizeDegrees, ...
    ageRange, pupilDiameterMm, lambdaMaxShiftRange, fractionBleached)
% [contrastMap, nominalLambdaMax, ageRange, lambdaMaxShiftRange] = CalculateSplatter(S, backgroundSpd, modulationSpd, photoreceptorClasses, fieldSizeDegrees, ...
%    ageRange, pupilDiameterMm, lambdaMaxShiftRange, fractionBleached)
%
% Calculates splatter on the passed photopigments. Assumes that lambda-max
% and age will be primary variables of interest. To explore effect of
% field size, pupil size or fraction bleached you need to call CalculateSplatter multiple
% times and compare the output.
%
% Let n be the number of lambda-max values and m be the number of age
% values.  If there are k photoreceptors, then contrastMap is a cell array
% with k entries, one for each photoreceptors. Each entry is an n by m
% matrix whose elements contain the contrast corresponding to each
% lambda-max/age pair.
%
% If empty variables are passed for any of the following variables,
% defaults will be assumed.
%
% Input:
%   S (1x3)                         - Wavelength spacing.
%                                     Default: [380 2 201]
%   backgroundSpd (1xlength(S(3))   - Spectral power distribution of
%                                     background.
%   modulationSpd (1xlength(S(3))   - Spectral power distribution of
%                                     modulation.
%   photoreceptorClasses (cell)     - Cell with names of photoreceptor
%                                     classes.
%   fieldSizeDegrees (1x1)          - Field size in degrees.
%                                     Default: 10
%   ageRange (1xM)                  - Vector of observer age
%                                     Default: 20:60
%   pupilDiameterMm (1x1)           - Pupil diameter in mm.
%                                     Default: 3
%   lambdaMaxShiftRange (1xN)       - Vector of lambda-max shifts
%                                     Default: -5:0.5:5
%   fractionBleached (1xN)          - Fraction of pigment bleached.
%                                     Default: 0
%
% Output:
%   contrast (cell, {K}(MxN)        - 2D contrast map, for 1...K
%                                     photoreceptor classes
%   nominalLambdaMax (1xK)          - Contains the nominal lambda max for
%                                     the k classes
%   ageRange (1xM)                  - Age range
%   lambdaMaxShiftRange (1xN)       - Shift range for lambda-max
%
% 1/21/14   ms    Wrote it based on old code.
% 5/26/14   dhb   New calling form for GetHumanPhotoreceptorSS.
% 11/21/14  ms    Cleaned up and commented

NPhotoreceptorClasses = length(photoreceptorClasses);

%%  Check if all variables have been passed with a value, otherwise assign defaults.
if isempty(S)
    S = [380 2 201];
end

if size(backgroundSpd, 2) ~= 1
    error('Passed background spd has wrong dimension.')
end

if size(backgroundSpd, 1) ~= S(3)
    error('Passed background uses different wavelength sampling than the passed S vector.')
end

if size(modulationSpd, 2) ~= 1
    error('Passed modulation spd has wrong dimension.')
end

if size(modulationSpd, 1) ~= S(3)
    error('Passed modulation uses different wavelength sampling than the passed S vector.')
end

if isempty(photoreceptorClasses)
    error('No photoreceptor classes to calculate contrast for have been passed.');
end

if isempty(lambdaMaxShiftRange)
    lambdaMaxShiftRange = -10:0.5:10;
end

if isempty(ageRange)
    ageRange = 20:60;
end

if isempty(fieldSizeDegrees)
    fieldSizeDegrees = 10;
end

if isempty(pupilDiameterMm)
    pupilDiameterMm = 3;
end

if (isempty(fractionBleached))
    fractionBleached = zeros(NPhotoreceptorClasses,1);
end


%% Iterate over photoreceptor classes (k), then over lambda-max shift (m), then over age (n)
for k = 1:NPhotoreceptorClasses
    [X, Y] = meshgrid(ageRange, lambdaMaxShiftRange);
    X = X(:);
    Y = Y(:);
    NObservers = length(X)
    
    % Preassign memory
    T_energy = zeros(NObservers, S(3));
    for i = 1:NObservers
        [T_energy(i, :), ~, nominalLambdaMax(k)] = GetHumanPhotoreceptorSS(S, {photoreceptorClasses{k}}, fieldSizeDegrees, X(i), pupilDiameterMm, Y(i), fractionBleached(k));
    end
    tmp = (T_energy*(modulationSpd-backgroundSpd))./(T_energy*backgroundSpd);
    tmp2 = reshape(tmp, length(lambdaMaxShiftRange), length(ageRange));
    
    % Finally, assign the contrast map.
    contrastMap{k} = tmp2;
end