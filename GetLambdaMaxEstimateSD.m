function lambdaMaxSDnm = GetLambdaMaxEstimateSD(photopigmentClass, source)
% lambdaMaxSDnm = GetLambdaMaxEstimateSD(coneClass, source)
%
% Returns lambda-max variability as one standard deviation. The 'source'
% argument species whether only the data from W&M (1988) listed under "This
% study" will be returned, or an average across the cited studies.
%
% By default, for melanopsin, the 'L' variability is used, since it is the
% highest and puts an upper bound on the photopigment lambda-max
% variability.
%
% Taken from Webster & MacLeod (1988), Factors underlying individual
% differences in the color matches of normal observers, Table 4.
%
% Input:
%   photopigmentClass(str)  - Photopigment class of interest
%   source (str)            - Source of lambda-max SD estimates.
%
% Output:
%   lambdaMaxSDnm           - Standard deviation of lambda-max.
%
% See also:
%   GetChronicalAgeSDFromLensSD
%
% References:
%   Webster MA & MacLeod DI (1988) Factors underlying individual
%   differences in the color matches of normal observers. J Opt Soc Am A
%   5(10):1722-1735. ('M&W')
%
% 7/15/14       ms          Wrote it.
% 11/21/14      ms          Cleaned up and commented.

if strcmp(photopigmentClass, 'melanopsin')
   source = 'M&W';
   photopigmentClass = 'LCone';
end

switch source
    case 'M&W'
        switch photopigmentClass
            case {'LCone' 'LConeHemo'}
                lambdaMaxSDnm = 1.5;
            case {'MCone' 'MConeHemo'}
                lambdaMaxSDnm = 0.9;
            case {'SCone' 'SConeHemo'}
                lambdaMaxSDnm = 0.8;
            otherwise
                lambdaMaxSDnm = 1.5;
        end
    case 'Average'
        switch photopigmentClass
            case {'LCone' 'LConeHemo'}
                lambdaMaxSDnm = mean([1.5 1.4 1.3]);
            case {'MCone' 'MConeHemo'}
                lambdaMaxSDnm = mean([0.9 2.0]);
            case {'SCone' 'SConeHemo'}
                lambdaMaxSDnm = mean([0.8]);
            otherwise
                lambdaMaxSDnm = 1.5;
        end
    otherwise
        fprintf('Cone type unsupported, returning the maximum (LCone, 1.5 nm).');
        lambdaMaxSDnm = 1.5;
end