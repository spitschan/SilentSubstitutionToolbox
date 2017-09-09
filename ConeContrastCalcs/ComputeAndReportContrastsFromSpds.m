function [receptorContrasts, postreceptoralContrasts, postreceptoralStrings] = ComputeAndReportContrastsFromSpds(prefixString,receptorStrings,T_receptors,backgroundSpd,modulationSpd)
% ComputeAndReportContrastsFromSpds
%
% Usage:
%     contrasts = ComputeAndReportContrastsFromSpds(prefixString,photoreceptorClasses,T_receptors,backgroundSpd,modulationSpd)
%
% Description:
%    Report out contrasts function. Assumes that the modulationSpd is the spd
%    measured around the background.
%
%    The first three photoreceptorClasses must be the LMS cones, otherwise
%    the postreceptoralContrasts won't make sense.  If you just want
%    receptor contrasts, pass in 'doPostreceptoral' as false. The postreceptoralContrasts
%    are computed by ComputePostreceptoralContrastsFromLMSContrasts.
%
% Input:
%    prefixString - A string that gets printed in the MATLAB window before contrast values are reported.
%    photoreceptorClasses - Names of the photoreceptor classes in a cell array.
%    T_receptors - Array with the spectral sensitivities, of size nReceptors by nWls,
%                  where nReceptors is the number of receptors and N the number of 
%                  wavelength samples.
%    backgroundSpd - Background spd in column vector of length nWls.
%    modulationSpd - Modulation spd in column vector of length nWls.
%
% Output:
%    contrasts - Contrast values seen by the receptors in T_receptors
%    postreceptoralContrasts - Contrast values seen by the standard postreceptoral combinations.
%    postreceptoralStrings - Labels for the standard postreceptoral combinations.
%
% Optional key/value pairs:
%    'verbose' (logical)            - Whether to print out stuff or not.  Default true.
%    'doPostreceptoral' (logical)   - Whether to compute and report post-receptoral contrasts.  Default true.
%
% See also: ComputePostreceptoralContrastsFromLMSContrasts, ComputeAndReportContrastsFromOLPrimaries.

% 7/21/17  dhb   Put in comment placeholders and did my best.

%% Parse input
p = inputParser;
p.addParameter('verbose',true,@ischar);
p.addParameter('doPostreceptoral',true,@ischar);
p.parse;

%% Check
if p.Results.doPostreceptoral
    % Throw an error if the dimensions are inconsistent
    if (size(T_receptors, 1) < 3)
        error('Cannot compute postreceptoral contrasts with fewer than 3 receptor types');
    end
end

%% Print out some leading information
if (p.Results.verbose),  fprintf('\n<strong>%s</strong>\n', prefixString); end

%% Calculate the contrasts
backgroundReceptors = T_receptors*backgroundSpd;
differenceReceptors = T_receptors*(modulationSpd-backgroundSpd);
receptorContrasts = differenceReceptors ./ backgroundReceptors;
if (p.Results.verbose)
    for j = 1:size(T_receptors,1)
        fprintf('  * <strong>%s</strong>: contrast = %0.1f%%\n',receptorStrings{j},100*receptorContrasts(j));
    end
end

% Postreceptoral contrasts.  Assuming that first three receptor contrasts are L, M and S.
if p.Results.doPostreceptoral
    [postreceptoralContrasts, postreceptoralStrings] = ComputePostreceptoralContrastsFromLMSContrasts(receptorContrasts(1:3));
    
    % Print out postreceptoral contrasts
    if (p.Results.verbose)
        NCombinations = size(postreceptoralContrasts, 1);
        fprintf('\n');
        for ii = 1:NCombinations
            fprintf('* <strong>%s</strong>: contrast = %0.1f%%\n',postreceptoralStrings{ii},100*postreceptoralContrasts(ii));
        end
    end
end