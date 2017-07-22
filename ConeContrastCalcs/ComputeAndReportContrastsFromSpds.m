function [contrasts, postreceptoralContrasts, postreceptoralStrings] = ComputeAndReportContrastsFromSpds(prefixString,photoreceptorClasses,T_receptors,backgroundSpd,modulationSpd,postreceptoralCombinations,print)
% ComputeAndReportContrastsFromSpds
%
% Usage:
%     contrasts = ComputeAndReportContrastsFromSpds(prefixString,photoreceptorClasses,T_receptors,backgroundSpd,modulationSpd,postreceptoralCombinations,print)
%
% Description:
%    Report out contrasts function. Assumes that the modulationSpd is the spd
%    measured around the background.
%
% Input:
%    prefixString - A string that gets printed in the MATLAB window before
%                   contrast values are reported.
%    photoreceptorClasses - Names of the photoreceptor classes
%    T_receptors - Array with the spectral sensitivities, of size KxN,
%                  where K is the number of receptors and N the number of 
%                  wavelength samples
%    backgroundSpd - Background spd of size N
%    modulationSpd - Modulation spd of size N
%    postreceptoralCombinations - Matrix of how the receptors should be
%                                 combined into postreceptoral combinations
%    print - Logical flag that determines whether the contrast should be
%            printed out
%
% Output:
%    contrasts - Contrast values seen by the receptors in T_receptors
%    postreceptoralContrasts - Contrast values seen by the postreceptoral
%                              combinations of the receptors in T_receptors
%    postreceptoralStrings - Labels for the postreceptoral combinations
%
% Optional key/value pairs:
%   None.
%
% See also: ComputeAndReportContrastsFromOLPrimaries.

% 7/21/17  dhb   Put in comment placeholders and did my best.

% Do some check on the input
if (nargin < 7 | isempty(print))
    print = true;
end

if (nargin < 6 | isempty(postreceptoralCombinations))
    postreceptoralContrasts = [];
    postreceptoralStrings = {''};
    DO_POSTRECEPTORAL = false;
else
    DO_POSTRECEPTORAL = true;
end

if DO_POSTRECEPTORAL
    % Throw an error if the dimensions are inconsistent
    if size(postreceptoralCombinations, 2) ~= size(T_receptors, 1)
        error('Postreceptoral combinations are not well specified. Check dimensions.');
    end
end

% Print out some information
fprintf('\n<strong>%s</strong>\n', prefixString);

% Calculate the contrasts
backgroundReceptors = T_receptors*backgroundSpd;
modulationReceptors = T_receptors*(modulationSpd-backgroundSpd);
contrasts = modulationReceptors ./ backgroundReceptors;
if (print)
    for j = 1:size(T_receptors,1)
        fprintf('  * <strong>%s</strong>: contrast = %0.1f%%\n',photoreceptorClasses{j},100*contrasts(j));
    end
end

% Let this breathe
fprintf('\n');

% Postreceptoral contrasts.
if DO_POSTRECEPTORAL
    % Print out postreceptoral contrasts
    NCombinations = size(postreceptoralCombinations, 1);
    postreceptoralContrasts = postreceptoralCombinations' \ contrasts;
    for ii = 1:NCombinations
        % Assemble the string that describes the psotreceptoral combination
        theReceptors = find(postreceptoralCombinations(ii, :));
        theString = [];
        for jj = 1:length(theReceptors)
            if jj > 1
                if sign(postreceptoralCombinations(ii, theReceptors(jj))) == 1
                    theString = [theString ' + ' photoreceptorClasses{theReceptors(jj)}];
                elseif sign(postreceptoralCombinations(ii, theReceptors(jj))) == -1
                    theString = [theString ' - ' photoreceptorClasses{theReceptors(jj)}];
                end
            else
                theString = [theString photoreceptorClasses{theReceptors(jj)}];
            end
        end
        
        % Print out the contrasts
        fprintf('* <strong>%s</strong>: contrast = %0.1f%%\n',theString,100*postreceptoralContrasts(ii));
        postreceptoralStrings{ii} = theString;
    end
end