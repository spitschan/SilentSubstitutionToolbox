function [contrasts postreceptoralContrasts postreceptoralStrings] = ComputeAndReportContrastsFromSpds(string,photoreceptorClasses,T_receptors,backgroundSpd,modulationSpd,postreceptoralCombinations,print)
% contrasts = ComputeAndReportContrastsFromSpds(string,photoreceptorClasses,T_receptors,backgroundSpd,modulationSpd,postreceptoralCombinations,print)
%
% Report out contrasts function. Assumes that the modulationSpd is the spd
% measured around the background.
%
% See also ComputeAndReportContrastsFromOLPrimaries.

% Do some check on the input
if (nargin < 7 | isempty(print))
    print = true;
end

if (nargin < 6 | isempty(postreceptoralCombinations))
    postreceptoralCombinations = [];
end

% Throw an error if the dimensions are inconsistent
if size(postreceptoralCombinations, 2) ~= size(T_receptors, 1)
    error('Postreceptoral combinations are not well specified. Check dimensions.');
end

% Calculate the contrasts
backgroundReceptors = T_receptors*backgroundSpd;
modulationReceptors = T_receptors*(modulationSpd-backgroundSpd);
contrasts = modulationReceptors ./ backgroundReceptors;
if (print)
    for j = 1:size(T_receptors,1)
        fprintf('\t%s, <strong>%s</strong>: contrast = %0.1f%%\n',string,photoreceptorClasses{j},100*contrasts(j));
    end
end

% Let this breathe
fprintf('\n');

% Print out postreceptoral contrasts
NCombinations = size(postreceptoralCombinations, 1);
postreceptoralContrasts = postreceptoralCombinations' \ contrasts;
for ii = 1:NCombinations
    % Assemble the string
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
    fprintf('\t%s, <strong>%s</strong>: contrast = %0.1f%%\n',string,theString,100*postreceptoralContrasts(ii));
    postreceptoralStrings{ii} = theString;
end