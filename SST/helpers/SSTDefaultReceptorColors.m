function theRGB = SSTDefaultReceptorColors(whichReceptors)
% SSTDefaultReceptorColors(whichReceptors)
%
% Usage:
%     theRGB = SSTDefaultReceptorColors(whichReceptors)
%
% Description:
%     This returns RGB triplets to plot spectral sensitivities and other
%     quantities associated with receptors. If an unknown receptor is
%     passed, [0 0 0] is returned.
%
% Input:
%     whichReceptors - String or cell containing the receptor labels
%                      Possible options are: 'LCone', 'MCone', 'SCone',
%                      'Melanopsin', 'Rod', 'LConePenumbral',
%                      'MConePenumbral' and 'SConePenumbral'.
%
% Output:
%     theRGB - RGB triplets.
%
% Optional key/value pairs:
%     None.

% 9/9/17  ms  Added header comments.

% Assume some default receptors
if ~exist('whichReceptors', 'var') || isempty(whichReceptors)
    whichReceptors  = {'LCone' 'MCone' 'SCone' 'Mel' 'Rod'};
end

% Reformat the whichReceptors variable if only has been passed
if ischar(whichReceptors)
    whichReceptors = {whichReceptors};
end

% Iterate over the receptors and assign the RGB values
for ii = 1:length(whichReceptors)
    switch whichReceptors{ii}
        case 'LCone'
            theRGB(ii, :) = [0.8941 0.1020 0.1098];
        case 'MCone'
            theRGB(ii, :) = [0.3020 0.6863 0.2902];
        case 'SCone'
            theRGB(ii, :) = [0.2157 0.4941 0.7216];
        case 'Mel'
            theRGB(ii, :) = [0.1255 0.6549 0.7686];
        case 'Rod'
            theRGB(ii, :) = [0.0157 0.4118 0.3333];
            
            % Special types
        case 'LConePenumbral'
            theRGB(ii, :) = [0.6706 0.0765 0.0823];
        case 'MConePenumbral'
            theRGB(ii, :) = [0.2265 0.5147 0.2177];
        case 'SConePenumbral'
            theRGB(ii, :) = [0.1618 0.3706 0.5412];
        otherwise
            % Just use black for unknown receptor types.
            theRGB(ii, :) = [0 0 0];
    end
end