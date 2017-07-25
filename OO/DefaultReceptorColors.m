function theRGB = DefaultReceptorColors(whichReceptors)
% theRGB = DefaultReceptorColors(whichReceptors)
%
% Returns RGB values in [0 1] for each of the receptors classes.
%
% 7/25/17   ms   Written.

% Assume some default receptors
if ~exist('whichReceptors', 'var') || isempty(whichReceptors)
   whichReceptors  = {'LCone' 'MCone' 'SCone' 'Melanopsin' 'Rod'}; 
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
        case 'Melanopsin'
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