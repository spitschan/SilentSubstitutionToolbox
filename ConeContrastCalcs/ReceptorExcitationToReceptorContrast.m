function [contrasts, excitationDiff] = ReceptorExcitationToReceptorContrast(receptorExcitations)
% Calculates contrast on photoreceptors from receptor excitations
%
% Syntax:
%   contrasts = ReceptorExcitationToReceptorContrast(receptorExcitations)
%
% Description:
%    Takes in vectors of receptor excitations and calculates the contrast
%    on each receptor between all pairs of excitations.
%
% Inputs:
%    receptorExcitations - RxN matrix of excitations of the R receptors for
%                          each of the N vectors of excitations.
%
% Outputs:
%    contrasts           - NxNxR matrix of contrasts in percent (one NxN
%                          matrix per receptor type), where
%                          contrasts(i,j,R) = excitationDiff(i,j,R) /
%                          excitation(R,i) * 100%
%    excitationDiff      - NxNxR matrix of differences in excitations (one
%                          NxN matrix per receptor type), where
%                          excitationDiff(i,j,R) = excitation(R,j) -
%                          excitation(R,i).
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    In the case that only 2 vectors of excitations are passed (e.g., a
%    background SPD and a direction SPD), the outputs are simplified as
%    follows:
%       excitationDiff - Rx2 matrix, where the first column is the
%                        excitation(R,j) - reponse(R,i), and the second
%                        column the inverse
%       contrasts      - Rx2 matrix, where the first column is the contrast
%                        relative to the first vector of excitations, and
%                        the second column is the contrast relative to the
%                        second vector of excitations.
%
%    In the case that only 1 vector of excitations is passed,
%    excitationDiff and contrasts are returned as nWlsx1 columnvectors of
%    NaNs.
%
% See also:
%    SPDToReceptorExcitation, SPDToReceptorContrast

% History:
%    03/02/18  jv  extracted from SPDToReceptorContrast

%% Input validation
parser = inputParser;
parser.addRequired('receptorExcitations',@isnumeric);
parser.parse(receptorExcitations);

%% Calculate contrasts
if size(receptorExcitations,2) <= 1
    excitationDiff = NaN(size(receptorExcitations));
    contrasts = NaN(size(receptorExcitations));
else
    % Calculate difference in receptor excitations between all SPDs
    temp = reshape(receptorExcitations',[1,size(receptorExcitations,2),size(receptorExcitations,1)]);
    temp = repmat(temp,[size(receptorExcitations,2),1,1]);
    excitationDiff = temp - permute(temp,[2 1 3]);
    
    % Squeeze, if only 2 SPDs were passed
    if size(receptorExcitations,2) == 2
        excitationDiff = squeeze([excitationDiff(1,2,:) excitationDiff(2,1,:)])';
        % denominator for contrast is the excitations matrix when N = 2
        temp = receptorExcitations;
    else
        % denominator for contrast is a permutation of the temp matrix
        temp = permute(temp,[2 1 3]);
    end
    
    % Calculate contrasts
    contrasts = excitationDiff ./ temp * 100;
end


end