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
%    contrasts           - NxNxR matrix of contrasts (one NxN
%                          matrix per receptor type), where
%                          contrasts(i,j,R) = excitationDiff(i,j,R) /
%                          excitation(R,i)
%    excitationDiff      - NxNxR matrix of differences in excitations (one
%                          NxN matrix per receptor type), where
%                          excitationDiff(i,j,R) = excitation(R,j) -
%                          excitation(R,i).
%
% Optional key/value pairs:
%    None.
%
% Examples are provided in the source code.
%
% Notes:
%    In the case that only 2 vectors of excitations are passed (e.g., a
%    background SPD and a direction SPD), the outputs are simplified as
%    follows:
%       excitationDiff - Rx2 matrix, where the first column is the
%                        excitation(R,j) - excitation(R,i), and the second
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

% Examples:
%{
    %% Two vectors of excitations on 3 receptors:
    % With just two vectors of excitations, output is simplified
    excitations = [.5, 3;...
                   .5, 3;...
                   .5, 3;];
    contrasts = ReceptorExcitationToReceptorContrast(excitations);
    
    % first column is contrast expressed relative to excitations(:,1):
    % this should be (3-.5)/.5 = 4.5/.5 = 5:
    contrasts(:,1);
    
    % second column is contrast expressed relative to excitations(:,2):
    % this should be (.5-3)/3 = -2.5/3 = -.8333:
    contrasts(:,2);
%}
%{
    %% Two vectors of excitations on 8 receptors:
    excitations = [2.6495    2.6495;...
                   2.1163    2.1163;...
                   0.2272    0.2272;...
                   0.6161    1.5224;...
                   1.0450    1.8510;...
                   2.6488    2.7850;...
                   2.1552    2.2503;...
                   0.2385    0.2530;];
    contrasts = ReceptorExcitationToReceptorContrast(excitations);
%}
%{
    %% Three columnvectors of excitations on 3 receptors:
    excitations = [1 1 1;...
                   1 2 4;...
                   1 1 3;];
    contrasts = ReceptorExcitationToReceptorContrast(excitations);

    % Each row (i,:,r) is the contrast on receptor r relative to the vector
    % of excitations i.
    % For the first receptor, this should be [0 0 0] (since the excitation
    % stays constant on this receptor):
    contrasts(1,:,1);
    % For the second receptor and relative to i = 1, this should be [0 1 3]
    contrasts(1,:,2);
    % To get all contrasts relative to excitations(:,1):
    reshape(contrasts(1,:,:),[3,3,1])';
    % where the first column is excitations(:,1) compared to
    % excitations(:,1), the second column is excitations(:,2) compared to
    % excitations(:,1), etc.

    % Each (:,:,i) panel is the full contrast matrix on receptor i. The
    % diagonal should always be 0's, since it represents contrast between a
    % a vector of excitations and itself.
    % For the first receptor, this panel should be all zeros:
    contrasts(:,:,1);
    % For the second receptor, it should be: [  0    1   3;
    %                                         -.5    0   1;
    %                                         -.75 -.50  0;]
    contrasts(:,:,2);
%}


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
        excitationDiff = squeeze([excitationDiff(1,2,:), excitationDiff(2,1,:)])';
        if size(receptorExcitations,1) == 1
            excitationDiff = excitationDiff';
        end
        % denominator for contrast is the excitations matrix when N = 2
        temp = receptorExcitations;
    else
        % denominator for contrast is a permutation of the temp matrix
        temp = permute(temp,[2 1 3]);
    end
    
    % Calculate contrasts
    contrasts = excitationDiff ./ temp;
end


end