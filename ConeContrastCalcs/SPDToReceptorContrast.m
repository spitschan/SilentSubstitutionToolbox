function [contrasts, response, responseDiff] = SPDToReceptorContrast(SPDs,SSTReceptors)
% Calculates the contrast on photoreceptors between SPDs
%
% Syntax:
%   contrasts = SPDtoReceptorContrast(SPDs,SSTReceptors)
%
% Description:
%    SPDToReceptorContrast takes in spectra, as a set of columnvectors of
%    power at wavelength for each spectrum, and takes in a set of
%    receptors, as a SSTReceptor object, and calculates the contrast on
%    each receptor between all pairs of spectra.
%
% Inputs:
%    SPDs         - nWlsxN Matrix of spectral power distributions, where
%                   nWls is the wavelength specification, and N is the
%                   number of spectra to calculate contrasts across
%    SSTReceptors - SSTReceptor object that specifies the receptors on
%                   which to calculate contrasts
%
% Outputs:
%    contrasts    - NxNxR matrix of contrasts in % (one NxN
%                   matrix per receptor type), where contrasts(i,j,R) =
%                   responseDiff(i,j,R) / response(R,i) * 100%
%    response     - RxN matrix of responses of each receptor type to each 
%                   SPD
%    responseDiff - NxNxR matrix of differences in responses (one NxN 
%                   matrix per receptor type), where responseDiff(i,j,R) =
%                   response(R,j) - response(R,i).
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    In the case that only 2 SPDs are passed (e.g., a background SPD and a
%    direction SPD), the outputs are simplified as follows:
%       responseDiff - Rx2 matrix, where the first column is the
%                      response(R,j) - reponse(R,i), and the second column 
%                      the inverse
%       contrasts    - Rx2 matrix, where the first column is the contrast
%                      relative to the first SPD, and the second column is
%                      the contrast relative to the second SPD.
%
%    In the case that only 1 SPD is passed, response is returned as normal,
%    but responseDiff and contrasts are returned as nWlsx1 columnvectors of
%    NaNs.

% History:
%    12/01/17  jv  created based on ComputeAndReportContrastsFromSpds     
%
    
% Calculate receptor response
response = SSTReceptors.T.T_energyNormalized * SPDs;

if size(SPDs,2) <= 1
    responseDiff = NaN(size(SPDs));
    contrasts = NaN(size(SPDs));
else
    % Calculate difference in receptor responses between all SPDs
    temp = reshape(response',[1,size(SPDs,2),size(SSTReceptors.labels,2)]);
    temp = repmat(temp,[size(SPDs,2),1,1]);
    responseDiff = temp - permute(temp,[2 1 3]);
    
    % Squeeze, if only 2 SPDs were passed
    if size(SPDs,2) == 2
        responseDiff = squeeze([responseDiff(1,2,:) responseDiff(2,1,:)])';
        % denominator for contrast is the responses matrix when N = 2
        temp = response;
    else
        % denominator for contrast is a permutation of the temp matrix
        temp = permute(temp,[2 1 3]);
    end
    
    % Calculate contrasts
    contrasts = responseDiff ./ temp * 100;
end
end