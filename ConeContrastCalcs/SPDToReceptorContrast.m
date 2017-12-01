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
%    contrasts    - NxNxR matrix of contrasts (one NxN matrix per receptor 
%                   type)
%    response     - RxN matrix of responses of each receptor type to each 
%                   PSD
%    responseDiff - NxNxR matrix of differences in responses (one NxN 
%                   matrix per receptor type)
%
% Optional key/value pairs:
%   None.
%
% Notes:
%
% See also:

% History:
%    12/01/17  jv  created. based on ComputerAndReportContrastsFromSpds     
%
     
    % Calculate receptor response
    response = SSTReceptors.T.T_energyNormalized * SPDs;

    % Calculate difference in receptor responses between all SPDs
    temp = reshape(response',[1,size(SPDs,2),size(SSTReceptors.labels,2)]);
    temp = repmat(temp,[size(SPDs,2),1,1]);
    responseDiff = temp - permute(temp,[2 1 3]);
    
    % Calculate contrasts
    contrasts = responseDiff ./ permute(temp,[2 1 3]) * 100;
end