function [contrasts, excitation, excitationDiff] = SPDToReceptorContrast(SPDs,receptors)
% Calculates contrast on photoreceptors between SPDs
%
% Syntax:
%   contrasts = SPDtoReceptorContrast(SPDs,T_Receptors)
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
%    receptors    - either:
%                   - RxnWls matrix (T_receptors) of R receptor
%                   sensitivities sampled at nWls wavelength bands
%                   - SSTReceptor-object, in which case the
%                     T.T_energyNormalized matrix will be used
%
% Outputs:
%    contrasts    - NxNxR matrix of contrasts (one NxN
%                   matrix per receptor type), where contrasts(i,j,R) =
%                   excitationDiff(i,j,R) / excitation(R,i)
%    excitation     - RxN matrix of excitations of each receptor type to each
%                   SPD
%    excitationDiff - NxNxR matrix of differences in excitations (one NxN
%                   matrix per receptor type), where excitationDiff(i,j,R) =
%                   excitation(R,j) - excitation(R,i).
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    In the case that only 2 SPDs are passed (e.g., a background SPD and a
%    direction SPD), the outputs are simplified as follows:
%       excitationDiff - Rx2 matrix, where the first column is the
%                      excitation(R,j) - reponse(R,i), and the second column
%                      the inverse
%       contrasts    - Rx2 matrix, where the first column is the contrast
%                      relative to the first SPD, and the second column is
%                      the contrast relative to the second SPD.
%
%    In the case that only 1 SPD is passed, excitation is returned as normal,
%    but excitationDiff and contrasts are returned as nWlsx1 columnvectors of
%    NaNs.

% History:
%    12/01/17  jv  created based on ComputeAndReportContrastsFromSpds
%    03/02/18  jv  extracted excitation calculation to
%                  SPDToReceptorReceptorExcitation, and the contrast
%                  calculations to ReceptorExcitationToContrast.

excitation = SPDToReceptorExcitation(SPDs, receptors);
[contrasts, excitationDiff] = ReceptorExcitationToReceptorContrast(excitation);
end