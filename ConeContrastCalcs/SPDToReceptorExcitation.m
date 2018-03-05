function excitation = SPDToReceptorExcitation(SPD,receptors)
% Calculates receptor excitations from spectral power distributions
%
% Syntax:
%   excitation = SPDToReceptorExcitation(SPD,T_receptors);
%   excitation = SPDToReceptorExcitation(SPD,SSTReceptors);
%
% Description:
%    Takes in one or more spectral power distributions, and a set of
%    receptors (sensitivities), and returns the excitation of each
%    receptor type to each SPD.
%
% Inputs:
%    SPD        - nWlsxN matrix of N spectral power distributions at each
%                 of the nWls wavelength bands.
%    receptors  - either:
%                 - RxnWls matrix (T_receptors) of R receptor sensitivities
%                   sampled at nWls wavelength bands
%                 - SSTReceptor-object, in which case the
%                   T.T_energyNormalized matrix will be used
%
% Outputs:
%    excitation - RxN matrix of excitations of the R receptors for each of
%                 the N SPDs.
%
% Optional key/value pairs:
%    None.
%
% Examples are provided in the source code.
%
% See also:
%    ReceptorExcitationToReceptorContrast, SPDToReceptorContrast

% History:
%    03/02/18  jv  wrote it, extracted from
%                  ComputeAndReportContrastsFromSPDs

% Examples:
%{
    receptors = SSTReceptorHuman('verbosity','low'); % default human receptors

    % use monochromatic SPDs, with peak at peak wavelength for a receptor:
    SPDs = receptors.T.T_energyNormalized' .* (receptors.T.T_energyNormalized' == max(receptors.T.T_energyNormalized'));
    
    % input argument can be an SSTReceptor object:
    excitations = SPDToReceptorExcitation(SPDs,receptors);

    % input argument can also be a T_receptors matrix
    excitations = SPDToReceptorExcitation(SPDs,receptors.T.T_energyNormalized);
%}

%% Input validation
parser = inputParser;
parser.addRequired('SPD',@isnumeric);
parser.addRequired('receptors',@(x) isnumeric(x) || isa(x,'SSTReceptor'));
parser.parse(SPD,receptors);

%% If SSTReceptor, extract T_receptors
if isa(receptors,'SSTReceptor')
    receptors = receptors.T.T_energyNormalized;
end

%% Calculate excitation
excitation = receptors * SPD;

end