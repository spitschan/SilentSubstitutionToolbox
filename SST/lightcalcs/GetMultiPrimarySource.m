function B_primary = GetMultiPrimarySource(SOrWls, whichSource, rectifyYesOrNo)
% GetMultiPrimarySource(SOrWls, whichSource)
%
% Usage:
%     B_primary = GetMultiPrimarySource(SOrWls, whichSource);
%
% Description:
%     This function returns the spectral bases for a set of multi-primary
%     light sources.
%
% Input:
%     SOrWls - Wavelength column vector or S triplet to specify wavelength
%              sampling.
%     whichSource - The name of the light source to be returned. Possible
%                   options are:
%                       - 'SmithsonFourPrimary'
%                       - 'iQLED'
%                       - 'LEDCube'
%                       - 'LIFXColour1000'
%                       - 'OneLight'
%     rectifyYesOrNo - Logical determining whether to rectify negative 
%                      values
%
% Output:
%     B_primary bases of the light source
%
% Optional key/value pairs:
%     None

% 9/30/17  ms  Added header comments.

% In case of doubt, turn into 'S' representation
S = MakeItS(SOrWls);

% Search for the projects folder
basePath = tbGetPref('projectRoot', []);

% See if the folder is there
if ~isdir(fullfile(basePath, 'LightAndReceptorCalculations'))
    error('Cannot find private LightAndReceptorCalculations repo');
end

% Load the light source
tmp = load(fullfile(basePath, 'LightAndReceptorCalculations', 'xLightSources', whichSource, ['spd_' whichSource '.mat']));

% Rename the variables
eval(['wls_orig = tmp.wls_' whichSource ';']);
eval(['spd = tmp.spd_' whichSource ';']);

% Spline it to specified wavelength basis
B_primary = SplineSpd(wls_orig, spd, S);

% Rectify if required
if rectifyYesOrNo
   B_primary(B_primary < 0) = 0;
end