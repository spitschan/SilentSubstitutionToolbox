% SilentSubstitutionToolbox - A toolbox for silent substitution
% 
% SilentSubstitutionFunctions
%
% ContrastSplatter functions
% -----------------------------
% 
%   CalculateSplatter   - Calculates contrast splatter maps
%   SaveSplatter        - Saves contrast splatter maps
%   SaveSplatterConfidenceBounds    
%                       - Saves splatter statistics in CIs
%   PlotSplatter        - Plots contrast splatter maps
%   GetHumanPhotopigmentSS
%                       - Wrapper around PTB functions to obtain standard
%                         spectral sensitivities
%   GetChronologicalAgeSDFromLensSD     
%                       - Converts age standard deviation into lens density
%                         standard deviation
%   GetLambdaMaxEstimateSD
%                       - Returns lambda-max standard deviation for a given
%                         photopigment
%   GetConeFractionBleachedFromSpectrum
%                       - Calculates proportion of pigment bleached for a
%                         given background spectrum
%   GetHemoglobinTransmittance
%                       - Returns hemoglobin transmittance
%   GetHemoglobin       - Wrapper function to get hemoglobin
%                         extinction/absorptivity/absorption functions
%
% Data
% ----
%
%   ContrastSplatter/data/XuPokornySmith1997_Fig4A.csv
%                       - Contains digitized chronological age vs. lens
%                         density plot, for use in
%                         GetChronologicalAgeSDFromLensSD
%   ContrastSplatter/data/spd_background.mat
%                       - Contains a sample background spectrum
%   ContrastSplatter/data/spd_melIsolatingSpd.mat
%                       - Contains a sample melanopsin-isolating spectrum
%
% Third-party functions
% ---------------------
%
%   error_ellipse       - Plots an error ellipse. No known license.
%                         http://www.mathworks.com/matlabcentral/fileexchange/4705-error-ellipse
%   lbmap               - Color maps. BSD license.
%                         http://www.mathworks.com/matlabcentral/fileexchange/17555-light-bartlein-color-maps/content/lbmap.m
% 
%
% Demo
% ----
%
%   ContrastSplatterDemo - Demos the splatter calculations
%
%
% Requirements
% ------------
%   Psychtoolbox - http://psychtoolbox.org/
%       SToWls
%       SplineSrf
%       SplineCmf
%       QuantaToEnergy
%       EnergyToQuanta
%       ComputeCIEConeFundamentals
%       DefaultPhotoreceptors
%       FillInPhotoreceptors
%       PhotonAbsortionRate
%       RadianceToRetIrradiance
%       RetIrradianceToTrolands
%       RetinalMMToDegrees
%
%   Statistics Toolbox
%       mvnpdf