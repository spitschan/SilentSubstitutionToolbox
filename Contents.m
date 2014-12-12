% SilentSubstitutionToolbox - A toolbox for silent substitution calculations.
%
% Manuel Spitschan, Geoffrey Aguirre, David Brainard (c) 2014.
%
% If you make use of this software in a publication, please cite
%   Spitschan, M.,  Aguirre, G.K. & Brainard D.H. (2014), Manuscript in
%   preparation.
%
% Demos
% -------------------------
%   ReceptorIsolateDemo  - Demos the ReceptorIsolate, PhotoreceptorSensitivity,
%                          and HemoglobinTransmivitity function
%   ContrastSplatterDemo - Demos the contrast splatter calculations
%
% ReceptorIsolate Functions
% -------------------------
%   ReceptorIsolate     - Function that finds isolating device primary
%                         settings
%
% PhotoreceptorSensitivity Functions
% -----------------------------
%   GetHumanPhotopigmentSS - Wrapper around PTB functions to obtain
%                         photoreceptor spectral sensitivities
%   GetChronologicalAgeSDFromLensSD - Converts age standard deviation into lens density
%                         standard deviation
%   GetLambdaMaxEstimateSD - Returns lambda-max standard deviation for a given
%                         photopigment
%   GetConeFractionBleachedFromSpectrum - Calculates proportion of pigment bleached for a
%                         given background spectrum
%
% HemoglobinTransmisivity Functions
% -----------------------------
%   GetHemoglobinTransmittance
%                       - Returns hemoglobin transmittance
%   GetHemoglobin       - Wrapper function to get hemoglobin
%                         extinction/absorptivity/absorption functions
%
% ContrastSplatter Functions
% -----------------------------
%   CalculateSplatter   - Calculates contrast splatter maps
%   SaveSplatter        - Saves contrast splatter maps SaveSplatterConfidenceBounds    
%                       - Saves splatter statistics in CIs
%   PlotSplatter        - Plots contrast splatter maps
%
% Data
% ----
%   PhotoreceptorSensitivities /xRawData/XuPokornySmith1997_Fig4A.csv
%                       - Contains digitized chronological age vs. lens
%                         density plot, for use in GetChronologicalAgeSDFromLensSD
%   ContrastSplatter/ContrastSplatterDemoData/spd_background.mat
%                       - Contains a sample background spectrum
%   ContrastSplatter/ContrastSplatterDemoData/spd_melIsolatingSpd.mat
%                       - Contains a sample melanopsin-isolating spectrum
%   ContrastSplatter/ContrastSplatterDemoOutput
%                       - Data and plots produced by ContrastSplatterDemo
%   ReceptorIsolate/ReceptorIsolateDemoOutput
%                       - Data and plots produced by ReceptorIsolateDemo
%
% External Functions
% ---------------------
%     error_ellipse       - Plots an error ellipse. No known license.
%                           http://www.mathworks.com/matlabcentral/fileexchange/4705-error-ellipse
%     lbmap               - Color maps. BSD license.
%                           http://www.mathworks.com/matlabcentral/fileexchange/17555-light-bartlein-color-maps/content/lbmap.m
%
% Requirements
% ------------
%   Open-source Psychtoolbox - http://psychtoolbox.org/
%   Matlab Statistics Toolbox
%
% License
% -------
%   Licensed under the MIT Open-Source License.