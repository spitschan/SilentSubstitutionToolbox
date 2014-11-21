% SplatterToolbox - A toolbox to calculate uncertainty in receptor
% isolation
% 
% Main functions
% --------------
% 
%   CalculateSplatter   - Calculates contrast splatter maps
%   SaveSplatter        - Saves contrast splatter maps
%   SaveSplatterConfidenceBounds    
%                       - Saves splatter statistics in CIs
%
% Plotting functions
% ------------------
%
%   PlotSplatter        - Plots contrast splatter maps
%
% Auxiliary functions
% -------------------
%
%   GetChronologicalAgeSDFromLensSD     
%                       - Converts age standard deviation into lens density
%                         standard deviation
%   GetLambdaMaxEstimateSD
%                       - Returns lambda-max standard deviation for a given
%                         photopigment
%   GetConeFractionBleachedFromSpectrum
%                       - Calculates proportion of pigment bleached for a
%                         given background spectrum
%
% Data
% ----
%
%   data/XuPokornySmith1997_Fig4A.csv
%                       - Contains digitized chronological age vs. lens
%                         density plot, for use in
%                         GetChronologicalAgeSDFromLensSD
%   data/spd_background.mat
%                       - Contains a sample background spectrum
%   data/spd_melIsolatingSpd.mat
%                       - Contains a sample melanopsin-isolating spectrum
%
% Demo
% ----
%
%   SplatterToolboxDemo - Demos the splatter calculations
%
%
% Requirements
% ------------
%   Psychtoolbox - http://psychtoolbox.org/
%   error_ellipse - http://www.mathworks.com/matlabcentral/fileexchange/4705-error-ellipse