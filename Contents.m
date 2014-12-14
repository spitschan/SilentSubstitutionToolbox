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
%   GetHumanPhotoreceptorSS - Wrapper around PTB functions to obtain
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
%   Matlab Optimization Toolbox
%   Matlab Statistics Toolbox
%
% License
% -------
%   Except as noted above for files in the External directory, the
%   SilentSubstitutionToolbox is covered by:
% 
%  The MIT License (MIT)
% 
%   Copyright (c) Manuel Spitschan, Geoffrey Aguirre, David Brainard and the
%   Trustees of the University of Pennsylvania.
% 
%   Permission is hereby granted, free of charge, to any person obtaining a
%   copy of this software and associated documentation files (the
%   "Software"), to deal in the Software without restriction, including
%   without limitation the rights to use, copy, modify, merge, publish,
%   distribute, sublicense, and/or sell copies of the Software, and to permit
%   persons to whom the Software is furnished to do so, subject to the
%   following conditions:
% 
%   The above copyright notice and this permission notice shall be included
%   in all copies or substantial portions of the Software.
% 
%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
%   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
%   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
%   NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
%   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
%   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
%   USE OR OTHER DEALINGS IN THE SOFTWARE.
% 
%   If you use this code in support of work in a published paper, please cite us.
%   We are working on publishing a paper that includes a description of the logic used in this code as
%   well as examples of its use.  Once we have published our paper, that will be the appropriate
%   work to cite.  For now, please use:
%     Sptitschan, M., Aguirre, G.K., & Brainard D.H. (2015), The silent substitution toolbox,
%     https://github.com/spitschan/SilentSubstitutionToolbox.
