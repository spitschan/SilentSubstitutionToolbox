function [postreceptoralContrasts, postreceptoralStrings] = ComputePostreceptoralContrastsFromLMSContrasts(LMSContrasts,varargin)
% ComputePostreceptoralContrastsFromLMSContrasts  Conver LMS cone contrasts to post-receptoral contrasts
%
% Syntax:
%     [postreceptoralContrasts, postreceptoralStrings] = ComputePostreceptoralContrastsFromLMSContrasts(LMSContrasts)
%
% Description:
%     Compute standard post-receptoral contrasts from LMS cone contrasts.
%
%     T LMS cone contrasts to a post-receptoral opponent representation
%     assuming mechanism sensitivities to cone contrast for luminance
%     (isochromatic), red-green, and blue-yellow mechanisms of [0.5 0.5 0],
%     [0.5 -0.5 0], and [-0.5 -0.5 1] respectively. This transformation
%     corresponds to the DKL opponent color space representation when the
%     background produces equal excitations in the L, M and S cones, for
%     the case where the L and M cone spectral sensitivities are scaled so
%     that they sum to produce the luminous efficiency curve. We regard
%     this as a a reasonable choice of reference conditions to define the
%     transformation, as it leads to intuitively straightforward mechanism
%     properties.
%
%     This choice of conversion has the advantage that it seems about as
%     reasonable as anything else, matches what we have been doing in the
%     Brainard lab recently, and is explainable (as above). The property
%     that we expect iso-response contours to better line up with the axes
%     of an opponent space makes this a reasonable space for looking at
%     splatter, even if it is more complicated than just examining LMS
%     contrast.

%     When the input is expressed as cone contrast, the result does not
%     depend on the scaling chosen for the cone fundamentals, so agreement
%     for one reasonable choice of scaling is sufficient.
%
%     For other backgrounds, this transformation will describe the opponent
%     mechanism responses to the extent that those responses are the same
%     for modulations seen against different backgrounds when their LMS
%     cone contrasts are matched.
%
% Input:
%    LMSContrasts - 3 by N matrix of LMS cone contrasts
%
% Output:
%    postreceptoralContrasts - 3 by N matrix of post-receptoral contrasts.
%    postreceptoralStrings   - String labels for the postreceptoral combinations
%
% Optional key/value pairs:
%   None.
%
% See also: 

% 09/08/17  dhb   Wrote this.

%% Set transformation matrix
%
% This is expressed in terms of the contrast sensitivities
% of the mechanisms.
xformMatrix = [0.5 0.5 0; 0.5 -0.5 0; -0.5 -0.5 1];

%% Do the transform
postreceptoralContrasts = xformMatrix*LMSContrasts;

%% Set the strings
postreceptoralStrings = {'L+M+S', 'L-M', 'S-(L+M)'};

%% Another way to do the xform, for checking
%
% The columns of this matrix are the directions that 
% isolate each of the three mechanisms.
postreceptoralContrastsCheck = [1 1 1 ; 1 -1 0; 0 0 1]' \ LMSContrasts;
checkTolerance = 1e-8;
if (any(abs(postreceptoralContrasts(:) - postreceptoralContrastsCheck(:)) > checkTolerance))
    error('Two ways of doing the computation do not match');
end

%% To check, try this:
%
% ComputePostreceptoralContrastsFromLMSContrasts([1 1 1 ; 1 -1 0 ; 0 0 1 ; 0.5 0.2 0.3]')
%
% Which should produce this:
%
% ans =
% 
%             1            0            0         0.35
%             0            1            0         0.15
%             0            0            1        -0.05


