function indDiffParams = DefaultIndDiffParams
% DefaultIndDiffParams
%
% Usage:
%     indDiffParams = DefaultIndDiffParams()
%
% Description:
%     This returns the default individual difference parameters (no change
%     relative to the nominal ones).
%
% Input:
%     None
%
% Output:
%     indDiffParams - struct containing the default individual difference
%                     parameters in its field. These are:
%                           indDiffParams.dlens - lens density
%                           indDiffParams.dmac - macular pigment density
%                           indDiffParams.dphotopigment - 1x3 vector of
%                                                         photopigment
%                                                         densities for LMS
%                                                         cones
%                           indDiffParams.lambdaMaxShift - 1x3 vector of
%                                                          lambda max
%                                                          shifts for LMS
%                                                          cones
%                           indDiffParams.shiftType - shift type for
%                                                     lambda-max adjustment
%
% Optional key/value pairs:
%     None.

% 9/8/17  ms  Added header comments.
indDiffParams.dlens = 0;
indDiffParams.dmac = 0;
indDiffParams.dphotopigment = [0 0 0];
indDiffParams.lambdaMaxShift = [0 0 0];
indDiffParams.shiftType = 'linear';
