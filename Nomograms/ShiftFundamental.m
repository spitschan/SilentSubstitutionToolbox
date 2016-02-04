function out = ShiftFundamental(wls, in, shift);
% out = ShiftFundamental(wls, in, shift);
%
% This function takes a spectral sensitivity (or any function as a a
% function of wavelength) and shifts it according to shift.
%
% 2/4/16    ms      Wrote it.

out = interp1(wls+shift, in, wls, 'spline');