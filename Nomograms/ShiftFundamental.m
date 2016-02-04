function out = ShiftFundamental(wls, in, shift);
out = interp1(wls+shift, in, wls, 'spline');