function B_primary = LEDPrimaries(S, peakWls, FWHM, maxPower)
%LEDPRIMARIES constructs B_primary for any N LEDs
%
%   S = spectral resolution to define primaries, in standard PTB format ([start delta N]) 
%
%   peakWls = vector of wavelength (in nm) of peak output for each LED
%
%   FWHM = vector of full-width half maximum for each LED
%
%   maxPower = vector of power at peak output for each LED
%
if (S ~= MakeItS(S))
    warning('S not in standard PTB format ( [start delta N] )! Casting...');
end

assert(numel(peakWls) == numel(FWHM) && numel(peakWls) == numel(maxPower),...
    'Inconsistent number of LEDs specified in input arguments');

wls = SToWls(S);

spd = normpdf(wls,peakWls,FWHM);
B_primary = spd ./ max(spd) .* maxPower;


end

