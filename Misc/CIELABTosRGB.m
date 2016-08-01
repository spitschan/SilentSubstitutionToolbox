function [XYZs bgXYZ RGBs bgRGB] = CIELABTosRGB(angles, theLstar, theSaturation)
% function [XYZs bgXYZ, RGBs, bgRGB] = generateColorsCIELAB(angles, theLstar, theSaturation);
%
% 12/4/2012     spitschan      Wrote it.
% 12/8/2012     spitschan      Clean up, added to SVN
% 12/9/2012     spitschan      Specific version for L*a*b*
% 7/2/2013      spitschan      Unearthed.
%
% Input: angles - Angles in degree around color wheel
%        theLstar - L* coordinate in CIELAB space
%        theSaturation - saturation value (i.e. radial eccentricity)
%
% Output: XYZ coordinates for passed colors and background (half-on)
%         sRGB coordinates for passed colors and background (half-on)
%
% Demo:
%
% theAngles = 1:360;
% theLstar = 25;
% theSaturation = 10;
%
% [theXYZ, ~, theRGBs] = CIELABTosRGB(theAngles, theLstar, theSaturation);
% imagesc(permute(theRGBs, [3, 2, 1])/255); % Note the scaling 

% Load in the XYZ fundamentals
load T_xyz1931
T_xyz = 683*T_xyz1931;
S_xyz = S_xyz1931;

% Find out what brightest thing we can display is (here, in sRGB), in XYZ
maxXYZ = SRGBPrimaryToXYZ([1 1 1]');

% Compute L*a*b* coordinates
Lab = [ones(1, length(angles))*theLstar ; theSaturation*sind(angles) ; theSaturation*cosd(angles)]';
bgLab = [theLstar/2 0 0]';

% For sanity, check the delta-E metric
% for i = 2:360;
%    deltaE(i) = ComputeDE(Lab(i, :)',Lab(i-1, :)')
% end

% Use D65 as white/gray point
load spd_D65
load T_xyz1931
XYZ_D65 = T_xyz1931*spd_D65;

% Convert L*a*b* to XYZ tristimulus coordinates
XYZs = LabToXYZ(Lab',XYZ_D65);
bgXYZ = LabToXYZ(bgLab,XYZ_D65);

RGBs =  XYZToSRGBPrimary(XYZs);
bgRGB = XYZToSRGBPrimary(bgXYZ);
