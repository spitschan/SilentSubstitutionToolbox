function out = CalculateLightProperties(S, spd, pupilDiameterMm)

out.radianceWattsPerM2Sr = spd;
out.radianceWattsPerM2Sr(out.radianceWattsPerM2Sr < 0) = 0;
out.radianceWattsPerCm2Sr = (10.^-4)*out.radianceWattsPerM2Sr;
out.radianceQuantaPerCm2SrSec = EnergyToQuanta(S,out.radianceWattsPerCm2Sr);

%% Load CIE functions.   
load T_xyz1931
T_xyz = SplineCmf(S_xyz1931,683*T_xyz1931,S);
out.photopicLuminanceCdM2 = T_xyz(2,:)*out.radianceWattsPerM2Sr;
out.chromaticityXY = T_xyz(1:2,:)*out.radianceWattsPerM2Sr/sum(T_xyz*out.radianceWattsPerM2Sr);

%% Load cone spectral sensitivities
load T_cones_ss2
T_cones = SplineCmf(S_cones_ss2,T_cones_ss2,S);

%% Load in a directly measured sunlight through window
% and off piece of white paper towel on floor for comparison.
% Surely that is safe to look at.
load spd_phillybright
spd_phillybright = SplineSpd(S_phillybright,spd_phillybright,S);
photopicLuminancePhillyBrightCdM2 = T_xyz(2,:)*spd_phillybright;
OLSLratio = out.radianceWattsPerM2Sr./spd_phillybright;

%% Compute irradiance, trolands, etc.
if ~exist('pupilDiameterMm', 'var') | isempty(pupilDiameterMm)
    pupilDiameterMm = 4.7;
    pupilDiameterMm = GetWithDefault('Enter observer pupil diameter in mm',pupilDiameterMm);
end
pupilAreaMm2 = pi*((pupilDiameterMm/2)^2);
eyeLengthMm = 17;
degPerMm = RetinalMMToDegrees(1,eyeLengthMm);
out.irradianceWattsPerUm2 = RadianceToRetIrradiance(out.radianceWattsPerM2Sr,S,pupilAreaMm2,eyeLengthMm);
out.irradianceScotTrolands = RetIrradianceToTrolands(out.irradianceWattsPerUm2, S, 'Scotopic', [], num2str(eyeLengthMm));
out.irradiancePhotTrolands = RetIrradianceToTrolands(out.irradianceWattsPerUm2, S, 'Photopic', [], num2str(eyeLengthMm));
out.irradianceQuantaPerUm2Sec = EnergyToQuanta(S,out.irradianceWattsPerUm2);
out.irradianceWattsPerCm2 = (10.^8)*out.irradianceWattsPerUm2;
out.irradianceQuantaPerCm2Sec = (10.^8)*out.irradianceQuantaPerUm2Sec;
out.irradianceQuantaPerDeg2Sec = (degPerMm^2)*(10.^-2)*out.irradianceQuantaPerCm2Sec;

%% Pupil adjustment factor for Ansi MPE 
mpepupilDiameterMm = 3;
%mpepupilDiameterMm  = GetWithDefault('Enter ANSI 2007 MPE caclulations assumed pupil diameter in mm',mpepupilDiameterMm );
pupilAdjustFactor = (pupilDiameterMm/mpepupilDiameterMm).^2;

%% Get trolands another way.  For scotopic trolands, this just uses scotopic vlambda (in PTB as T_rods)
% and the magic factor of 1700 scotopic lumens per Watt from Wyszecki & Stiles (2cd edition),
% p. 257.  (This is the analog of 683 photopic lumens per Watt.  Then apply formula from
% page 103 of same book.
%
% Same idea for photopic trolands, although we already have luminance in cd/m2 from above so
% we can short cut a little.
%
% The agreement is good to integer scotopic trolands and I'm will to write off the rest
% as round off error.
load T_rods
T_scotopicVlambda = SplineCmf(S_rods,T_rods,S);
out.irradianceScotTrolands_check = pupilAreaMm2*1700*(T_scotopicVlambda*out.radianceWattsPerM2Sr);
out.irradiancePhotTrolands_check = pupilAreaMm2*out.photopicLuminanceCdM2;

%% Get cone coordinates from radiance, and also adjust by pupil area.
% Useful for comparing to light levels produced by monochromatic lights
% in other papers
theLMS = T_cones*out.radianceWattsPerM2Sr;
theLMSTimesPupilArea = pupilAreaMm2*theLMS;

%% Compute irradiance arriving at cornea
%
% According to OSA Handbook of Optics, 2cd Edition, Chaper 24 (vol 2), pp. 24.13-24.15, the
% conversion is (assuming some approximations), irradiance = radiance*stimulusArea/distance^2.
% This is implemented in RadianceAndDistanceAreaToCornIrradiance
stimulusRadiusMm = 6;
stimulusDistanceMm = 25;
stimulusRadiusM = stimulusRadiusMm/1000;
stimulusAreaM2 = pi*(stimulusRadiusM^2);
stimulusDistanceM = stimulusDistanceMm/1000;
stimulusRadiusDeg = rad2deg(stimulusRadiusMm/stimulusDistanceMm);
stimulusAreaDegrees2 = pi*(stimulusRadiusDeg^2);
out.cornealIrradianceWattsPerM2 = RadianceAndDistanceAreaToCornIrradiance(out.radianceWattsPerM2Sr,stimulusDistanceM,stimulusAreaM2);
out.cornealIrradianceWattsPerCm2 = (10.^-4)*out.cornealIrradianceWattsPerM2;
out.cornealIrradianceQuantaPerCm2Sec = EnergyToQuanta(S,out.cornealIrradianceWattsPerCm2);
    
%% Let's convert to melanopic units, as well as to equivalent stimulation at specific wavelengths
%
% We have retinal and corneal spectral irradiance
%  S
%  irradianceQuantaPerCm2Sec
%  cornealIrradianceWattsPerCm2
melanopsinAssumedFieldSizeDeg = 10;
melanopsonAssumeAgeYears = 32;
[~,T_melanopsinQuantal] = GetHumanPhotoreceptorSS(S, {'Melanopsin'},melanopsinAssumedFieldSizeDeg,melanopsonAssumeAgeYears,pupilDiameterMm,[],[],[],[]);
T_melanopsinQuantal = T_melanopsinQuantal/max(T_melanopsinQuantal(1,:));
out.melIrradianceQuantaPerCm2Sec = T_melanopsinQuantal*out.irradianceQuantaPerCm2Sec;
out.melCornealIrradianceQuantaPerCm2Sec = T_melanopsinQuantal*out.cornealIrradianceQuantaPerCm2Sec;
fprintf('\n');
fprintf('  * Melanopic retinal irradiance %0.1f log10 melanopic quanta/[cm2-sec]\n',log10(out.melIrradianceQuantaPerCm2Sec));
fprintf('  * Melanopic corneal irradiance %0.1f log10 melanopic quanta/[cm2-sec]\n',log10(out.melCornealIrradianceQuantaPerCm2Sec));

% Convert Dacey reference retinal irradiances to melanopic units for
% comparison.  Dacey et al. 2005 gives 11-15 log quanta/[cm2-sec] as
% the range over which they measured sustained melanopsin responses.  It
% isn't clear that things were saturating at the high end, though.
index = find(SToWls(S) == 470);
if (isempty(index))
    error('Oops.  Need to find closest wavelength match as exact is not in sampled wls');
end
tempSpd = zeros(size(out.irradianceQuantaPerCm2Sec));
tempSpd(index) = 10^11;
fprintf('  * Dacey 2005 low, 11 log10 quanta/[cm2-sec] at 470 nm, is %0.1f is log10 melanopic quanta/[cm2-sec]\n',log10(T_melanopsinQuantal*tempSpd));
tempSpd = zeros(size(out.irradianceQuantaPerCm2Sec));
tempSpd(index) = 10^15;
fprintf('  * Dacey 2005 high, 15 log10 quanta/[cm2-sec] at 470 nm, is %0.1f is log10 melanopic quanta/[cm2-sec]\n',log10(T_melanopsinQuantal*tempSpd));

% Lucas (2012) says that a mouse retina has an area of 18 mm2 and that the
% mouse pupil varies between 9 to 0.1 mm2 in are.  For a fully
% dialated pupil and a full field stimulus, this givea a correction between
% corneal and retinal irradiance is that retinal is
% corneal*(pupilArea/retinalArea).  So for fully dialated pupil, retinal
% illuminance is about about half of corneal.
%
% As a check of this formula, we could ask whether our retinal and corneal
% irradiances are related in this way.  We have our pupil area, and we are
% assuming that the stimulus is 6 mm in radius at 25 mm from the eye.
% Given an eye length of 17 mm, which is about right, we can compute the
% retina radius of the stimulus as 6*17/25 as about 4 mm, and its area as
% about 50 mm2.  Given a pupil diameter of 4.7 mm, the pupil area is 17.34.
% So we should have retinal illuminance is 17.34/50*corneal illuminance, or
% 0.34 * corneal illuminance
%   sum(irradianceQuantaPerCm2Sec)
%   0.34*sum(cornealIrradianceQuantaPerCm2Sec)
% These numbers agree within 3 percent, which we're taking to be good
% enough for now.
%
% Lucas gives something like 11 log quanta/[cm2-sec] in mice as the low end of the
% melanopsin operating range for light between 480 and 500 nm.  We convert
% to retinal illuminance by multiplying by 0.5, and then to melanopic units
tempSpd = zeros(size(out.irradianceQuantaPerCm2Sec));
index = find(SToWls(S) == 480);
tempSpd(index) = (0.5*10^11)/3;
index = find(SToWls(S) == 490);
tempSpd(index) = (0.5*10^11)/3;
index = find(SToWls(S) == 500);
tempSpd(index) = (0.5*10^11)/3;
fprintf('  * Lucas low is %0.1f log10 melanopic quanta/[cm2-sec]\n',log10(T_melanopsinQuantal*tempSpd));

% At the high end, Lucas gives estimate of 10^15 quanta/[cm2-sec] in same
% wl range as melanopsin saturation, at the retina.
index = find(SToWls(S) == 480);
tempSpd(index) = (10^15)/3;
index = find(SToWls(S) == 490);
tempSpd(index) = (10^15)/3;
index = find(SToWls(S) == 500);
tempSpd(index) = (10^15)/3;
fprintf('  * Lucas high is %0.1f log10 melanopic quanta/[cm2-sec]\n',log10(T_melanopsinQuantal*tempSpd));

%% Get MPE from as a function of wavelength.  For each wavelength,
% take minimum radiance over specified sizes and durations.

% Specify what parameters to test
minLogSize = -1; maxLogSize = 2;
minLogDuration = -1; maxLogDuration = 4;
minLogYRad = -3; maxLogYRad = 2;
minLogYIrrad = -5; maxLogYIrrad = 0;
minLogYIntRad = 0; maxLogYIntRad = 3;
minLogYRadExp = -4; maxLogYRadExp = -1;
measuredWls = SToWls(S);
index = find(measuredWls >= 400);
stimulusWavelengthsNm = measuredWls(index);
stimulusSizesDeg = logspace(minLogSize,maxLogSize,5);
stimulusDurationsSec = logspace(minLogDuration,maxLogDuration,5);
%fprintf('Computing MPE over wavelengths from %0.1f to %0.1f deg\n',min(stimulusWavelengthsNm),max(stimulusWavelengthsNm));
clear MPELimitIntegratedRadiance_JoulesPerCm2Sr MPELimitRadiance_WattsPerCm2Sr MPELimitCornealIrradiance_WattsPerCm2 MPELimitCornealRadiantExposure_JoulesPerCm2
for w = 1:length(stimulusWavelengthsNm)
    stimulusWavelengthNm = stimulusWavelengthsNm(w);
    if (rem(w,10) == 0)  
        %fprintf('\tComputing minimum MPE for wavelength %d nm\n',stimulusWavelengthNm);
    end
    MPELimitIntegratedRadiance_JoulesPerCm2Sr(w) = Inf;
    MPELimitRadiance_WattsPerCm2Sr(w) = Inf;
    MPELimitCornealIrradiance_WattsPerCm2(w) = Inf;
    MPELimitCornealRadiantExposure_JoulesPerCm2(w) = Inf;
    for s = 1:length(stimulusSizesDeg)
        stimulusSizeDeg = stimulusSizesDeg(s);
        stimulusSizeMrad = DegToMrad(stimulusSizeDeg);
        for t = 1:length(stimulusDurationsSec)
            stimulusDurationSec = stimulusDurationsSec(t);
            
            % Compute MPE.  We don't understand how the cone limit computations fit in with
            % the standard, or not.  So, we run it both ways and take the lower limit returned.
            [temp1, temp2, temp3, temp4] = ...
                AnsiZ136MPEComputeExtendedSourceLimit(stimulusDurationSec,stimulusSizeDeg,stimulusWavelengthNm,0);
            [temp5, temp6, temp7, temp8] = ...
                AnsiZ136MPEComputeExtendedSourceLimit(stimulusDurationSec,stimulusSizeDeg,stimulusWavelengthNm,1);
            if (temp5 < temp1)
                temp1 = temp5;
            end
            if (temp6 < temp2)
                temp2 = temp6;
            end
            if (temp7 < temp3);
                temp3 = temp7;
            end
            if (temp8 < temp4)
                temp4 = temp8;
            end
            clear temp5 temp6 temp7 temp8
            
            % Store minimum at each wavelength.
            if (temp1 < MPELimitIntegratedRadiance_JoulesPerCm2Sr(w))
                MPELimitIntegratedRadiance_JoulesPerCm2Sr(w) = temp1;
            end
            if (temp2 < MPELimitRadiance_WattsPerCm2Sr(w))
                MPELimitRadiance_WattsPerCm2Sr(w) = temp2;
            end
            if (temp3 < MPELimitCornealIrradiance_WattsPerCm2(w))
                MPELimitCornealIrradiance_WattsPerCm2(w) = temp3;
            end
            if (temp4 < MPELimitCornealRadiantExposure_JoulesPerCm2(w))
                MPELimitCornealRadiantExposure_JoulesPerCm2(w) = temp4;
            end
        end
    end
end

%% Find how much total radiance we could tolerate if all our power was at the 
% wavelength with minimum MPE.
minMPERadiance = min(MPELimitRadiance_WattsPerCm2Sr(:));
fprintf('\n');
fprintf('  * Compute ANSI 2007 MPE as a function of wavelength.  For each wavelength, took minimum over size and duration\n');
fprintf('    * Size range: %0.1f to %0.1f degrees\n',min(stimulusSizesDeg),max(stimulusSizesDeg));
fprintf('    * Duration range: %0.1f to %0.1f seconds\n',min(stimulusDurationsSec),max(stimulusDurationsSec));
fprintf('  * Minimum ANSI MPE value over wavelengths: radiance %0.1f log W/[cm2-sr]\n',log10(minMPERadiance));
fprintf('    * Compare with total stimulus radiance %0.1f log  W/[cm2-sr]\n',log10(sum(out.radianceWattsPerCm2Sr)));
fprintf('    * Compare with total pupil adjusted radiance %0.1f log  W/[cm2-sr]\n',log10(sum(out.radianceWattsPerCm2Sr))+log10(pupilAdjustFactor));
fprintf('    * Pupil adjustment assumes observer pupil diameter of %0.1f mm, MPE standard diameter of %0.1f mm\n',pupilDiameterMm,mpepupilDiameterMm);

%% Sum over wavelength of power divided by MPE
% Could put this back in, but would have to think
% a bit harder about wavelength spacing adjustment.
%
% index = find(stimulusWavelengthsNm >= 400);
% deltaMeasuredWls = measuredWls(2)-measuredWls(1);
% deltaMPEWls = stimulusWavelengthsNm(2)-stimulusWavelengthsNm(1);
% MPERatioSum = 0;
% for i = 1:length(stimulusWavelengthsNm)
%     index = find(measuredWls == stimulusWavelengthsNm(i));
%     MPERatioSum = MPERatioSum + radianceWattsPerCm2Sr(index)/MPELimitRadiance_WattsPerCm2Sr(i);
% end
% fprintf('MPERatioSum = %0.4f\n',MPERatioSum*deltaMPEWls/deltaMeasuredWls);

%% Now compare to the ISO Standard
stimulusDurationForISOMPESecs = 60*60;
[IsOverLimit,ISO2007MPEStruct] = ISO2007MPECheckType1ContinuousRadiance(S,out.radianceWattsPerM2Sr,stimulusDurationForISOMPESecs,stimulusAreaDegrees2,eyeLengthMm);
fprintf('  * ISO MPE Analysis\n');
ISO2007MPEPrintAnalysis(IsOverLimit,ISO2007MPEStruct);
fprintf('  * Assumed duration seconds %0.1f, hours %0.1f\n',stimulusDurationForISOMPESecs,stimulusDurationForISOMPESecs/3600);
