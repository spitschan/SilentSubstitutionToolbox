% LightSafetyExample
%
% Example calculation comparing a broadband spectrum to light safety limits.
%
% This relies on underlying routines in the Psychophysics toolbox.
%
% We think the relevant standard for braodband light is:
%   Ansi ISO 15004-2 (2007). Ophthalmic instruments - Fundamental requirements and test methods -
%   Part 2: Light hazard protection.
% We rely on the Psychtoolbox implementation of this standard. Please see
% extensive notes there.
%
% There is also an ANSI standard for laser light.  This is a little harder
% to connect to broadband light, but gives us another check.  Again, we
% rely on the Psychtoolbox implementation.  The PTB currently implements
% the 2007 version of this standard.  There is now an updated version which
% is more conservative and for which we do not currently have an
% implementation.
%
% This routine compares an example broadband spectrum to the ISO standard.
% It prompts for the pupil size of the observer and for that implicit in
% the standard (since the latter is not completely clear) and adjusts
% retinal illuminance according to these two input sizes.
%
% It also compares the example spectrum to the light measured off a white
% piece of paper in sunlight in Philadelphia, as a simple check.
%
% *****************************************************************
%
% IMPORTANT
%
%   Individuals using these routines must accept full responsibility 
%   for light exposures they implement. We recommend that values computed
%   with these routines be carefully checked against independent calculations.
%   We have done our best to follow the standards, but they are very complex and
%   there may be errors.  If you find an error in any of our calculations,
%   please let us know about it.
%
% *****************************************************************
%
% 1/28/15  dhb, ms  Extracted example from the larger workhorse program used in our lab.

%% Clear and close
clear; close all

%% Load in an example spectrum.
% This has units of radiance in Watts/[sr-M2-wlband]
load(fullfile('LightSafetyExampleData','spd_lightsafetyexample.mat'));
S = S_lightsafetyexample;
radianceWattsPerM2Sr = spd_lightsafetyexample;
clear spd_lightsafetyexample S_lightsafetyexample

%% Unit coversion
radianceWattsPerCm2Sr = (10.^-4)*radianceWattsPerM2Sr;
radianceQuantaPerCm2SrSec = EnergyToQuanta(S,radianceWattsPerCm2Sr);

%% Load CIE functions.   
load T_xyz1931
T_xyz = SplineCmf(S_xyz1931,683*T_xyz1931,S);
photopicLuminanceCdM2 = T_xyz(2,:)*radianceWattsPerM2Sr;
chromaticityXY = T_xyz(1:2,:)*radianceWattsPerM2Sr/sum(T_xyz*radianceWattsPerM2Sr);

%% Load cone spectral sensitivities
load T_cones_ss2
T_cones = SplineCmf(S_cones_ss2,T_cones_ss2,S);

%% Load in a directly measured sunlight through window
% and off piece of white paper towel on floor for comparison.
% That seems like something safe to look at and gives us a reality check.
load spd_phillybright
spd_phillybright = SplineSpd(S_phillybright,spd_phillybright,S);
photopicLuminancePhillyBrightCdM2 = T_xyz(2,:)*spd_phillybright;
OLSLratio = radianceWattsPerM2Sr./spd_phillybright;

%% Compute irradiance, trolands, etc.
pupilDiamMm = 4.7;
pupilDiamMm = GetWithDefault('Enter observer pupil diameter in mm',pupilDiamMm);
pupilAreaMm2 = pi*((pupilDiamMm/2)^2);
eyeLengthMm = 17;
degPerMm = RetinalMMToDegrees(1,eyeLengthMm);
irradianceWattsPerUm2 = RadianceToRetIrradiance(radianceWattsPerM2Sr,S,pupilAreaMm2,eyeLengthMm);
irradianceScotTrolands = RetIrradianceToTrolands(irradianceWattsPerUm2, S, 'Scotopic', [], num2str(eyeLengthMm));
irradiancePhotTrolands = RetIrradianceToTrolands(irradianceWattsPerUm2, S, 'Photopic', [], num2str(eyeLengthMm));
irradianceQuantaPerUm2Sec = EnergyToQuanta(S,irradianceWattsPerUm2);
irradianceWattsPerCm2 = (10.^8)*irradianceWattsPerUm2;
irradianceQuantaPerCm2Sec = (10.^8)*irradianceQuantaPerUm2Sec;
irradianceQuantaPerDeg2Sec = (degPerMm^2)*(10.^-2)*irradianceQuantaPerCm2Sec;

%% Pupil adjustment factor for Ansi MPE 
mpePupilDiamMm = 3;
mpePupilDiamMm  = GetWithDefault('Enter ANSI 2007 MPE caclulations assumed pupil diameter in mm',mpePupilDiamMm );
pupilAdjustFactor = (pupilDiamMm/mpePupilDiamMm).^2;

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
irradianceScotTrolands_check = pupilAreaMm2*1700*(T_scotopicVlambda*radianceWattsPerM2Sr);
irradiancePhotTrolands_check = pupilAreaMm2*photopicLuminanceCdM2;

%% Get cone coordinates from radiance, and also adjust by pupil area.
% Useful for comparing to light levels produced by monochromatic lights
% in other papers
theLMS = T_cones*radianceWattsPerM2Sr;
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
cornealIrradianceWattsPerM2 = RadianceAndDistanceAreaToCornIrradiance(radianceWattsPerM2Sr,stimulusDistanceM,stimulusAreaM2);
cornealIrradianceWattsPerCm2 = (10.^-4)*cornealIrradianceWattsPerM2;
cornealIrradianceQuantaPerCm2Sec = EnergyToQuanta(S,cornealIrradianceWattsPerCm2);

%% Report on stimulus
fprintf('\n');
fprintf('  * Stimulus diameter mm %0.1f, degrees %0.1f\n',2*stimulusRadiusMm,2*stimulusRadiusDeg);
fprintf('  * Stimulus radiance %0.1f log10 watts/[m2-sr], %0.1f log10 watts/[cm2-sr]\n',log10(sum(radianceWattsPerM2Sr)),log10(sum(radianceWattsPerCm2Sr)));
fprintf('  * Stimulus luminance %0.1f candelas/m2\n',photopicLuminanceCdM2);
fprintf('  * Stimulus chromaticity x=%0.4f, y=%0.4f\n',chromaticityXY(1), chromaticityXY(2));
fprintf('    * For comparison, sunlight in Philly: %0.1f cd/m2\n',photopicLuminancePhillyBrightCdM2);
fprintf('  * Stimulus %0.0f (check val %0.0f) scotopic trolands, %0.0f photopic trolands (check val %0.0f)\n',irradianceScotTrolands,irradianceScotTrolands_check,...
    irradiancePhotTrolands,irradiancePhotTrolands_check);
fprintf('  * Stimulus %0.1f log10 scotopic trolands, %0.1f log10 photopic trolands\n',log10(irradianceScotTrolands),log10(irradiancePhotTrolands));
fprintf('  * Stimulus retinal irradiance %0.1f log10 watts/cm2\n',log10(sum(irradianceWattsPerCm2)));
fprintf('  * Stimulus retinal irradiance %0.1f log10 quanta/[cm2-sec]\n',log10(sum(irradianceQuantaPerCm2Sec)));
fprintf('  * Stimulus retinal irradiance %0.1f log10 quanta/[deg2-sec]\n',log10(sum(irradianceQuantaPerDeg2Sec)));
fprintf('  * Stimulus corneal irradiance %0.1f log10 watts/cm2\n',log10(sum(cornealIrradianceWattsPerCm2)));
fprintf('  * Stimulus corneal irradiance %0.1f log10 quanta/[cm2-sec]\n',log10(sum(cornealIrradianceQuantaPerCm2Sec)));
fprintf('  * Pupil area times LMS: %0.2f, %0.2f, %0.2f\n',...
        theLMSTimesPupilArea(1),theLMSTimesPupilArea(2),theLMSTimesPupilArea(3));

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
clear MPELimitIntegratedRadiance_JoulesPerCm2Sr MPELimitRadiance_WattsPerCm2Sr MPELimitCornealIrradiance_WattsPerCm2 MPELimitCornealRadiantExposure_JoulesPerCm2
for w = 1:length(stimulusWavelengthsNm)
    stimulusWavelengthNm = stimulusWavelengthsNm(w);
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
fprintf('    * Compare with total stimulus radiance %0.1f log  W/[cm2-sr]\n',log10(sum(radianceWattsPerCm2Sr)));
fprintf('    * Compare with total pupil adjusted radiance %0.1f log  W/[cm2-sr]\n',log10(sum(radianceWattsPerCm2Sr))+log10(pupilAdjustFactor));
fprintf('    * Pupil adjustment assumes observer pupil diameter of %0.1f mm, MPE standard diameter of %0.1f mm\n',pupilDiamMm,mpePupilDiamMm);

%% Now compare to the ISO Standard
stimulusDurationForISOMPESecs = 120*60;
[IsOverLimit,ISO2007MPEStruct] = ISO2007MPECheckType1ContinuousRadiance(S,pupilAdjustFactor*radianceWattsPerM2Sr,stimulusDurationForISOMPESecs,stimulusAreaDegrees2,eyeLengthMm);
fprintf('\n');
fprintf('  * ISO MPE Analysis\n');
ISO2007MPEPrintAnalysis(IsOverLimit,ISO2007MPEStruct);
fprintf('  * Assumed duration seconds %0.1f, hours %0.1f\n',stimulusDurationForISOMPESecs,stimulusDurationForISOMPESecs/3600);

%% Root name for plots
plotRoot = sprintf('MPEPlot_%d_%d',10*pupilDiamMm,10*mpePupilDiamMm);

%% Plot of stimulus radiance
%    Black solid, our spectrum
%    Black dashed, our spectrum bumped up by pupilAdjustFactor
%    Blue, sunlight measured off a piece of paper in my office in Philly
%    Red, ANSI MPE as a function of wavelength (power per band)
fig2 = figure; clf; hold on
set(gcf,'Position',[127         198         900 700]);
set(gca,'FontName','Helvetica','FontSIze',18);
log10radianceWattsPerCm2Sr = log10(radianceWattsPerCm2Sr);
log10radianceWattsPerCm2Sr(log10radianceWattsPerCm2Sr < -15) = NaN;
plot(SToWls(S),log10radianceWattsPerCm2Sr,'k','LineWidth',2);
plot(SToWls(S),log10radianceWattsPerCm2Sr+log10(pupilAdjustFactor),'k:','LineWidth',2);
plot(SToWls(S),log10(1e-4*spd_phillybright),'b:','LineWidth',2);
plot(stimulusWavelengthsNm,log10(MPELimitRadiance_WattsPerCm2Sr),'r','LineWidth',3);
xlabel('Wavelength');
ylabel('Radiance (W/[cm2-sr-wlband]');
theTitle{1} = sprintf('Luminance %0.1f cd/m2, total radiance %0.1f log10 watts/cm2-sr',photopicLuminanceCdM2,log10(sum(radianceWattsPerCm2Sr)));
theTitle{2} = sprintf('Pupil %0.1f mm, MPE Assumed Pupil %0.1f mm',pupilDiamMm,mpePupilDiamMm);
theTitle{3} = 'Black - spectrum, Black dashed - pupil adjusted, Blue - Philly sunlight, Red - Ansi MPE';
title(theTitle,'FontSize',16);


