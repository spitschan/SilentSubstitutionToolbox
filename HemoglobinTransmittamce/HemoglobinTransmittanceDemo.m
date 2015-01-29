% HemoglobinTransmittanceDemo.m
%
% Demo of how we compute hemoglobin absorptivity for retinal blood vessels.
%
% The molar excitation coefficients were obtained in tabulated fashion from
%   http://www.npsg.uwaterloo.ca/data/hemoglobin.php.
% They list the original data source as Scott Prahl, Optical Absorption of Hemoglobin, Oregon
% Medical Laser Center (http://omlc.ogi.edu/spectra/hemoglobin/index.html).
%
% The following website provides information on how to calculate into
% usable units:
%   http://omlc.ogi.edu/spectra/hemoglobin/index.html
%
% The demo loads the tabulated molar excitation coefficients from
% oxygenated and deoxygenated blood and plot them.
%
% 3/12/13   ms      Wrote it.
% 6/4/14    dhb     Review and comment.
%                   Change blood vessel thickness range computed for
%                   Make Prahl plot show Prahl data (it was still van Kampen)
%                   Sorted out and documented how Manuel's original calculation corresponds
%                     to a unit-based calculation from Prahl data.
% 12/15/14 dhb      A few cleanups.
% 1/29/15  dhb, ms  Fix up and comment about issue of log base 10 versus
%                   natural log and what the factor of 2.303 is.

%% Clear and close
clear; close all;

%% Run in the directory that contains this file.
cd(fileparts(mfilename('fullpath')));

%% Figure parameters
figureHeight = 600;
figureWidth = 600;
longFigureWidth = 1200;

%% Define folder to save figures into:
figFolder = fullfile(fileparts(mfilename), 'HemoglobinTransmittanceDemoOutput');
if ~isdir(figFolder)
    mkdir(figFolder);
end
origcheckDataFolder = fullfile(fileparts(mfilename), 'HemoglobinTransmittanceDemoData');
if ~isdir(origcheckDataFolder)
    mkdir(origcheckDataFolder);
end

%% Oxyhemoglobin calculations (Packer & Williams)
[~, wls_vanKampen, absorptivityLitersPerMmolPerCm_vanKampen] = GetHemoglobin('vanKampen');

% Packer & Williams (2003) say that the concentration of hemoglobin in
% blood is about 10 mmol liter^-1. We have the absorptivity of hemoglobin
% in units of liters mmol^-1 cm^-1. At max (415 nm), this absorbitivity is
% 131. Multiplying by 10 gives absorptivity as 1310 per cm.  {The term absorptivity
% is being a bit overloaded here, but with units specified it should be clear
% enough.)
absorptivityPerCm_vanKampen = absorptivityLitersPerMmolPerCm_vanKampen*10;

% According to Packer & Williams, blood vessels of 2 um lead to
% transmissivity of 0.55 at 415 nm.  2 um is 0.0002 cm. So this gives an absorptance
% of 1310*(0.0002) = 0.2620.
absorptance_vanKampen = absorptivityPerCm_vanKampen*0.0002;

% To get transmissivity, we take 10^(-0.2620) and get 0.547. This
% corresponds to the number given by Packer & Williams (2003). Yay!
% Other points on our transmissivity plot match up pretty well with
% the numberical values given by Packer and Williams for the transmissivity
transmissivity_vanKampen = 10.^(-absorptance_vanKampen);

% Plot out
%
% Absorptivity.  This matches up in shape with Figure 2.8 of Packer and Williams,
% although the units are different (they say theirs are arbitrary).
%
% The one mysetery is where Packer and Williams got data that go down below 415 nm,
% since it isn't in the table in van Kampen and Zijistra.
theVanKampenFig = figure;
set(gcf,'Position',[100 100 longFigureWidth figureHeight]);
subplot(1, 3, 1);
set(gca,'FontName','Helvetica','FontSize',16);
plot(wls_vanKampen, absorptivityPerCm_vanKampen, '-k', 'LineWidth', 2);
xlabel('Wavelength [nm]'); ylabel('Absorptivity per cm'); xlim([380 700]); pbaspect([1 1 1]);
title({'HbO2' 'Absorptivity [cm]' 'van Kampen & Zijlstra (1983)'});

% Absorptance.  This is just a scaled verion of what's in the first subpanel.
subplot(1, 3, 2);
set(gca,'FontName','Helvetica','FontSize',16);
plot(wls_vanKampen, absorptance_vanKampen, '-k', 'LineWidth', 2);
xlabel('Wavelength [nm]'); ylabel('Absorptance'); xlim([380 700]); pbaspect([1 1 1]);
title({'HbO2' 'Absorptance for 2 um vessel' 'van Kampen & Zijlstra (1983)'});

subplot(1, 3, 3);
set(gca,'FontName','Helvetica','FontSize',16);
plot(wls_vanKampen, transmissivity_vanKampen, '-k', 'LineWidth', 2);
xlabel('Wavelength [nm]'); ylabel('Transmissivity'); xlim([380 700]); ylim([0 1]); pbaspect([1 1 1]);
title({'HbO2' 'Transmissivity for 2 um vessel' 'van Kampen & Zijlstra (1983)'});
saveas(theVanKampenFig, fullfile(figFolder, 'Oxyhemoglobin_vanKampen'), 'pdf');

%% Explore now the effect of diameter on transmissivity in a reasonable
% range for retinal blood vessels. Snodderly, Weinhaus & Choi (1992, p.
% 1188) note that 3 um is the smallest opening through which a human red
% blood cell can pass undamaged. In their sample of macaque retina, the
% median capillary diameter for one subject (F46) was 5 um.
%
% Packer & Williams (2003) tell us that 'light incident on the
% photoreceptors is filtered by a vasculature that, if it were uniform,
% would be equivalent to a layer of blood about 2 um in thickness at
% eccentricities exceeding 1 mm (cf. their Fig. 2.7, p. 49). Their range of
% mean blood thinkness is 0.0-3.0 um. We calculate the tranmissivity in
% this range.  Note, however, that this is different from the thickness
% of a blood vessel itself -- it is an estimate of how thick a uniformly
% spread layer of blood would be.
minThicknessUm = 1;
maxThicknessUm = 10;
UmToCm = 1e-4;
theThicknessCm = linspace(minThicknessUm*UmToCm, maxThicknessUm*UmToCm, 10);
for t = 1:length(theThicknessCm)
    absorptanceAcrossDiameter_vanKampen(:, t) = absorptivityPerCm_vanKampen*theThicknessCm(t);
    transmissivityAcrossDiameter_vanKampen(:, t) = 10.^(-absorptanceAcrossDiameter_vanKampen(:, t));
end

% Plot this out
theVariationAcrossDiameterFig = figure; hold on
set(gcf,'Position',[100 100 figureWidth figureHeight]);
set(gca,'FontName','Helvetica','FontSize',16);
plot(wls_vanKampen, transmissivityAcrossDiameter_vanKampen(:, 1), '-b', 'LineWidth', 2);
plot(wls_vanKampen, transmissivityAcrossDiameter_vanKampen(:, end), '-r', 'LineWidth', 2);
plot(wls_vanKampen, transmissivityAcrossDiameter_vanKampen, '-k');
legend([num2str(theThicknessCm(1)/UmToCm) ' um'], [num2str(theThicknessCm(end)/UmToCm) ' um'], 'Location', 'SouthEast'); legend boxoff;
xlabel('Wavelength [nm]'); ylabel('Transmissivity'); xlim([380 700]); ylim([0 1]); pbaspect([1 1 1]);
title({'HbO2' 'Transmissivity' 'van Kampen & Zijlstra (1983)'});
saveas(theVariationAcrossDiameterFig,fullfile(figFolder,'Oxyhemoglobin_Thickness_vanKampen'),'pdf');

%% Obtain the Prahl hemoglobin estimates.
%
% We call our underlying routine with oxy fraction 0 and 1 and 2 uM
% thickness for this example
vesselThicknessUm = 2;
vesselThicknessCm = vesselThicknessUm*1e-4;
[transmittance_deoxy_Prahl,absorptance_deoxy_Prahl,S_Prahl] = GetHemoglobinTransmittance([],0,vesselThicknessUm,'Prahl');
[transmittance_oxy_Prahl,absorptance_oxy_Prahl] = GetHemoglobinTransmittance([],1,vesselThicknessUm,'Prahl');
wls_Prahl = SToWls(S_Prahl);

% For looking at the underlying data, we convert from the returned
% absorptance back to absorptivity by dividing out the blook vessel
% thickness in cm.
absorptivityPerCm_oxy_Prahl = absorptance_oxy_Prahl/vesselThicknessCm;
absorptivityPerCm_deoxy_Prahl = absorptance_deoxy_Prahl/vesselThicknessCm;

% Again, for plotting it is useful to have the molar extinction coefficients (cm-1/M), whatever an 'M' is.
[~,~,oxyMolarExtinction_Prahl, deoxyMolarExtinction_Prahl] = GetHemoglobin('Prahl');

%
% % Prahl writes:
% %     To convert this data to absorption coefficient in (cm-1), multiply by the molar concentration and 2.303,
% %         ua = (2.303) e (x g/liter)/(64,500 g Hb/mole)
% %     where x is the number of grams per liter. A typical value of x for whole blood is x=150 g Hb/liter.
% absorptivityPerCm_oxy_Prahl = (oxyMolarExtinction_Prahl * 150)/64500;
% absorptivityPerCm_deoxy_Prahl = (deoxyMolarExtinction_Prahl * 150)/64500;
%
% % Let's calculate now the absorptance given a 2 um blood vessel.

%
% % Get transmissivity. See above in the Packer & Williams/van Kampen calculations for the logic.
% transmittance_oxy_Prahl = 10.^(-absorptance_oxy_Prahl);
%  = 10.^(-absorptance_deoxy_Prahl);

% Plot out the Prahl data.
thePrahlFig = figure;
set(gcf,'Position',[100 100 longFigureWidth figureHeight]);
subplot(1, 3, 1); hold on
set(gca,'FontName','Helvetica','FontSize',16);
plot(wls_Prahl, absorptivityPerCm_oxy_Prahl, '-k', 'LineWidth', 2);
plot(wls_Prahl, absorptivityPerCm_deoxy_Prahl, '-r', 'LineWidth', 2);
xlabel('Wavelength [nm]'); ylabel('Absorptivity per cm');
xlim([380 700]);
pbaspect([1 1 1]);
title({'HbO2' 'Absorptivity [cm]' 'Prahl, Oregon Medical Laser Center'});
legend({'Oxyhemoglobin', 'Deoxyhemoglobin'});

subplot(1, 3, 2); hold on
set(gca,'FontName','Helvetica','FontSize',16);
plot(wls_Prahl, absorptance_oxy_Prahl, '-k', 'LineWidth', 2);
plot(wls_Prahl, absorptance_deoxy_Prahl, '-r', 'LineWidth', 2);
xlabel('Wavelength [nm]'); ylabel('Absorptance'); xlim([380 700]); pbaspect([1 1 1]);
title({'HbO2' 'Absorptance for 2 um vessel' 'Prahl, Oregon Medical Laser Center'});
legend({'Oxyhemoglobin', 'Deoxyhemoglobin'});

subplot(1, 3, 3); hold on
set(gca,'FontName','Helvetica','FontSize',16);
plot(wls_Prahl, transmittance_oxy_Prahl, '-k', 'LineWidth', 2);
plot(wls_Prahl, transmittance_deoxy_Prahl, '-r', 'LineWidth', 2);
xlabel('Wavelength [nm]'); ylabel('Transmissivity'); xlim([380 700]); ylim([0 1]); pbaspect([1 1 1]);
title({'HbO2' 'Transmissivity for 2 um vessel' 'Prahl, Oregon Medical Laser Center'});
legend({'Oxyhemoglobin', 'Deoxyhemoglobin'},'Location','SouthEast');
saveas(thePrahlFig, fullfile(figFolder, 'Oxyhemoglobin_Prahl'),  'pdf');

%% Compare now the two calculations (Packer & Williams/van Kampen and Prahl)
theComparisonFig = figure; hold on;
set(gcf,'Position',[100 100 figureWidth figureHeight]);
set(gca,'FontName','Helvetica','FontSize',16);
plot(wls_Prahl, transmittance_oxy_Prahl, '-k', 'LineWidth', 2);
plot(wls_Prahl, transmittance_deoxy_Prahl, '-r', 'LineWidth', 2);
plot(wls_vanKampen, transmissivity_vanKampen, '-b', 'LineWidth', 2);
xlabel('Wavelength [nm]'); ylabel('Transmissivity'); xlim([380 700]); ylim([0 1]); pbaspect([1 1 1]);
title({'HbO2' 'Transmissivity' 'Comparison, 2 um vessel'});
legend('Prahl Oxy, Oregon Medical Laser Center', 'Prahl Deoxy, Oregon Medical Laser Center', 'van Kampen & Zijlstra (1983)', ...
    'Location', 'SouthEast');
saveas(theComparisonFig, fullfile(figFolder, 'Oxyhemoglobin_Comparison_Prahl_vanKampen'), 'pdf');

%% Horiguchi et al. calculation
%
% Our original estimates were based on the Horiguchi et al. calculations.
%
% They start with tabulated data about oxy and deoxy hemoglobin absorptivity.
% Manuel imported these into files T_oxyhemoglobin and T_deoxyhemoglobin in
% this directory.  At the time, we didn't really understand the units.  But
% these functions turn out to be the Prahl data, which we show by comparing
% the two sets explicitly in the plot below.
%
% We also think they used an extra factor of 2.303 in their caluclations,
% see more extensive comments about this factor in GetHemoglobinTransmittance.
checkInputOxy = load(fullfile(origcheckDataFolder,'T_oxyhemoglobin'));
checkInputDeoxy = load(fullfile(origcheckDataFolder,'T_deoxyhemoglobin'));
S_Horiguchi_oxy = checkInputOxy.S_oxyhemoglobin;
S_Horiguchi_deoxy = checkInputDeoxy.S_deoxyhemoglobin;
oxyMolarExtinction_Horiguchi = checkInputOxy.T_oxyhemoglobin;
deoxyMolarExtinction_Horiguchi = checkInputDeoxy.T_deoxyhemoglobin;
clear checkInputOxy checkInputDeoxy

theFigure = figure; hold on;
set(gcf,'Position',[100 100 figureWidth figureHeight]);
set(gca,'FontName','Helvetica','FontSize',16);plot(SToWls(S_Horiguchi_oxy), oxyMolarExtinction_Horiguchi, '-k', 'LineWidth', 4);
plot(SToWls(S_Horiguchi_deoxy), deoxyMolarExtinction_Horiguchi, '-r', 'LineWidth', 4);
plot(wls_Prahl, oxyMolarExtinction_Prahl, ':y', 'LineWidth', 2);
plot(wls_Prahl, deoxyMolarExtinction_Prahl, ':b', 'LineWidth', 2);
xlabel('Wavelength [nm]');
ylabel('Molar Extinction');
legend('Oxyhemoglobin, Horiguchi et al.', 'Deoxyhemoglobin, Horiguchi et al.', 'Oxyhemoglobin, Prahl','Deoxyhemoglobin, Prahl');
pbaspect([1 1 1]);
saveas(theFigure, fullfile(figFolder, 'Horiguchi_Prahl_MolarExctinctionCompare'), 'pdf');

%% Geoff says that oxy hemoglobin is about 95% of hemoglobin in arteries and about 75% of
% hemoglobin in veins at room air oxygenation, based on the known partial pressure in these
% vessels and the hemoglobin-oxygen dissociation function, which he or any good med student
% could explain to us if we ask.
%
% Taking the average gives 85%.
%
% Who knows what the right thickness to use is for purposes of correcting the CIE
% cone sensitivities.  These probably have some amount of blood absorption effectively
% built in, via the fitting to measured 10-deg color matching functions.  And, penumbral
% cones may only have their light passed through partial thickness.  On the third hand,
% it is the bigger ones that we see when we flicker the melanopsin. We'll plot out how
% much it changes with a range of values, and then we'll just have to pick a set.
theFigure = figure; hold on;
set(gcf,'Position',[100 100 figureWidth figureHeight]);
set(gca,'FontName','Helvetica','FontSize',16);
thicknessesUm = [3 5 7];
oxyFractions = [0.75 0.85 0.95];
S = [380 4 101];
theColors = ['r' 'k' 'g' 'b' 'y' 'r' 'k' 'g' 'b' 'y'];
nextColor = 1;
for t = 1:length(thicknessesUm);
    for o = 1:length(oxyFractions)
        transmittance{t,o} = GetHemoglobinTransmittance(S,oxyFractions(o),thicknessesUm(t),'Prahl');
        normTransmittance{t,o} = transmittance{t,o}/max(transmittance{t,o}(:));
        plot(SToWls(S), normTransmittance{t,o}, theColors(nextColor), 'LineWidth', 4);
    end
    nextColor = nextColor + 1;
end
xlabel('Wavelength [nm]');
ylabel('Transmissivity');
title('Effect of vessel thickness and oxy fraction on transmittance');
pbaspect([1 1 1]);
saveas(theFigure, fullfile(figFolder, 'EffectOfOxyFractionAndThickness'),  'pdf');


%% **************************************************************
%
% This code makes some checks about how we used to compute hemoglogin
% transmittance with how we are doing it now.  We are keeping it here
% in case we ever want to go back, but it is not of general interest.
%
% In the end what this code demonstrates is that everying checked out
% when we had the erroneous extra factor of 2.303 in our calculations,
% and that now that we've fixed that the old estimates were a bit off.
%
% If you go into GetHemoglobinTranmittance and put back the factor 2.303
% then the curves in the graph that results will in line up again. 
DO_HISTORICAL_CHECKS = false;
if (DO_HISTORICAL_CHECKS)
    %% This is how we used the data to get blood vessel transmisivity
    %
    % The Wandell Lab has a function called s_CreateSkinAbsorptionFunctions, which
    % contains the appropriate conversion coefficients to go from molar
    % extinction to what they call absorption coefficients.  This consists
    % of multiplying by 5.4e-3.  Note that this equals 2.303*150/64500 to
    % a few places.  The latter is what we derived above for the Prahl data.  We are calling
    % this absorptivityPerCm....
    molExt2absCoef_Horiguchi = 5.4e-3;
    if (abs(molExt2absCoef_Horiguchi - (2.303*150/64500)) > 1e-4)
        error('Oops, conversion factor mismatch');
    end
    absorptivityPerCm_oxy_Horiguchi = molExt2absCoef_Horiguchi * oxyMolarExtinction_Horiguchi;
    absorptivityPerCm_deoxy_Horiguchi = molExt2absCoef_Horiguchi * deoxyMolarExtinction_Horiguchi;
    
    % Normalize.  We are not sure why Manuel decided to do this.  It moves us out of well understood
    % units into unknown units.  On the other hand, he didn't multiply by a vessel thickness.
    %
    % But, it is equivalent to multipying by some factor, and that factor is the thickness of the blood
    % vessels we were assuming, in cm.  We can print out that factor.
    absorptivityPerCm_oxy_HoriguchiNormalized = absorptivityPerCm_oxy_Horiguchi/max(absorptivityPerCm_oxy_Horiguchi);
    absorptivityPerCm_deoxy_HoriguchiNormalized = absorptivityPerCm_deoxy_Horiguchi/max(absorptivityPerCm_deoxy_Horiguchi);
    theEffectiveCmOxy = 1/max(absorptivityPerCm_oxy_Horiguchi);
    theEffectiveCmDeoxy = 1/max(absorptivityPerCm_deoxy_Horiguchi);
    fprintf('In Manuel''s original (April 2014) calculations, effective vessel thickness assumed was %0.1f um (oxy), %0.1f um (deoxy)\n', ...
        theEffectiveCmOxy/UmToCm,theEffectiveCmDeoxy/UmToCm);
    
    %% Incorporate hemoglobin transmittance into spectral sensitivities
    %
    % Take the mean of oxy and deoxy.  This is equivalent to assuming blood consists of oxy/deoxy in the ratio given
    % by the two effective vessel thicknesses above.
    if (any(S_Horiguchi_oxy ~= S_Horiguchi_deoxy))
        error('We think these two wavelength samplings should match');
    end
    absorptivityPerCm_oxyAndDeoxy_HoriguchiNormalized = mean([absorptivityPerCm_oxy_HoriguchiNormalized absorptivityPerCm_deoxy_HoriguchiNormalized], 2);
    hemoglobinTransmittance_oxyAndDeoxy_HoriguchiNormalized = 10.^(-absorptivityPerCm_oxyAndDeoxy_HoriguchiNormalized);
    %plot(SToWls(S_oxyhemoglobin), hemoglobinTransmittance_oxyAndDeoxy_HoriguchiNormalized)
    fprintf('So after taking the mean of oxy/deoxy absorptivity,\nthe overall transmisivity should be equal to an oxy layer of thickness %0.1f um and a deoxy layer of thickness %0.1f um\n',...
        (theEffectiveCmOxy/2)/UmToCm,(theEffectiveCmDeoxy/2)/UmToCm);
    fprintf('Oxy fraction: %0.2f\n',theEffectiveCmOxy/(theEffectiveCmOxy+theEffectiveCmDeoxy));
    
    % Here is what we wrote out originally, with the den_Hemoglobin data file
    % saved over in the BrainardLabToolbox at the time.
    S_Hemoglobin_origcheck = S_Horiguchi_oxy;
    den_Hemoglobin_origcheck = absorptivityPerCm_oxyAndDeoxy_HoriguchiNormalized;
    trans_Hemoglobin_origcheck = 10.^(-den_Hemoglobin_origcheck);
    save(fullfile(origcheckDataFolder,'den_Hemoglobin_origcheck'),'S_Hemoglobin_origcheck','den_Hemoglobin_origcheck');
    CHECK_ORIG = true;
    if (CHECK_ORIG)
        if (~exist(fullfile(origcheckDataFolder,'den_Hemoglobin_origmatfile.mat'),'file'))
            error('No den_Hemoglobin_origmatfile file available to compare against');
        end
        if (exist('den_Hemoglobin','var') || exist('S_Hemoglobin','var'))
            error('There is a variable called den_Hemoglobin or S_Hemoglobin already.')
        end
        origVersion = load(fullfile(origcheckDataFolder,'den_Hemoglobin_origmatfile.mat'));
        if (any(origVersion.S_Hemoglobin ~= S_Hemoglobin_origcheck))
            error('Mismatch between our current calculations of orignal values and what is in the loaded .mat file');
        else
            fprintf('Wavelength sampling matches original den_Hemoglobin file\n');
        end
        if (any(origVersion.den_Hemoglobin ~= den_Hemoglobin_origcheck))
            error('Mismatch between our current calculations of orignal values and what is in the loaded .mat file');
        else
            fprintf('den_Hemoglobin_origcheck computed now sampling matches original den_Hemoglobin computation\n');
        end
    end
    
    %% Compute transmisivity from Prahl data with effective thicknesses and make sure it matches up to what we expect
    % based on the above analysis
    checkOrigAbsorptance_oxy = absorptivityPerCm_oxy_Prahl*(theEffectiveCmOxy/2);
    checkOrigAbsorptance_deoxy = absorptivityPerCm_deoxy_Prahl*(theEffectiveCmDeoxy/2);
    checkOrigAbsorptance_oxyAndDeoxy = checkOrigAbsorptance_oxy + checkOrigAbsorptance_deoxy;
    checkOrigTrans_oxyAndDeoxy = 10.^(-checkOrigAbsorptance_oxyAndDeoxy);
    
    theFigure = figure;
    set(gcf,'Position',[100 100 round(2*longFigureWidth/3) figureHeight]);
    subplot(1,2,1); hold on;
    set(gca,'FontName','Helvetica','FontSize',16);
    plot(wls_Prahl,checkOrigAbsorptance_oxyAndDeoxy,'-r','LineWidth',4);
    plot(SToWls(S_Hemoglobin_origcheck),den_Hemoglobin_origcheck,'k:','LineWidth',2);
    title('Hemoglobin Absoprtance');
    xlabel('Wavelength [nm]');
    ylabel('Absorptance');
    legend('Prahl-based calc', 'Manuel''s original','Location','NorthEast');
    pbaspect([1 1 1]);
    subplot(1,2,2); hold on;
    set(gca,'FontName','Helvetica','FontSize',16);
    plot(wls_Prahl,checkOrigTrans_oxyAndDeoxy,'-r','LineWidth',4);
    plot(SToWls(S_Hemoglobin_origcheck),trans_Hemoglobin_origcheck,'k:','LineWidth',2);
    title('Hemoglobin Transmittance');
    xlabel('Wavelength [nm]');
    ylabel('Transmittance');
    legend('Prahl-based calc', 'Manuel''s original','Location','SouthEast');
    pbaspect([1 1 1]);
    saveas(theFigure, fullfile(origcheckDataFolder, 'RederivationOfApril2014Absorptance'),  'pdf');
end





