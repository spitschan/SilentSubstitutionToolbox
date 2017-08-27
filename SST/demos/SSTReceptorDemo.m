%% SSTReceptorDemo.m
%
% This script demonstrates the use of the SSTReceptorHuman object and its
% ability to investigate the splatter either (A) parametrically or using (B) a
% resampling approach.
%
% (A) The SSTReceptorHuman object can create variants of the spectral
%     sensitivities which differ in the 8-parameter Asano et al. model for
%     cone individual differences. The parameters are lens density, macular
%     pigment  density, the optical density for the LMS pigments, and the
%     shift in lambda-max of the LMS pigments. Once the SSTReceptor object
%     has been  defined, this is done through a call to
%     makeSpectralSensitivitiesParametricVariation, which is a method of
%     the object.
%
% (B) The SSTReceptorHuman object can also create variants of the spectral
%     sensitivities by resampling them using the known (from the
%     literature) standard deviations of the individual differences
%     parameters. In any given resampled set of LMS cone fundamentals, the
%     parameters have been changed. Of course, the resampling is set up
%     such that the lens and macular pigment density parameters affect the
%     LMS cone fundamentals, while the optical density and the lambda-max
%     shift parameters affect the LMS cones independently.
%
% Under the hood, the SSTReceptorHuman object uses machinery from
% Psychtoolbox-3 to generate the spectral sensitivities of the cones.
%
% The demo loads background and modulation spectra which are part of the
% SilentSubstitutionToolbox. These were added when we set up SST originally
% and are not of specific interest beyond showing how to load spectra and
% make splatter calculations.
%
% 7/22/17   ms      Commented.

%% Clear work space
clearvars; close all; clc;

% Check if the SST is installed properly
CheckSSTInstalled();

% Infer the SST root
sstRoot0 = mfilename('fullpath');
sstRoot1 = cd(fullfile(fileparts(sstRoot0), '../..'));
sstRoot = pwd;

%% Get plot colors
theRGB = DefaultReceptorColors;

%% Load backgound and modulation spectra from some demo data in SST
tmp = load(fullfile(sstRoot, 'ContrastSplatter/ContrastSplatterDemoData/spd_contrastsplatterdemo_bg.mat'));
bgSpd = tmp.spd;
tmp = load(fullfile(sstRoot, 'ContrastSplatter/ContrastSplatterDemoData/spd_contrastsplatterdemo_mod.mat'));
modSpd = tmp.spd;

%% Set up receptor object
% This creates a static version of the spectral sensitivities which
% correspond to the parameters passed into the SSTReceptorHuman function.
receptorObj = SSTReceptorHuman('verbosity', 'high', 'obsAgeInYrs', 32, 'doPenumbralConesTrueFalse', true);

%% Plot the fundamentals in various formats
%receptorObj.plotSpectralSensitivities('whichFormat', 'T_quantalIsomerizations', 'saveFigPath', fullfile(sstRoot, 'SST', 'plots'));
%receptorObj.plotSpectralSensitivities('whichFormat', 'T_quantalAbsorptions', 'saveFigPath', fullfile(sstRoot, 'SST', 'plots'));
%receptorObj.plotSpectralSensitivities('whichFormat', 'T_quantalAbsorptionsNormalized', 'saveFigPath', fullfile(sstRoot, 'SST', 'plots'));
%receptorObj.plotSpectralSensitivities('whichFormat', 'T_energy', 'saveFigPath', fullfile(sstRoot, 'SST', 'plots'));
receptorObj.plotSpectralSensitivities('whichFormat', 'T_energyNormalized', 'saveFigPath', fullfile(sstRoot, 'SST', 'plots'));

%% Parametric variation of individual difference parameters
% Define the individual difference parameters that we want to look at here.
% The possible ones are:
%   dlens - Lens density
%   dmac - Macular pigment density
%   dphotopigment{L|M|S} - Optical density of the LMS photopgiments,
%                          respectively
%   lambdaMaxShift{L|M|S} - Shift of lambda-max for the LMS photopigments
%   obsPupilDiameterMm - Pupil diameter of the observer
theIndDiffParams = {'dlens' 'dmac' 'dphotopigmentL' 'dphotopigmentM' 'dphotopigmentS' ...
    'lambdaMaxShiftL' 'lambdaMaxShiftM' 'lambdaMaxShiftS' 'obsPupilDiameterMm'};

% Set the number sampling spacing for each of the individual difference
% parameters that we're interested in.
NTitrations = 16;

% Iterate over the parameters. This creates versions of the spectral
% sensitivities which differ in the value of the individual difference
% parameters. This gets saved out in the "Tp" field of the receptor object.
for ss = 1:length(theIndDiffParams)
    [~, parv{ss}, parvlabel{ss}, parvlabellong{ss}, parvreal{ss}] = makeSpectralSensitivitiesParametricVariation(receptorObj, ...
        'WhichParameter', theIndDiffParams{ss}, 'NTitrations', NTitrations);
end

%% Stochastic variation of individual difference parameters by resampling
% Instead of varying only one parameter, as we did above, we can also
% resample the spectral sensitivities, drawing from the distribution of the
% individual difference parameters.

% Define the number of samples.
NSamples = 200;

% Resample! This gets saved out in the "Ts" field of the receptor object.
receptorObj.makeSpectralSensitivitiesStochastic('NSamples', NSamples);

% Save the receptor object hash
receptorObj.setMD5Hash();

%% Calculate contrast and plot - Parametric
% Plot the contrast on the LMS cones as a function of the individual
% parameters which we changed above  under "Parametric variation of
% individual difference parameters".
figParv = figure;
yAxLims = [-0.06 0.06];
for ss = 1:length(theIndDiffParams)
    if parv{ss}(1) < 0
        xAxLims = [parv{ss}(1)*1.1 parv{ss}(end)*1.1];
    else
        xAxLims = [parv{ss}(1)*0.9 parv{ss}(end)*1.1];
    end
    
    % Calculate contrast on the LMS cones and the postreceptoral
    % combinations
    for ii = 1:size(receptorObj.Tp, 2)
        T_receptors = receptorObj.Tp{ss, ii}.T_energyNormalized;
        for jj = 1:size(T_receptors, 1)
            contrasts(jj, ii) = (T_receptors(jj, :)*(modSpd-bgSpd))./(T_receptors(jj, :)*bgSpd);
        end
        lmContrast(:, ii) = [1 1 0  0 0]' \ contrasts(:, ii);
        postRecepContrasts(:, ii) = [1 1 1 0 0 0 0 0 ; 1 -1 0 0 0 0 0 0 ; 0 0 1 0 0 0 0 0]' \ contrasts(:, ii);
    end
    
    % Plot contrast as a function of the individual difference parameters
    subplot(1, length(theIndDiffParams), ss);
    hold on;
    plot(xAxLims, [0 0], ':k');
    for ii = 1:size(contrasts, 1)
        plot(parv{ss}, contrasts(ii, :)', '-o', 'Color', theRGB(ii, :), 'MarkerEdgeColor', 'k', ...
            'MarkerFaceColor', theRGB(ii, :)); hold on;
    end
    
    % Add title and tweak plots
    title(parvlabellong{ss});
    box off;
    xlabel(parvlabel{ss});
    xlim(xAxLims);
    ylim(yAxLims);
    pbaspect([1 1 1]);
    set(gca, 'TickDir', 'out');
end

% Save out the figure
set(figParv, 'PaperPosition', [0 0 13 3.5]);
set(figParv, 'PaperSize', [13 3.5]);
set(figParv, 'Color', 'w');
set(figParv, 'InvertHardcopy', 'off');
saveas(figParv, fullfile(sstRoot, 'SST', 'plots', 'SSTReceptorDemo_ParvFig.png'), 'png');

%% Calculate contrast and plot - Resampling approach
% Calculate contrast for each of the spectral sensitivities in the Ts
% field.
clear contrasts postRecepContrasts;
for ii = 1:NSamples
    T_receptors = receptorObj.Ts{ii}.T_energyNormalized;
    for jj = 1:size(receptorObj.Ts{ii}.T_energyNormalized, 1)
        contrasts(jj, ii) = (T_receptors(jj, :)*(modSpd-bgSpd))./(T_receptors(jj, :)*bgSpd);
    end
    postRecepContrasts(:, ii) = [1 1 1 0 0 0 0 0; 1 -1 0 0 0 0 0 0; 0 0 1 0 0 0 0 0]' \ contrasts(:, ii);
end

%% [1] L+M+S vs. L-M contrast
figLMSvsLMinusM = figure;
XNominalContrast = 0;
YNominalContrast = 0;
XAxLims = XNominalContrast+[-0.03 0.03];
YAxLims = YNominalContrast+[-0.03 0.03];
XBinWidth = max(XAxLims)/10;
YBinWidth = max(YAxLims)/10;

% Throw it in a scatter plot
ScatterplotWithHistogram(postRecepContrasts(1, :), postRecepContrasts(2, :), ...
    'XLim', XAxLims, 'YLim', YAxLims, 'XBinWidth', XBinWidth, 'YBinWidth', YBinWidth, ...
    'XLabel', 'L+M+S contrast', 'YLabel', 'L-M contrast', ...
    'XRefLines', [XAxLims ; YNominalContrast YNominalContrast], ...
    'YRefLines', [XNominalContrast XNominalContrast ; YAxLims], ...
    'XNominalContrast', XNominalContrast, ...
    'YNominalContrast', YNominalContrast, ...
    'Color', [theRGB(1, :) ; theRGB(2, :)]);

% Plot mean
plot(mean(postRecepContrasts(1, :)), mean(postRecepContrasts(2, :)), '+k');

% Get the error ellipse
[X, Y] = get_error_ellipse([postRecepContrasts(1, :) ; postRecepContrasts(2, :)]', 0.95);
plot(X, Y, '-k');

% Save out the figures
set(figLMSvsLMinusM, 'PaperPosition', [0 0 6 6]);
set(figLMSvsLMinusM, 'PaperSize', [6 6]);
set(figLMSvsLMinusM, 'Color', 'w');
set(figLMSvsLMinusM, 'InvertHardcopy', 'off');
saveas(figLMSvsLMinusM, fullfile(sstRoot, 'SST', 'plots', 'SSTReceptorDemo_LMSvsLMinusM_Contrast.png'), 'png');

%% [2] L+M+S vs. S contrast
figLMSvsS = figure;
XNominalContrast = 0;
YNominalContrast = 0;
XAxLims = XNominalContrast+[-0.03 0.03];
YAxLims = YNominalContrast+[-0.1 0.1];
XBinWidth = max(XAxLims)/10;
YBinWidth = max(YAxLims)/10;

% Throw it in a scatter plot
ScatterplotWithHistogram(postRecepContrasts(1, :), postRecepContrasts(3, :), ...
    'XLim', XAxLims, 'YLim', YAxLims, 'XBinWidth', XBinWidth, 'YBinWidth', YBinWidth, ...
    'XLabel', 'L+M+S contrast', 'YLabel', 'S contrast', ...
    'XRefLines', [XAxLims ; YNominalContrast YNominalContrast], ...
    'YRefLines', [XNominalContrast XNominalContrast ; YAxLims], ...
    'XNominalContrast', XNominalContrast, ...
    'YNominalContrast', YNominalContrast, ...
    'Color', [theRGB(1, :) ; theRGB(3, :)]);

% Plot mean
plot(mean(postRecepContrasts(1, :)), mean(postRecepContrasts(3, :)), '+k');

% Get the error ellipse
[X, Y] = get_error_ellipse([postRecepContrasts(1, :) ; postRecepContrasts(3, :)]', 0.95);
plot(X, Y, '-k');

% Save out the figures
set(figLMSvsS, 'PaperPosition', [0 0 6 6]);
set(figLMSvsS, 'PaperSize', [6 6]);
set(figLMSvsS, 'Color', 'w');
set(figLMSvsS, 'InvertHardcopy', 'off');
saveas(figLMSvsS, fullfile(sstRoot, 'SST', 'plots', 'SSTReceptorDemo_LMSvsS_Contrast.png'), 'png');

%%
wls = SToWls(receptorObj.S);
wlIdxLens = find(wls == 440);
wlIdxMac = find(wls == 450);
for ii = 1:NSamples
    lensTransmittance(ii) = receptorObj.Ts{ii}.indDiffParams.dlens;
    macTransmittance(ii) = receptorObj.Ts{ii}.indDiffParams.dmac;
    dphotopigmentL(ii) = receptorObj.Ts{ii}.indDiffParams.dphotopigment(1);
    dphotopigmentM(ii) = receptorObj.Ts{ii}.indDiffParams.dphotopigment(2);
    dphotopigmentS(ii) = receptorObj.Ts{ii}.indDiffParams.dphotopigment(3);
    lambdaMaxShiftL(ii) = receptorObj.Ts{ii}.indDiffParams.lambdaMaxShift(1);
    lambdaMaxShiftM(ii) = receptorObj.Ts{ii}.indDiffParams.lambdaMaxShift(2);
    lambdaMaxShiftS(ii) = receptorObj.Ts{ii}.indDiffParams.lambdaMaxShift(3);
end

%% Plot the contrasts as a function of the parameter values from the resampling
whichContrastToPlot = contrasts;
%whichContrastToPlot = postRecepContrasts;

%% Lens
subplot(3, 8, 1);
plot(lensTransmittance, whichContrastToPlot(1, :), '.', 'Color', theRGB(1, :));
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.03 0.03]);
ylabel('Contrast');

subplot(3, 8, 9);
plot(lensTransmittance, whichContrastToPlot(2, :), '.', 'Color', theRGB(2, :));
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.03 0.03]);
ylabel('Contrast');

subplot(3, 8, 17);
plot(lensTransmittance, whichContrastToPlot(3, :), '.', 'Color', theRGB(3, :));
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.1 0.1]);
ylabel('Contrast');

% Macula
subplot(3, 8, 2);
plot(macTransmittance, whichContrastToPlot(1, :), '.', 'Color', theRGB(1, :));
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.03 0.03]);

subplot(3, 8, 10);
plot(macTransmittance, whichContrastToPlot(2, :), '.', 'Color', theRGB(2, :));
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.03 0.03]);

subplot(3, 8, 18);
plot(macTransmittance, whichContrastToPlot(3, :), '.', 'Color', theRGB(3, :));
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.1 0.1]);

% dphotopigmentL
subplot(3, 8, 3);
plot(dphotopigmentL, whichContrastToPlot(1, :), '.', 'Color', theRGB(1, :));
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.03 0.03]);

subplot(3, 8, 11);
plot(dphotopigmentL, whichContrastToPlot(2, :), '.', 'Color', theRGB(2, :));
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.03 0.03]);

subplot(3, 8, 19);
plot(dphotopigmentL, whichContrastToPlot(3, :), '.', 'Color', theRGB(3, :));
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.1 0.1]);

% dphotopigmentM
subplot(3, 8, 4);
plot(dphotopigmentM, whichContrastToPlot(1, :), '.', 'Color', theRGB(1, :));
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.03 0.03]);

subplot(3, 8, 12);
plot(dphotopigmentM, whichContrastToPlot(2, :), '.', 'Color', theRGB(2, :));
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.03 0.03]);

subplot(3, 8, 20);
plot(dphotopigmentM, whichContrastToPlot(3, :), '.', 'Color', theRGB(3, :));
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.1 0.1]);

% dphotopigmentS
subplot(3, 8, 5);
plot(dphotopigmentS, whichContrastToPlot(1, :), '.', 'Color', theRGB(1, :));
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.03 0.03]);

subplot(3, 8, 13);
plot(dphotopigmentS, whichContrastToPlot(2, :), '.', 'Color', theRGB(2, :));
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.03 0.03]);

subplot(3, 8, 21);
plot(dphotopigmentS, whichContrastToPlot(3, :), '.', 'Color', theRGB(3, :));
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.1 0.1]);

% lambdaMaxShiftL
subplot(3, 8, 6);
plot(lambdaMaxShiftL, whichContrastToPlot(1, :), '.', 'Color', theRGB(1, :));
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.03 0.03]);

subplot(3, 8, 14);
plot(lambdaMaxShiftL, whichContrastToPlot(2, :), '.', 'Color', theRGB(2, :));
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.03 0.03]);

subplot(3, 8, 22);
plot(lambdaMaxShiftL, whichContrastToPlot(3, :), '.', 'Color', theRGB(3, :));
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.1 0.1]);

% lambdaMaxShiftM
subplot(3, 8, 7);
plot(lambdaMaxShiftM, whichContrastToPlot(1, :), '.', 'Color', theRGB(1, :));
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.03 0.03]);

subplot(3, 8, 15);
plot(lambdaMaxShiftM, whichContrastToPlot(2, :), '.', 'Color', theRGB(2, :));
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.03 0.03]);

subplot(3, 8, 23);
plot(lambdaMaxShiftM, whichContrastToPlot(3, :), '.', 'Color', theRGB(3, :));
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.1 0.1]);

% lambdaMaxShiftS
subplot(3, 8, 8);
plot(lambdaMaxShiftS, whichContrastToPlot(1, :), '.', 'Color', theRGB(1, :));
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.03 0.03]);

subplot(3, 8, 16);
plot(lambdaMaxShiftS, whichContrastToPlot(2, :), '.', 'Color', theRGB(2, :));
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.03 0.03]);

subplot(3, 8, 24);
plot(lambdaMaxShiftS, whichContrastToPlot(3, :), '.', 'Color', theRGB(3, :));
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; ylim([-0.1 0.1]);

%%
%F = [lensTransmittance ; macTransmittance ; dphotopigmentL ; dphotopigmentM ; dphotopigmentS ; lambdaMaxShiftL ; lambdaMaxShiftM ; lambdaMaxShiftS]