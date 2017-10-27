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

% Infer the SST root
sstRoot0 = mfilename('fullpath');
sstRoot1 = cd(fullfile(fileparts(sstRoot0), '../..'));
sstRoot = pwd;

% Check if the SST is installed properly
CheckSSTInstalled();

%% Load backgound and modulation spectra from some demo data in SST
%
% We magically know that the wavelength sampling of these spectra is
% 380nm to 780 nm in 2nm increments.
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

%% Get plot colors
theRGB = SSTDefaultReceptorColors(receptorObj.labels);

%% Print out nominal contrast
fprintf('* Contrast values:\n');
for ii = 1:size(receptorObj.T.T_energyNormalized, 1)
    theContrast = (receptorObj.T.T_energyNormalized(ii, :)*(modSpd-bgSpd)) ./ (receptorObj.T.T_energyNormalized(ii, :)*(bgSpd));
    if length(receptorObj.labels{ii}) < 6
        fprintf('  - <strong>%s</strong>: \t\t%+.2f%s\n', receptorObj.labels{ii}, 100*theContrast, '%');
    else
        fprintf('  - <strong>%s</strong>: \t%+.2f%s\n', receptorObj.labels{ii}, 100*theContrast, '%');
    end
end
fprintf('\n');

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

% Set the number of titrations for each of the individual difference
% parameter that we're interested in.
NTitrations = 16;

% Iterate over the parameters. This creates versions of the spectral
% sensitivities which differ in the value of the individual difference
% parameters. This gets saved out in the "Tp" field of the receptor object.
for ss = 1:length(theIndDiffParams)
    makeSpectralSensitivitiesParametricVariation(receptorObj, ...
        'whichParameter', theIndDiffParams{ss}, 'NTitrations', NTitrations);
end

%% Calculate contrast and plot - Parametric
% Plot the contrast on the LMS cones as a function of the individual
% parameters which we changed above  under "Parametric variation of
% individual difference parameters".
figParv = figure;
yAxLims = [-0.06 0.06];
for ss = 1:length(theIndDiffParams)
    if receptorObj.Tp_i{ss}.parameterVariation.value(1) < 0
        xAxLims = [receptorObj.Tp_i{ss}.parameterVariation.value(1)*1.1 receptorObj.Tp_i{ss}.parameterVariation.value(end)*1.1];
    else
        xAxLims = [receptorObj.Tp_i{ss}.parameterVariation.value(1)*0.9 receptorObj.Tp_i{ss}.parameterVariation.value(end)*1.1];
    end
    
    % Calculate contrast on the LMS cones and the postreceptoral
    % combinations
    for ii = 1:size(receptorObj.Tp, 2)
        T_receptors = receptorObj.Tp{ss, ii}.T_energyNormalized;
        for jj = 1:size(T_receptors, 1)
            contrastsParametricVariation{ss}(jj, ii) = (T_receptors(jj, :)*(modSpd-bgSpd))./(T_receptors(jj, :)*bgSpd);
        end
        postRecepContrastsParametricVariation{ss}(:, ii) = [1 1 1 0 0 0 0 0 ; 1 -1 0 0 0 0 0 0 ; 0 0 1 0 0 0 0 0]' \ contrastsParametricVariation{ss}(:, ii);
    end
    
    % Plot contrast as a function of the individual difference parameters
    subplot(1, length(theIndDiffParams), ss);
    hold on;
    plot(xAxLims, [0 0], ':k');
    for ii = 1:3 % LMS only
        plot(receptorObj.Tp_i{ss}.parameterVariation.value, contrastsParametricVariation{ss}(ii, :)', '-o', 'Color', theRGB(ii, :), 'MarkerEdgeColor', 'k', ...
            'MarkerFaceColor', theRGB(ii, :)); hold on;
    end
    
    % Add title and tweak plots
    title(receptorObj.Tp_i{ss}.parameterVariation.label);
    box off;
    xlabel(receptorObj.Tp_i{ss}.parameterVariation.labelShort);
    xlim(xAxLims);
    ylim(yAxLims);
    pbaspect([1 1 1]);
    set(gca, 'TickDir', 'out');
end

% Enlargen the figure
screenDims = get(0, 'Screensize');
set(figParv, 'Position', [1 screenDims(4) screenDims(3)*0.8 screenDims(4)/3]);

% Save out the figure
set(figParv, 'PaperPosition', [0 0 40 10]);
set(figParv, 'PaperSize', [40 10]);
set(figParv, 'Color', 'w');
set(figParv, 'InvertHardcopy', 'off');
saveas(figParv, fullfile(sstRoot, 'SST', 'plots', 'SSTReceptorDemo_ParvFig.png'), 'png');

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

%% Calculate contrast and plot - Resampling approach
% Calculate contrast for each of the spectral sensitivities in the Ts
% field.
clear contrasts postRecepContrasts;
for ii = 1:NSamples
    T_receptors = receptorObj.Ts{ii}.T_energyNormalized;
    for jj = 1:size(receptorObj.Ts{ii}.T_energyNormalized, 1)
        contrastsResampling(jj, ii) = (T_receptors(jj, :)*(modSpd-bgSpd))./(T_receptors(jj, :)*bgSpd);
    end
    postRecepContrastsResampling(:, ii) = [1 1 1 0 0 0 0 0; 1 -1 0 0 0 0 0 0; 0 0 1 0 0 0 0 0]' \ contrastsResampling(:, ii);
end

%% [1] L+M+S vs. L-M contrast
figLMSvsLMinusM = figure;
XNominalContrast = 0;
YNominalContrast = 0;
XAxLims = XNominalContrast+[-0.03 0.03];
yAxLims = YNominalContrast+[-0.03 0.03];
XBinWidth = max(XAxLims)/10;
YBinWidth = max(yAxLims)/10;

% Throw it in a scatter plot
ScatterplotWithHistogram(postRecepContrastsResampling(1, :), postRecepContrastsResampling(2, :), ...
    'XLim', XAxLims, 'YLim', yAxLims, 'XBinWidth', XBinWidth, 'YBinWidth', YBinWidth, ...
    'XLabel', 'L+M+S contrast', 'YLabel', 'L-M contrast', ...
    'XRefLines', [XAxLims ; YNominalContrast YNominalContrast], ...
    'YRefLines', [XNominalContrast XNominalContrast ; yAxLims], ...
    'XNominalContrast', XNominalContrast, ...
    'YNominalContrast', YNominalContrast, ...
    'Color', [theRGB(1, :) ; theRGB(2, :)]);

% Plot mean
plot(mean(postRecepContrastsResampling(1, :)), mean(postRecepContrastsResampling(2, :)), '+k');

% Get the error ellipse
[X, Y] = get_error_ellipse([postRecepContrastsResampling(1, :) ; postRecepContrastsResampling(2, :)]', 0.95);
plot(X, Y, '-k');

% Save out the figures
set(figLMSvsLMinusM, 'PaperPosition', [0 0 12 12]);
set(figLMSvsLMinusM, 'PaperSize', [12 12]);
set(figLMSvsLMinusM, 'Color', 'w');
set(figLMSvsLMinusM, 'InvertHardcopy', 'off');
saveas(figLMSvsLMinusM, fullfile(sstRoot, 'SST', 'plots', 'SSTReceptorDemo_LMSvsLMinusM_Contrast.png'), 'png');

%% [2] L+M+S vs. S contrast
figLMSvsS = figure;
XNominalContrast = 0;
YNominalContrast = 0;
XAxLims = XNominalContrast+[-0.03 0.03];
yAxLims = YNominalContrast+[-0.1 0.1];
XBinWidth = max(XAxLims)/10;
YBinWidth = max(yAxLims)/10;

% Throw it in a scatter plot
ScatterplotWithHistogram(postRecepContrastsResampling(1, :), postRecepContrastsResampling(3, :), ...
    'XLim', XAxLims, 'YLim', yAxLims, 'XBinWidth', XBinWidth, 'YBinWidth', YBinWidth, ...
    'XLabel', 'L+M+S contrast', 'YLabel', 'S contrast', ...
    'XRefLines', [XAxLims ; YNominalContrast YNominalContrast], ...
    'YRefLines', [XNominalContrast XNominalContrast ; yAxLims], ...
    'XNominalContrast', XNominalContrast, ...
    'YNominalContrast', YNominalContrast, ...
    'Color', [theRGB(1, :) ; theRGB(3, :)]);

% Plot mean
plot(mean(postRecepContrastsResampling(1, :)), mean(postRecepContrastsResampling(3, :)), '+k');

% Get the error ellipse
[X, Y] = get_error_ellipse([postRecepContrastsResampling(1, :) ; postRecepContrastsResampling(3, :)]', 0.95);
plot(X, Y, '-k');

% Save out the figures
set(figLMSvsS, 'PaperPosition', [0 0 12 12]);
set(figLMSvsS, 'PaperSize', [12 12]);
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
whichContrastToPlot = contrastsResampling;
%whichContrastToPlot = postRecepContrasts;

% Open a new figure
figPars = figure;

% Set some axis limits
yAxLims = [-0.06 0.06];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lens
xAxLims = [-50 50];
subplot(3, 8, 1);
plot([0 0], yAxLims, ':', xAxLims, [0 0], ':', 'Color', [0.1 0.1 0.1]); hold on;
plot(lensTransmittance, whichContrastToPlot(1, :), '.', 'Color', theRGB(1, :));

% Also add the values from the parametric variation approach
ss = 1; ii = 1;
plot(receptorObj.Tp_i{ss}.parameterVariation.value, contrastsParametricVariation{ss}(ii, :)', '-k', 'LineWidth', 2); hold on;

pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; xlim(xAxLims); ylim(yAxLims);
ylabel('Contrast');  title('Lens density');

subplot(3, 8, 9);
plot([0 0], yAxLims, ':', xAxLims, [0 0], ':', 'Color', [0.1 0.1 0.1]); hold on;
plot(lensTransmittance, whichContrastToPlot(2, :), '.', 'Color', theRGB(2, :));

% Also add the values from the parametric variation approach
ss = 1; ii = 2;
plot(receptorObj.Tp_i{ss}.parameterVariation.value, contrastsParametricVariation{ss}(ii, :)', '-k', 'LineWidth', 2); hold on;

pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; xlim(xAxLims);  ylim(yAxLims);
ylabel('Contrast');

subplot(3, 8, 17);
plot([0 0], yAxLims, ':', xAxLims, [0 0], ':', 'Color', [0.1 0.1 0.1]); hold on;
plot(lensTransmittance, whichContrastToPlot(3, :), '.', 'Color', theRGB(3, :));

% Also add the values from the parametric variation approach
ss = 1; ii = 3;
plot(receptorObj.Tp_i{ss}.parameterVariation.value, contrastsParametricVariation{ss}(ii, :)', '-k', 'LineWidth', 2); hold on;

pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; xlim(xAxLims);  ylim(yAxLims);
ylabel('Contrast'); xlabel('%\DeltaD_{lens}');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Macula
xAxLims = [-100 100];
subplot(3, 8, 2);
plot([0 0], yAxLims, ':', xAxLims, [0 0], ':', 'Color', [0.1 0.1 0.1]); hold on;
plot(macTransmittance, whichContrastToPlot(1, :), '.', 'Color', theRGB(1, :));

% Also add the values from the parametric variation approach
ss = 2; ii = 1;
plot(receptorObj.Tp_i{ss}.parameterVariation.value, contrastsParametricVariation{ss}(ii, :)', '-k', 'LineWidth', 2); hold on;

pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; xlim(xAxLims);  ylim(yAxLims);
title('Macular density');

subplot(3, 8, 10);
plot([0 0], yAxLims, ':', xAxLims, [0 0], ':', 'Color', [0.1 0.1 0.1]); hold on;
plot(macTransmittance, whichContrastToPlot(2, :), '.', 'Color', theRGB(2, :));

% Also add the values from the parametric variation approach
ss = 2; ii = 2;
plot(receptorObj.Tp_i{ss}.parameterVariation.value, contrastsParametricVariation{ss}(ii, :)', '-k', 'LineWidth', 2); hold on;

pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; xlim(xAxLims);  ylim(yAxLims);

subplot(3, 8, 18);
plot([0 0], yAxLims, ':', xAxLims, [0 0], ':', 'Color', [0.1 0.1 0.1]); hold on;
plot(macTransmittance, whichContrastToPlot(3, :), '.', 'Color', theRGB(3, :));

% Also add the values from the parametric variation approach
ss = 2; ii = 3;
plot(receptorObj.Tp_i{ss}.parameterVariation.value, contrastsParametricVariation{ss}(ii, :)', '-k', 'LineWidth', 2); hold on;

pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; xlim(xAxLims);  ylim(yAxLims);
xlabel('%\DeltaD_{mac}');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dphotopigmentL
xAxLims = [-50 50];
subplot(3, 8, 3);
plot([0 0], yAxLims, ':', xAxLims, [0 0], ':', 'Color', [0.1 0.1 0.1]); hold on;
plot(dphotopigmentL, whichContrastToPlot(1, :), '.', 'Color', theRGB(1, :));

% Also add the values from the parametric variation approach
ss = 3; ii = 1;
plot(receptorObj.Tp_i{ss}.parameterVariation.value, contrastsParametricVariation{ss}(ii, :)', '-k', 'LineWidth', 2); hold on;

pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; xlim(xAxLims); ylim(yAxLims);
title({'Photopigment' ; 'optical density [L]'});

subplot(3, 8, 11);
plot([0 0], yAxLims, ':', xAxLims, [0 0], ':', 'Color', [0.1 0.1 0.1]); hold on;
plot(dphotopigmentL, whichContrastToPlot(2, :), '.', 'Color', theRGB(2, :));

% Also add the values from the parametric variation approach
ss = 3; ii = 2;
plot(receptorObj.Tp_i{ss}.parameterVariation.value, contrastsParametricVariation{ss}(ii, :)', '-k', 'LineWidth', 2); hold on;

pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; xlim(xAxLims);  ylim(yAxLims);

subplot(3, 8, 19);
plot([0 0], yAxLims, ':', xAxLims, [0 0], ':', 'Color', [0.1 0.1 0.1]); hold on;
plot(dphotopigmentL, whichContrastToPlot(3, :), '.', 'Color', theRGB(3, :));

% Also add the values from the parametric variation approach
ss = 3; ii = 3;
plot(receptorObj.Tp_i{ss}.parameterVariation.value, contrastsParametricVariation{ss}(ii, :)', '-k', 'LineWidth', 2); hold on;

pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; xlim(xAxLims);  ylim(yAxLims);
xlabel('%\DeltaD_{pigment} [L]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dphotopigmentM
xAxLims = [-50 50];
subplot(3, 8, 4);
plot([0 0], yAxLims, ':', xAxLims, [0 0], ':', 'Color', [0.1 0.1 0.1]); hold on;
plot(dphotopigmentM, whichContrastToPlot(1, :), '.', 'Color', theRGB(1, :));

% Also add the values from the parametric variation approach
ss = 4; ii = 1;
plot(receptorObj.Tp_i{ss}.parameterVariation.value, contrastsParametricVariation{ss}(ii, :)', '-k', 'LineWidth', 2); hold on;

pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; xlim(xAxLims);  ylim(yAxLims);
title({'Photopigment' ; 'optical density [M]'});

subplot(3, 8, 12);
plot([0 0], yAxLims, ':', xAxLims, [0 0], ':', 'Color', [0.1 0.1 0.1]); hold on;
plot(dphotopigmentM, whichContrastToPlot(2, :), '.', 'Color', theRGB(2, :));

% Also add the values from the parametric variation approach
ss = 4; ii = 2;
plot(receptorObj.Tp_i{ss}.parameterVariation.value, contrastsParametricVariation{ss}(ii, :)', '-k', 'LineWidth', 2); hold on;

pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; xlim(xAxLims);  ylim(yAxLims);

subplot(3, 8, 20);
plot([0 0], yAxLims, ':', xAxLims, [0 0], ':', 'Color', [0.1 0.1 0.1]); hold on;
plot(dphotopigmentM, whichContrastToPlot(3, :), '.', 'Color', theRGB(3, :));

% Also add the values from the parametric variation approach
ss = 4; ii = 3;
plot(receptorObj.Tp_i{ss}.parameterVariation.value, contrastsParametricVariation{ss}(ii, :)', '-k', 'LineWidth', 2); hold on;

pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; xlim(xAxLims);  ylim(yAxLims);
xlabel('%\DeltaD_{pigment} [M]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dphotopigmentS
xAxLims = [-50 50];
subplot(3, 8, 5);
plot([0 0], yAxLims, ':', xAxLims, [0 0], ':', 'Color', [0.1 0.1 0.1]); hold on;
plot(dphotopigmentS, whichContrastToPlot(1, :), '.', 'Color', theRGB(1, :));

% Also add the values from the parametric variation approach
ss = 5; ii = 1;
plot(receptorObj.Tp_i{ss}.parameterVariation.value, contrastsParametricVariation{ss}(ii, :)', '-k', 'LineWidth', 2); hold on;

pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; xlim(xAxLims);  ylim(yAxLims);
title({'Photopigment' ; 'optical density [S]'});

subplot(3, 8, 13);
plot([0 0], yAxLims, ':', xAxLims, [0 0], ':', 'Color', [0.1 0.1 0.1]); hold on;
plot(dphotopigmentS, whichContrastToPlot(2, :), '.', 'Color', theRGB(2, :));

% Also add the values from the parametric variation approach
ss = 5; ii = 2;
plot(receptorObj.Tp_i{ss}.parameterVariation.value, contrastsParametricVariation{ss}(ii, :)', '-k', 'LineWidth', 2); hold on;

pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; xlim(xAxLims);  ylim(yAxLims);

subplot(3, 8, 21);
plot([0 0], yAxLims, ':', xAxLims, [0 0], ':', 'Color', [0.1 0.1 0.1]); hold on;
plot(dphotopigmentS, whichContrastToPlot(3, :), '.', 'Color', theRGB(3, :));

% Also add the values from the parametric variation approach
ss = 5; ii = 3;
plot(receptorObj.Tp_i{ss}.parameterVariation.value, contrastsParametricVariation{ss}(ii, :)', '-k', 'LineWidth', 2); hold on;

pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; xlim(xAxLims);  ylim(yAxLims);
xlabel('%\DeltaD_{pigment} [S]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lambdaMaxShiftL
xAxLims = [-5 5];
subplot(3, 8, 6);
plot([0 0], yAxLims, ':', xAxLims, [0 0], ':', 'Color', [0.1 0.1 0.1]); hold on;
plot(lambdaMaxShiftL, whichContrastToPlot(1, :), '.', 'Color', theRGB(1, :));

% Also add the values from the parametric variation approach
ss = 6; ii = 1;
plot(receptorObj.Tp_i{ss}.parameterVariation.value, contrastsParametricVariation{ss}(ii, :)', '-k', 'LineWidth', 2); hold on;

pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; xlim(xAxLims);  ylim(yAxLims);
title({'Peak sensitivity' '\lambda_{max} [L]'});

subplot(3, 8, 14);
plot([0 0], yAxLims, ':', xAxLims, [0 0], ':', 'Color', [0.1 0.1 0.1]); hold on;
plot(lambdaMaxShiftL, whichContrastToPlot(2, :), '.', 'Color', theRGB(2, :));

% Also add the values from the parametric variation approach
ss = 6; ii = 2;
plot(receptorObj.Tp_i{ss}.parameterVariation.value, contrastsParametricVariation{ss}(ii, :)', '-k', 'LineWidth', 2); hold on;

pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; xlim(xAxLims);  ylim(yAxLims);

subplot(3, 8, 22);
plot([0 0], yAxLims, ':', xAxLims, [0 0], ':', 'Color', [0.1 0.1 0.1]); hold on;
plot(lambdaMaxShiftL, whichContrastToPlot(3, :), '.', 'Color', theRGB(3, :));

% Also add the values from the parametric variation approach
ss = 6; ii = 3;
plot(receptorObj.Tp_i{ss}.parameterVariation.value, contrastsParametricVariation{ss}(ii, :)', '-k', 'LineWidth', 2); hold on;

pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; xlim(xAxLims);  ylim(yAxLims);
xlabel('\Delta\lambda_{max} [L]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lambdaMaxShiftM
xAxLims = [-5 5];
subplot(3, 8, 7);
plot([0 0], yAxLims, ':', xAxLims, [0 0], ':', 'Color', [0.1 0.1 0.1]); hold on;
plot(lambdaMaxShiftM, whichContrastToPlot(1, :), '.', 'Color', theRGB(1, :));

% Also add the values from the parametric variation approach
ss = 7; ii = 1;
plot(receptorObj.Tp_i{ss}.parameterVariation.value, contrastsParametricVariation{ss}(ii, :)', '-k', 'LineWidth', 2); hold on;

pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; xlim(xAxLims); ylim(yAxLims);
title({'Peak sensitivity' '\lambda_{max} [M]'});

subplot(3, 8, 15);
plot([0 0], yAxLims, ':', xAxLims, [0 0], ':', 'Color', [0.1 0.1 0.1]); hold on;
plot(lambdaMaxShiftM, whichContrastToPlot(2, :), '.', 'Color', theRGB(2, :));

% Also add the values from the parametric variation approach
ss = 7; ii = 2;
plot(receptorObj.Tp_i{ss}.parameterVariation.value, contrastsParametricVariation{ss}(ii, :)', '-k', 'LineWidth', 2); hold on;

pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; xlim(xAxLims);  ylim(yAxLims);

subplot(3, 8, 23);
plot([0 0], yAxLims, ':', xAxLims, [0 0], ':', 'Color', [0.1 0.1 0.1]); hold on;
plot(lambdaMaxShiftM, whichContrastToPlot(3, :), '.', 'Color', theRGB(3, :));

% Also add the values from the parametric variation approach
ss = 7; ii = 3;
plot(receptorObj.Tp_i{ss}.parameterVariation.value, contrastsParametricVariation{ss}(ii, :)', '-k', 'LineWidth', 2); hold on;

pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; xlim(xAxLims);  ylim(yAxLims);
xlabel('\Delta\lambda_{max} [M]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lambdaMaxShiftS
xAxLims = [-5 5];
subplot(3, 8, 8);
plot([0 0], yAxLims, ':', xAxLims, [0 0], ':', 'Color', [0.1 0.1 0.1]); hold on;
plot(lambdaMaxShiftS, whichContrastToPlot(1, :), '.', 'Color', theRGB(1, :));

% Also add the values from the parametric variation approach
ss = 8; ii = 1;
plot(receptorObj.Tp_i{ss}.parameterVariation.value, contrastsParametricVariation{ss}(ii, :)', '-k', 'LineWidth', 2); hold on;

pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; xlim(xAxLims);  ylim(yAxLims);
title({'Peak sensitivity' '\lambda_{max} [S]'});

subplot(3, 8, 16);
plot([0 0], yAxLims, ':', xAxLims, [0 0], ':', 'Color', [0.1 0.1 0.1]); hold on;
plot(lambdaMaxShiftS, whichContrastToPlot(2, :), '.', 'Color', theRGB(2, :));

% Also add the values from the parametric variation approach
ss = 8; ii = 2;
plot(receptorObj.Tp_i{ss}.parameterVariation.value, contrastsParametricVariation{ss}(ii, :)', '-k', 'LineWidth', 2); hold on;

pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; xlim(xAxLims);  ylim(yAxLims);

subplot(3, 8, 24);
plot([0 0], yAxLims, ':', xAxLims, [0 0], ':', 'Color', [0.1 0.1 0.1]); hold on;
plot(lambdaMaxShiftS, whichContrastToPlot(3, :), '.', 'Color', theRGB(3, :));

% Also add the values from the parametric variation approach
ss = 8; ii = 3;
plot(receptorObj.Tp_i{ss}.parameterVariation.value, contrastsParametricVariation{ss}(ii, :)', '-k', 'LineWidth', 2); hold on;

pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; xlim(xAxLims);  ylim(yAxLims);
xlabel('\Delta\lambda_{max} [S]');

% Enlargen the figure
screenDims = get(0, 'Screensize');
set(figPars, 'Position', [1 screenDims(4) screenDims(3)*0.9 screenDims(4)*0.9]);


% Save out the figure
set(figPars, 'PaperPosition', [0 0 40 20]);
set(figPars, 'PaperSize', [40 20]);
set(figPars, 'Color', 'w');
set(figPars, 'InvertHardcopy', 'off');
saveas(figPars, fullfile(sstRoot, 'SST', 'plots', 'SSTReceptorDemo_ParvResamplingFig.png'), 'png');