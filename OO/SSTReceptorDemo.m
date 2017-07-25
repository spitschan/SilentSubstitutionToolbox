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
% 7/22/17   ms      Commented.

%% Clear work space
clearvars; close all; clc;

% Check if the SST is installed properly
CheckSSTInstalled();

% Infer the SST root
sstRoot0 = mfilename('fullpath');
sstRoot1 = cd(fullfile(fileparts(sstRoot), '..'));
sstRoot = pwd;

%% Get some colors
theRGB = DefaultReceptorColors;

%% Load backgound and modulation spectra from some demo data in SST
tmp = load(fullfile(sstRoot, 'ContrastSplatter/ContrastSplatterDemoData/spd_contrastsplatterdemo_bg.mat'));
bgSpd = tmp.spd;
tmp = load(fullfile(sstRoot, 'ContrastSplatter/ContrastSplatterDemoData/spd_contrastsplatterdemo_mod.mat'));
modSpd = tmp.spd;

%% Set up receptor object
% This creates a static version of the spectral sensitivities which
% correspond to the parameters passed into the SSTReceptorHuman function.
receptorObj = SSTReceptorHuman('verbosity', 'high', 'obsAgeYrs', 32);
receptorObj.plotSpectralSensitivities;

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
        lmContrast(:, ii) = [1 1 0]' \ contrasts(:, ii);
        postRecepContrasts(:, ii) = [1 1 1 ; 1 -1 0 ; 0 0 1]' \ contrasts(:, ii);
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
saveas(figParv, fullfile(sstRoot, 'OO', 'plots', 'SSTReceptorDemo_ParvFig.png'), 'png');

%% Calculate contrast and plot - Resampling approach
% Calculate contrast for each of the spectral sensitivities in the Ts
% field.
for ii = 1:NSamples
    T_receptors = receptorObj.Ts{ii}.T_energyNormalized;
    for jj = 1:size(receptorObj.Ts{ii}.T_energyNormalized, 1)
        contrasts(jj, ii) = (T_receptors(jj, :)*(modSpd-bgSpd))./(T_receptors(jj, :)*bgSpd);
    end
    lmContrast(:, ii) = [1 1 0]' \ contrasts(:, ii);
    postRecepContrasts(:, ii) = [1 1 1 ; 1 -1 0 ; 0 0 1]' \ contrasts(:, ii);
end

% First, plot L and M contrast
axLims = [-0.03 0.03];
figLM = figure;
XNominalContrast = 0;
YNominalContrast = 0;
XAxLims = XNominalContrast+[-0.02 0.02];
YAxLims = YNominalContrast+[-0.02 0.02];

ScatterplotWithHistogram(contrasts(1, :), contrasts(2, :), ...
    'XLim', XAxLims, 'YLim', YAxLims, 'BinWidth', 0.004, ...
    'XLabel', 'L contrast', 'YLabel', 'M contrast', ...
    'XRefLines', [XAxLims ; YNominalContrast YNominalContrast], ...
    'YRefLines', [XNominalContrast XNominalContrast ; YAxLims], ...
    'XNominalContrast', XNominalContrast, ...
    'YNominalContrast', YNominalContrast, ...
    'Color', [theRGB(1, :) ; theRGB(2, :)]);

% Save out the figures
set(figLM, 'PaperPosition', [0 0 6 6]);
set(figLM, 'PaperSize', [6 6]);
set(figLM, 'Color', 'w');
set(figLM, 'InvertHardcopy', 'off');
saveas(figLM, fullfile(sstRoot, 'OO', 'plots', 'SSTReceptorDemo_LMContrast.png'), 'png');

% Then, plot L+M vs. S contrast
figS = figure;
XAxLims = [-0.02 0.02];
YAxLims = [-0.01 0.08];

XNominalContrast = 0;
YNominalContrast = 0;

ScatterplotWithHistogram(contrasts(1, :), contrasts(3, :), ...
    'XLim', XAxLims, 'YLim', YAxLims, 'BinWidth', 0.004, ...
    'XLabel', 'L+M contrast', 'YLabel', 'S contrast', ...
    'XRefLines', [XAxLims ; YNominalContrast YNominalContrast], ...
    'YRefLines', [XNominalContrast XNominalContrast ; YAxLims], ...
    'XNominalContrast', XNominalContrast, ...
    'YNominalContrast', YNominalContrast, ...
    'Color', [theRGB(1, :) ; theRGB(3, :)]);

% Save out the figures
set(figS, 'PaperPosition', [0 0 6 6]);
set(figS, 'PaperSize', [6 6]);
set(figS, 'Color', 'w');
set(figS, 'InvertHardcopy', 'off');
saveas(figS, fullfile(sstRoot, 'OO', 'plots', 'SSTReceptorDemo_SContrast.png'), 'png');