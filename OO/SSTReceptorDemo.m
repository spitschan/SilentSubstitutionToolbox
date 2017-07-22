% SSTReceptorDemo.m
%
% This script demonstrates the use of the SSTReceptor object.
%
% 7/22/17   ms      Commented.

%% Blank slate
clearvars; close all; clc;

%% Set up paths
r = tbUse('SilentSubstitutionToolbox');

% Infer the SST root
sstRoot = tbLocateToolbox('SilentSubstitutionToolbox');

%% Get some colors
theRGB = DefaultReceptorColors;

%% Load backgound and modulation spectra from some demo data in SST
tmp = load(fullfile(sstRoot, 'ContrastSplatter/ContrastSplatterDemoData/spd_contrastsplatterdemo_bg.mat'));
bgSpd = tmp.spd;
tmp = load(fullfile(sstRoot, 'ContrastSplatter/ContrastSplatterDemoData/spd_contrastsplatterdemo_mod.mat'));
modSpd = tmp.spd;

%% Set up receptor object
% If the receptor object has been pre-cached (this is sometimes useful for
% large resampling N), see if it exists.
cacheDir = fullfile(fileparts(which('SSTReceptorHuman')), 'cache');

% We're pulling out one specific pre-cached receptor obkect
theDesiredHash = '50be1fa91c9911cb5be28c2e2f3adad3';

% Determine the path
theFile = fullfile(cacheDir, ['SSTReceptorHuman_' theDesiredHash '.mat']);
% See if the desired hash exists in the cache dir
if exist(theFile, 'file')
    % Load the object if it already exists
    load(theFile);
else
    % If it does not exist, make the receptor object
    receptorObj = SSTReceptorHuman('verbosity', 'high', 'obsAgeYrs', 32);
    
    % Set the number of titrations per parameter
    NTitrations = 16;
    
    % Do the parametric "sampling"
    theIndDiffParams = {'dlens' 'dmac' 'dphotopigmentL' 'dphotopigmentM' 'dphotopigmentS' ...
        'lambdaMaxShiftL' 'lambdaMaxShiftM' 'lambdaMaxShiftS' 'obsPupilDiameterMm'};
    for ss = 1:length(theIndDiffParams)
        % Vary the parameter
        [~, parv, parvlabel, parvlabellong, parvreal] = makeSpectralSensitivitiesParametricVariation(receptorObj, ...
            'WhichParameter', theIndDiffParams{ss}, 'NTitrations', NTitrations);
    end
    
    % Stochastic resampling
    NSamples = 1000;
    receptorObj.makeSpectralSensitivitiesStochastic('NSamples', NSamples);
    
    % Save the receptor object
    receptorObj.setMD5Hash();
    save(fullfile(cacheDir, ['SSTReceptorHuman_' receptorObj.MD5Hash]), 'receptorObj');
end

%% Plot the parametric varation
figParv = figure;
yAxLims = [-0.06 0.06];
for ss = 1:length(theIndDiffParams)
    if parv(1) < 0
        xAxLims = [parv(1)*1.1 parv(end)*1.1];
    else
        xAxLims = [parv(1)*0.9 parv(end)*1.1];
    end
    
    % Calculate contrast
    for ii = 1:size(receptorObj.Ts, 2)
        T_receptors = receptorObj.Ts{ii}.T_energyNormalized;
        for jj = 1:size(receptorObj.Ts{ii}.T_energyNormalized, 1)
            contrasts(jj, ii) = (T_receptors(jj, :)*(modSpd-bgSpd))./(T_receptors(jj, :)*bgSpd);
        end
        lmContrast(:, ii) = [1 1 0]' \ contrasts(:, ii);
        postRecepContrasts(:, ii) = [1 1 1 ; 1 -1 0 ; 0 0 1]' \ contrasts(:, ii);
    end
    
    subplot(1, length(theIndDiffParams), ss);
    hold on;
    plot(xAxLims, [0 0], ':k');
    % Plot parametric variations
    for ii = 1:size(contrasts, 1)
        plot(parv, contrasts(ii, :)', '-o', 'Color', theRGB(ii, :), 'MarkerEdgeColor', 'k', ...
            'MarkerFaceColor', theRGB(ii, :)); hold on;
    end
    
    % Add title and tweak plots
    title(parvlabellong);
    box off;
    xlabel(parvlabel);
    xlim(xAxLims);
    ylim(yAxLims);
    pbaspect([1 1 1]);
    set(gca, 'TickDir', 'out');
end
set(figParv, 'PaperPosition', [0 0 13 3.5]);
set(figParv, 'PaperSize', [13 3.5]);
set(figParv, 'Color', 'w');
set(figParv, 'InvertHardcopy', 'off');
saveas(figParv, 'ParvFig.png', 'png');


%%
% Calculate contrast
for ii = 1:NSamples
    T_receptors = receptorObj.Ts{ii}.T_energyNormalized;
    for jj = 1:size(receptorObj.Ts{ii}.T_energyNormalized, 1)
        contrasts(jj, ii) = (T_receptors(jj, :)*(modSpd-bgSpd))./(T_receptors(jj, :)*bgSpd);
    end
    lmContrast(:, ii) = [1 1 0]' \ contrasts(:, ii);
    postRecepContrasts(:, ii) = [1 1 1 ; 1 -1 0 ; 0 0 1]' \ contrasts(:, ii);
end


axLims = [-0.03 0.03];

%% L vs. M
figLM = figure;
XNominalContrast = 300;
YNominalContrast = 0;
XAxLims = XNominalContrast+[-0.02 0.02];
YAxLims = YNominalContrast+[-0.02 0.02];

ScatterplotWithHistogram(300+contrasts(1, :), contrasts(2, :), ...
    'XLim', XAxLims, 'YLim', YAxLims, 'BinWidth', 0.004, ...
    'XLabel', 'L contrast', 'YLabel', 'M contrast', ...
    'XRefLines', [XAxLims ; YNominalContrast YNominalContrast], ...
    'YRefLines', [XNominalContrast XNominalContrast ; YAxLims], ...
    'XNominalContrast', XNominalContrast, ...
    'YNominalContrast', YNominalContrast, ...
    'Color', [theRGB(1, :) ; theRGB(2, :)]);

set(figLM, 'PaperPosition', [0 0 6 6]);
set(figLM, 'PaperSize', [6 6]);
set(figLM, 'Color', 'w');
set(figLM, 'InvertHardcopy', 'off');
saveas(figLM, 'LMContrast.png', 'png');

%% L+M vs. S
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

set(figS, 'PaperPosition', [0 0 6 6]);
set(figS, 'PaperSize', [6 6]);
set(figS, 'Color', 'w');
set(figS, 'InvertHardcopy', 'off');
saveas(figS, 'SContrast.png', 'png');