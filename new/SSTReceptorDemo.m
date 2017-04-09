tbUse('SilentSubstitutionToolbox');
clearvars; close all; clc;

%% Get some colors
theRGB = DefaultReceptorColors;

%% Load data
tmp = load('/Users/spitschan/Documents/MATLAB/toolboxes/SilentSubstitutionToolbox/ContrastSplatter/ContrastSplatterDemoData/spd_contrastsplatterdemo_bg.mat');
bgSpd = tmp.spd;
tmp = load('/Users/spitschan/Documents/MATLAB/toolboxes/SilentSubstitutionToolbox/ContrastSplatter/ContrastSplatterDemoData/spd_contrastsplatterdemo_mod.mat');
modSpd = tmp.spd;

%% Set up receptor object
receptorObj = SSTReceptorHuman('obsAgeYrs', 32);

%% Generate the space of cones 
NTitrations = 16;

% Do the parametric "sampling"
theIndDiffParams = {'dlens' 'dmac' 'dphotopigmentL' 'dphotopigmentM' 'dphotopigmentS' ...
    'lambdaMaxShiftL' 'lambdaMaxShiftM' 'lambdaMaxShiftS' 'obsPupilDiameterMm'};
for ss = 1:length(theIndDiffParams)
    % Vary the parameter
    [~, parv, parvlabel, parvlabellong, parvreal] = makeSpectralSensitivitiesParametricVariation(receptorObj, ...
        'WhichParameter', theIndDiffParams{ss}, 'NTitrations', NTitrations);
end

% Stochastic sampling
NSamples = 1000;
receptorObj.makeSpectralSensitivitiesStochastic('NSamples', NSamples);

% Save the receptor object
receptorObj.MD5Hash = DataHash(receptorObj);
outSaveDir = fullfile(which('SSTReceptorHuman'), 'cache');
save(fullfile(outSaveDir, ['SSTReceptor_' receptorObj.MD5Hash], 'receptorObj'));

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