tbUse('SilentSubstitutionToolbox');
%%
receptorObj = SSTReceptorHuman('obsAgeYrs', 30);

%%
theRGB = DefaultReceptorColors;


%%
% Load data
tmp = load('/Users/spitschan/Documents/MATLAB/toolboxes/SilentSubstitutionToolbox/ContrastSplatter/ContrastSplatterDemoData/spd_contrastsplatterdemo_bg.mat');
bgSpd = tmp.spd;
tmp = load('/Users/spitschan/Documents/MATLAB/toolboxes/SilentSubstitutionToolbox/ContrastSplatter/ContrastSplatterDemoData/spd_contrastsplatterdemo_mod.mat');
modSpd = tmp.spd;

%% Parametric variation
NTitrations = 50;
receptorObj.makeSpectralSensitivitiesParametricVariation('WhichParameter', 'dlens', 'NTitrations', NTitrations);

% Calculate contrast
for ii = 1:size(receptorObj.Ts, 2)
    T_receptors = receptorObj.Ts{ii}.T_energyNormalized;
    for jj = 1:size(receptorObj.Ts{ii}.T_energyNormalized, 1)
        contrasts(jj, ii) = (T_receptors(jj, :)*(modSpd-bgSpd))./(T_receptors(jj, :)*bgSpd);
    end
    lmContrast(:, ii) = [1 1 0]' \ contrasts(:, ii);
    postRecepContrasts(:, ii) = [1 1 1 ; 1 -1 0 ; 0 0 1]' \ contrasts(:, ii);
end


% Plot parametric variations
for ii = 1:size(contrasts, 1)
    plot(contrasts(ii, :)', '-o', 'Color', theRGB(ii, :)); hold on;
end

% Set up all other parameters

%% Stochastic sampling
NSamples = 1000;
receptorObj.makeSpectralSensitivitiesStochastic('NSamples', NSamples);

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
    'Color', [rgbL ; rgbM]);

set(figLM, 'PaperPosition', [0 0 6 6]);
set(figLM, 'PaperSize', [6 6]);
set(gcf, 'Color', 'w');
set(gcf, 'InvertHardcopy', 'off');
saveas(figLM, 'LMContrast.png', 'png');
saveas(figLM, 'LMContrast.pdf', 'pdf');

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
    'Color', [rgbL ; rgbM]);

set(figS, 'PaperPosition', [0 0 6 6]);
set(figS, 'PaperSize', [6 6]);
set(gcf, 'Color', 'w');
set(gcf, 'InvertHardcopy', 'off');
saveas(figS, 'SContrast.png', 'png');