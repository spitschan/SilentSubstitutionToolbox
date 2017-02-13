tbUse('SilentSubstitutionToolbox');

receptorObj = SSTReceptorHuman('obsAgeYrs', 30);
NSamples = 1000;
receptorObj.makeSpectralSensitivitiesStochastic('NSamples', NSamples);

tmp = load('/Users/spitschan/Documents/MATLAB/toolboxes/SilentSubstitutionToolbox/ContrastSplatter/ContrastSplatterDemoData/spd_contrastsplatterdemo_bg.mat');
bgSpd = tmp.spd;
tmp = load('/Users/spitschan/Documents/MATLAB/toolboxes/SilentSubstitutionToolbox/ContrastSplatter/ContrastSplatterDemoData/spd_contrastsplatterdemo_mod.mat');
modSpd = tmp.spd;

%%
% Calculate contrast
for ii = 1:NSamples
    T_receptors = receptorObj.Ts{ii}.T_energyNormalized;
    for jj = 1:size(receptorObj.Ts{ii}.T_energyNormalized, 1)
        contrasts(jj, ii) = (T_receptors(jj, :)*(modSpd-bgSpd))./(T_receptors(jj, :)*bgSpd);
        lmContrast(:, ii) = [1 1 0]' \ contrasts(:, ii);
        postRecepContrasts(:, ii) = [1 1 1 ; 1 -1 0 ; 0 0 1]' \ contrasts(:, ii);
    end
end

%%
rgbL = [0.8941    0.1020    0.1098];
rgbM = [0.3020    0.6863    0.2902];
rgbS = [0.2157    0.4941    0.7216];

axLims = [-0.03 0.03];

%%
figLM = figure;
XAxLims = [-0.025 0.025];
YAxLims = [-0.025 0.025];

ScatterplotWithHistogram(contrasts(1, :), contrasts(2, :), ...
    'XLim', XAxLims, 'YLim', YAxLims, 'BinWidth', 0.004, ...
    'XLabel', 'L contrast', 'YLabel', 'M contrast', ...
    'XRefLines', [XAxLims ; 0 0], ...
    'YRefLines', [0 0 ; YAxLims], ...
    'Color', [rgbL ; rgbM]);

set(figLM, 'PaperPosition', [0 0 6 6]);
set(figLM, 'PaperSize', [6 6]);
saveas(figLM, 'LMContrast.png', 'png');

%%
figS = figure;
XAxLims = [-0.02 0.02];
YAxLims = [-0.01 0.08];

ScatterplotWithHistogram(lmContrast(1, :), contrasts(3, :), ...
    'XLim', XAxLims, 'YLim', YAxLims, 'BinWidth', 0.004, ...
    'XLabel', 'L+M contrast', 'YLabel', 'S contrast', ...
    'XRefLines', [XAxLims ; 0 0], ...
    'YRefLines', [0 0 ; YAxLims], ...
    'Color', [0.5 0.5 0.5 ; rgbS]);

set(figS, 'PaperPosition', [0 0 6 6]);
set(figS, 'PaperSize', [6 6]);
saveas(figS, 'SContrast.png', 'png');

!open SContrast.png