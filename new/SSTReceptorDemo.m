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
    end
end

axLims = [-0.08 0.08];
subplot(1, 3, 1);
% Reference lines
hold on;
plot([0 0], axLims, '-k', axLims, [0 0], '-k');
plot(contrasts(1, :), contrasts(2, :), '.k');
xlim(axLims); ylim(axLims);
pbaspect([1 1 1]);
box off;
set(gca, 'TickDir', 'out');
xlabel('L cone contrast'); ylabel('M cone contrast');
e2 = FitEllipse(contrasts([1:2], :)', 2);

% Fit ellipse and plot
plot(e2(1,:), e2(2,:), '-r', 'LineWidth', 2);



%
subplot(1, 3, 2);
hold on;
plot([0 0], axLims, '-k', axLims, [0 0], '-k');
plot(contrasts(2, :), contrasts(3, :),'.k');
xlim(axLims); ylim(axLims);
pbaspect([1 1 1]);
box off;
set(gca, 'TickDir', 'out');
xlabel('M cone contrast'); ylabel('S cone contrast');

% Fit ellipse and plot
e2 = FitEllipse(contrasts([2:3], :)', 2);

% Fit ellipse and plot
plot(e2(1,:), e2(2,:), '-r', 'LineWidth', 2);


subplot(1, 3, 3);
hold on;
plot([0 0], axLims, '-k', axLims, [0 0], '-k');
plot(contrasts(1, :), contrasts(3, :), '.k');
xlim(axLims); ylim(axLims);
pbaspect([1 1 1]);
box off;
set(gca, 'TickDir', 'out');
xlabel('L cone contrast'); ylabel('S cone contrast');

% Fit ellipse and plot
e2 = FitEllipse(contrasts([1 3], :)', 2);

% Fit ellipse and plot
plot(e2(1,:), e2(2,:), '-r', 'LineWidth', 2);