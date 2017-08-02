% Set up parameters for the OneLight stuff
calPath = fullfile('/Users/spitschan/Documents/MATLAB/Toolboxes/SilentSubstitutionToolbox', 'ReceptorIsolate', 'ReceptorIsolateDemoData', []);
cal = LoadCalFile('OneLightDemoCal.mat',[],calPath);
S = cal.describe.S; wls = SToWls(S);
B_primary = cal.computed.pr650M;
ambientSpd = cal.computed.pr650MeanDark;
backgroundPrimary = 0.5*ones(size(B_primary,2),1);
whichPrimariesToPin = [];
primaryHeadRoom = 0.02;
maxPowerDiff = 10^-1.5;
whichReceptorsToTarget = [5];
whichReceptorsToIgnore = [];
whichReceptorsToMinimize = [];
desiredContrast = [];

% Make the receptors
T_receptors = receptorObj.T.T_energyNormalized;

% Do the optimization
modulationPrimary = ReceptorIsolate(T_receptors,whichReceptorsToTarget, whichReceptorsToIgnore, whichReceptorsToMinimize, ...
  B_primary, backgroundPrimary, backgroundPrimary, whichPrimariesToPin,...
  primaryHeadRoom, maxPowerDiff, desiredContrast, ambientSpd);

% Extract the background spd
backgroundSpd = B_primary*backgroundPrimary;
modulationSpd = B_primary*modulationPrimary;
T_receptors*(modulationSpd-backgroundSpd) ./ (T_receptors*backgroundSpd)

%% Iteration 0 
backgroundSpd0 = backgroundSpd;
modulationSpd0 = modulationSpd;
contrast0 = T_receptors*(backgroundSpd-backgroundSpd) ./ (T_receptors*backgroundSpd)
contrast4 = T_receptors*(modulationSpd-backgroundSpd) ./ (T_receptors*backgroundSpd)

%% Iteration 1
modulationSpd1 = modulationSpd;
modulationSpd1(1:42) = backgroundSpd(1:42);
modulationSpd1(79:end) = backgroundSpd(79:end);
contrast1 = T_receptors*(modulationSpd1-backgroundSpd) ./ (T_receptors*backgroundSpd)

%% Iteration 2
modulationSpd2 = modulationSpd;
modulationSpd2(79:end) = backgroundSpd(79:end);
contrast2 = T_receptors*(modulationSpd2-backgroundSpd) ./ (T_receptors*backgroundSpd)

%% Iteration 3
modulationSpd3 = modulationSpd;
modulationSpd3(105:end) = backgroundSpd(105:end);
contrast3 = T_receptors*(modulationSpd3-backgroundSpd) ./ (T_receptors*backgroundSpd)

%% Plot
theRGB = DefaultReceptorColors;

% Iteration 0
subplot(1, 2, 1);
plot(wls, backgroundSpd, '-k', 'LineWidth', 1.3); hold on;
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; xlim([380 780]); ylim([0 0.1]);
xlabel('Wavelength [nm]'); ylabel('Radiance'); hold off; hold off;

subplot(1, 2, 2);
for ii = 1:4
    h(ii) = bar(ii, 100*contrast0(ii)); hold on;
    set(h(ii), 'FaceColor', theRGB(ii, :), 'EdgeColor', theRGB(ii, :)); 
end
ylabel('Contrast [%]');
set(gca, 'XTick', 1:4, 'XTickLabel', {'L', 'M', 'S', 'Mel'});
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; xlim([0 5]); ylim([-70 70]);
pause(0.8); hold off;

% Save out the figure
set(gcf, 'PaperPosition', [0 0 6 4]);
set(gcf, 'PaperSize', [6 4]);
set(gcf, 'Color', 'w');
set(gcf, 'InvertHardcopy', 'off');
saveas(gcf, 'Explainer0.pdf', 'pdf');

% Iteration 0
subplot(1, 2, 1);
plot(wls, backgroundSpd, ':k', 'LineWidth', 1.3); hold on; plot(wls, modulationSpd1, '-r', 'LineWidth', 1.3);
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; xlim([380 780]); ylim([0 0.1]);
xlabel('Wavelength [nm]'); ylabel('Radiance'); hold off;

subplot(1, 2, 2);
for ii = 1:4
    h(ii) = bar(ii, 100*contrast1(ii)); hold on;
    set(h(ii), 'FaceColor', theRGB(ii, :), 'EdgeColor', theRGB(ii, :)); 
end
ylabel('Contrast [%]');
set(gca, 'XTick', 1:4, 'XTickLabel', {'L', 'M', 'S', 'Mel'});
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; xlim([0 5]); ylim([-70 70]);
pause(0.8); hold off;

% Save out the figure
set(gcf, 'PaperPosition', [0 0 6 4]);
set(gcf, 'PaperSize', [6 4]);
set(gcf, 'Color', 'w');
set(gcf, 'InvertHardcopy', 'off');
saveas(gcf, 'Explainer1.pdf', 'pdf');

% Iteration 1
subplot(1, 2, 1);
plot(wls, backgroundSpd, ':k', 'LineWidth', 1.3); hold on; plot(wls, modulationSpd2, '-r', 'LineWidth', 1.3);
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; xlim([380 780]); ylim([0 0.1]);
xlabel('Wavelength [nm]'); ylabel('Radiance'); hold off;

subplot(1, 2, 2);
for ii = 1:4
    h(ii) = bar(ii, 100*contrast2(ii)); hold on;
    set(h(ii), 'FaceColor', theRGB(ii, :), 'EdgeColor', theRGB(ii, :)); 
end
ylabel('Contrast [%]');
set(gca, 'XTick', 1:4, 'XTickLabel', {'L', 'M', 'S', 'Mel'});
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; xlim([0 5]); ylim([-70 70]);
pause(0.8); hold off;

% Save out the figure
set(gcf, 'PaperPosition', [0 0 6 4]);
set(gcf, 'PaperSize', [6 4]);
set(gcf, 'Color', 'w');
set(gcf, 'InvertHardcopy', 'off');
saveas(gcf, 'Explainer2.pdf', 'pdf');

% Iteration 2
subplot(1, 2, 1);
plot(wls, backgroundSpd, ':k', 'LineWidth', 1.3); hold on; plot(wls, modulationSpd3, '-r', 'LineWidth', 1.3);
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; xlim([380 780]); ylim([0 0.1]);
xlabel('Wavelength [nm]'); ylabel('Radiance'); hold off;

subplot(1, 2, 2);
for ii = 1:4
    h(ii) = bar(ii, 100*contrast3(ii)); hold on;
    set(h(ii), 'FaceColor', theRGB(ii, :), 'EdgeColor', theRGB(ii, :)); 
end
ylabel('Contrast [%]');
set(gca, 'XTick', 1:4, 'XTickLabel', {'L', 'M', 'S', 'Mel'});
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; xlim([0 5]); ylim([-70 70]);
pause(0.8); hold off;

% Save out the figure
set(gcf, 'PaperPosition', [0 0 6 4]);
set(gcf, 'PaperSize', [6 4]);
set(gcf, 'Color', 'w');
set(gcf, 'InvertHardcopy', 'off');
saveas(gcf, 'Explainer3.pdf', 'pdf');

% Iteration 3
subplot(1, 2, 1);
plot(wls, backgroundSpd, ':k', 'LineWidth', 1.3); hold on; plot(wls, modulationSpd, '-r', 'LineWidth', 1.3);
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; xlim([380 780]); ylim([0 0.1]);
xlabel('Wavelength [nm]'); ylabel('Radiance'); hold off;

subplot(1, 2, 2);
for ii = 1:4
    h(ii) = bar(ii, 100*contrast4(ii)); hold on;
    set(h(ii), 'FaceColor', theRGB(ii, :), 'EdgeColor', theRGB(ii, :)); 
end
ylabel('Contrast [%]');
set(gca, 'XTick', 1:4, 'XTickLabel', {'L', 'M', 'S', 'Mel'});
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off; xlim([0 5]); ylim([-70 70]);
pause(0.8); hold off

% Save out the figure
set(gcf, 'PaperPosition', [0 0 6 4]);
set(gcf, 'PaperSize', [6 4]);
set(gcf, 'Color', 'w');
set(gcf, 'InvertHardcopy', 'off');
saveas(gcf, 'Explainer4.pdf', 'pdf');