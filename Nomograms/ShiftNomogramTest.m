%% ShiftNomogramTest.m
%
% Program to test shifting nomograms.
%
% 2/3/16    ms      Wrote it.

%% Define some parameters
S = [380 2 201]; wls = SToWls(S);

%% Stockman-Sharpe tabulated vs. nomogram
% Check accordance of tabulated Stockman-Sharpe pigment absorbance and
% fitted nomogram

% Get Stockman-Sharpe tabulated absorbance
load T_log10coneabsorbance_ss
T_StockmanSharpeAbsorbance = 10.^SplineCmf(S_log10coneabsorbance_ss,T_log10coneabsorbance_ss,S,2);

% Get nomogram absorbance
lambdaMaxNominalStockmanSharpeNomogram = [558.9 530.3 420.7]';
T_nomogramAbsorbance = StockmanSharpeNomogram(S,lambdaMaxNominalStockmanSharpeNomogram);

%%
theLMSCols = [232 21 21 ; 33 143 51 ;  99 118 230]/255
figTabulatedVsNomogram = figure;
subplot(1, 2, 1);
hold on;
for ii = 1:3
    plot(wls, T_StockmanSharpeAbsorbance(ii, :), 'Color', theLMSCols(ii, :), 'LineWidth', 3);
    plot(wls, T_nomogramAbsorbance(ii, :), 'Color', 'k', 'LineWidth', 1.5);
end
pbaspect([1 1 1]);
set(gca, 'TickDir', 'out'); box off;
xlim([375 785]); ylim([-0.01 1.01]);
xlabel('Wavelength [nm]'); ylabel('Relative sensitivity');
title({'Tabulated vs. nomogram absorbance' 'Linear'});

subplot(1, 2, 2);
hold on;
for ii = 1:3
    h1(ii) = plot(wls, log10(T_StockmanSharpeAbsorbance(ii, :)), 'Color', theLMSCols(ii, :), 'LineWidth', 3);
    h2(ii) = plot(wls, log10(T_nomogramAbsorbance(ii, :)), 'Color', 'k', 'LineWidth', 1.5);
end
pbaspect([1 1 1]);
legend([h1 h2], 'L [tab.]', 'M [tab.]', 'S [tab.]', 'L [nomogram]', 'M [nomogram]', 'S [nomogram]', ...
    'Location', 'SouthWest'); legend boxoff;
set(gca, 'TickDir', 'out'); box off;
xlim([375 785]); ylim([-10 1]);
xlabel('Wavelength [nm]'); ylabel('log relative sensitivity');
title({'Tabulated vs. nomogram absorbance' 'Logarithmic'})

set(gcf, 'PaperPosition', [0 0 8 4]); %Position plot at left hand corner with width 8 and height 5.
set(gcf, 'PaperSize', [8 4]); %Set the paper to have width 8 and height 5.
saveas(gcf, 'TabulatedVsNomogramLMSFundamentals.png', 'png');
