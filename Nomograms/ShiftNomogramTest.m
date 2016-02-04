%% ShiftNomogramTest.m
%
% Program to test shifting nomograms.
%
% 2/3/16    ms      Wrote it.

%% Define some parameters
S = [380 2 201]; wls = SToWls(S);

%% Stockman-Sharpe tabulated vs. nomogram
% Check accordance of tabulated Stockman-Sharpe pigment absorbance and
% fitted nomogram. The second panel corresponds to Fig. 12 in Stockman &
% Sharpe (2000).

% Get Stockman-Sharpe tabulated absorbance
load T_log10coneabsorbance_ss
T_StockmanSharpeAbsorbance = 10.^SplineCmf(S_log10coneabsorbance_ss,T_log10coneabsorbance_ss,S,2);

% Get nomogram absorbance
lambdaMaxNominalStockmanSharpeNomogram = [558.9 530.3 420.7]';
T_nomogramAbsorbance = StockmanSharpeNomogram(S,lambdaMaxNominalStockmanSharpeNomogram);

% Make sure nomogram peaks at the expected places
%
% Each of these checks should come back as 1, which they
% do
check1 = StockmanSharpeNomogram([558.9 1 1],558.9);
check2 = StockmanSharpeNomogram([530.3 1 1],530.3);
check3 = StockmanSharpeNomogram([420.7 1 1],420.7);
if (abs(check1-1) > 1e-6 |  abs(check2-1) > 1e-6 | abs(check3-1) > 1e-6)
    error('Failed basic check');
end

% Plot the absorbances
theLMSCols = [232 21 21 ; 33 143 51 ;  99 118 230]/255;
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
xlabel('Wavelength [nm]'); ylabel('Relative spectral sensitivity');
title({'Tabulated vs. nomogram absorbance' 'Linear'});

% Compare this one with Stockman & Sharpe (2000), Fig. 12.  The deviation
% looks a little larger here but their plot has those huge points which
% make it hard to see the differences.  Note that in our tabular S cone
% absorbance, we linear extrapolate at long wavelengths just to keep code
% from crashing when it doesn't find the wavelength sampling across cones.
%
% Note that the axis in Stockman & Sharpe (2000), Fig. 12, is in log
% wavelength.
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
xlim([375 785]); ylim([-8 1]);
xlabel('Wavelength [nm]'); ylabel('log spectral sensitivity');
title({'Tabulated vs. nomogram absorbance' 'Logarithmic'})

set(gcf, 'PaperPosition', [0 0 8 4]); % Position plot at left hand corner with width 8 and height 4.
set(gcf, 'PaperSize', [8 4]); % Set the paper to have width 8 and height 4.
saveas(gcf, 'StockmanSharpe_TabulatedVsNomogramLMSFundamentals.png', 'png');

%% Plot the densities as a function of wavenumber normalized to the maximum wavenumber
% Neitz & Neitz (2011, http://www.ncbi.nlm.nih.gov/pubmed/21167193) write:
% "The variability in normal pigments makes understanding the photopigment
% complement of color anomalous individuals more complicated; however, the
% discovery that the absorption spectra of primate photoreceptors assume a
% common shape when plotted on a normalized wave number axis (wave number
% divided by wave number of maximum sensitivity) (Baylor et al., 1987,
% Lamb, 1995, Mansfield, 1985), has greatly simplified the characterization
% of variant cone pigments. The implication is that all primate
% photopigments have a common shape that can be used to completely
% characterize the spectral properties of any human photopigment variant
% just from knowing its wavelength maximum."
%
% We are doing this here, namely convert to a wavenumber axis, normalize by
% maximum sensitivity.

% Wavenumber is given as the reciprocal of wavelength.
wavenumber = 1./wls;

% Find the indices of maximum sensitivity
[~, maxIdx] = max(T_StockmanSharpeAbsorbance, [], 2);

% Normalize the wave number
for ii = 1:3
    wavenumberNorm(ii, :) = wavenumber/wavenumber(maxIdx(ii));
end


% Normalize the wave number
for ii = 1:3
    logFrequencyNorm(ii, :) = log10(wls) - log10(wls(maxIdx(ii)));
end

figNormalizedWaveNum = figure;
% Normalized wave number
subplot(1, 2, 1);
hold on;
% Plot in this normalized axis
for ii = 1:3
    h1(ii) = plot(wavenumberNorm(ii, :), T_StockmanSharpeAbsorbance(ii, :), 'Color', theLMSCols(ii, :), 'LineWidth', 2);
end
plot([1 1], [-0.01 1.01], '--k');
pbaspect([1 1 1]);
legend(h1, 'L [tab.]', 'M [tab.]', 'S [tab.]'); legend boxoff;
set(gca, 'TickDir', 'out'); box off;
xlim([0.5 1.5]); ylim([-0.01 1.01]);
xlabel('Normalized wavenumber [nm^{-1}]'); ylabel('Relative spectral sensitivity');
title({'Stockman-Sharpe (2000) pigment absorbances' 'Normalized wavenumber representation'});

% Log frequency
subplot(1, 2, 2);
hold on;
% Plot in this normalized axis
for ii = 1:3
    h1(ii) = plot(logFrequencyNorm(ii, :), T_StockmanSharpeAbsorbance(ii, :), 'Color', theLMSCols(ii, :), 'LineWidth', 2);
end
plot([0 0], [-0.01 1.01], '--k');
pbaspect([1 1 1]);
legend(h1, 'L [tab.]', 'M [tab.]', 'S [tab.]'); legend boxoff;
set(gca, 'TickDir', 'out'); box off;
xlim([-0.3 0.3]); ylim([-0.01 1.01]);
xlabel('Normalized log wavelength [log nm]'); ylabel('Relative spectral sensitivity');
title({'Stockman-Sharpe (2000) pigment absorbances' 'Normalized log wavelength representation'});

set(gcf, 'PaperPosition', [0 0 8 4]); % Position plot at left hand corner with width 8 and height 4.
set(gcf, 'PaperSize', [8 4]); % Set the paper to have width 8 and height 4.
saveas(gcf, 'StockmanSharpe_NormalizedWaveNumberAndLogFrequency.png', 'png');

%% Shift the tabulated spectral sensitivities in a log-wavelength, log-sensitivity plane
lambdaMaxShift = -30;
for ii = 1:3
    log10_T_StockmanSharpeAbsorbance_Shifted(ii, :) = ShiftFundamental(wls, log10(T_StockmanSharpeAbsorbance(ii, :)), lambdaMaxShift); % 2 nm shift
    [~, maxIdx1] = max(log10(T_StockmanSharpeAbsorbance(ii, :)));
    [~, maxIdx2] = max(log10_T_StockmanSharpeAbsorbance_Shifted(ii, :));
    fprintf('\n');
    fprintf('\t*** lambda-max shift: %.2f nm\n', lambdaMaxShift);
    fprintf('\t>>> Old lambda-max: %.2f nm\n', wls(maxIdx1));
    fprintf('\t>>> New lambda-max: %.2f nm\n', wls(maxIdx2));
end

for ii = 1:3
    subplot(1, 3, ii);
    hold on;
    
    % Plot the unshifted cone fundamental
    h1(ii) = plot(wls, log10(T_StockmanSharpeAbsorbance(ii, :)), '-', 'Color', 'k', 'LineWidth', 2);
    
    % Plot the shifted cone fundamental
    h2(ii) = plot(wls, log10_T_StockmanSharpeAbsorbance_Shifted(ii, :), 'Color', theLMSCols(ii, :), 'LineWidth', 2);
    
    xlabel('Normalized log wavelength [log nm]'); ylabel('Relative spectral sensitivity');
    set(gca, 'TickDir', 'out'); box off;
    xlim([375 785]); ylim([-8 1]);
    xlabel('Wavelength [nm]'); ylabel('log spectral sensitivity');
    title({'Shifted spectral sensitivities' [num2str(lambdaMaxShift) ' nm']});
    legend([h1(ii) h2(ii)], 'Unshifted', 'Shifted', 'Location', 'SouthWest'); legend boxoff;
    pbaspect([1 1 1]);
end


set(gcf, 'PaperPosition', [0 0 9 3]); % Position plot at left hand corner with width 8 and height 4.
set(gcf, 'PaperSize', [9 3]); % Set the paper to have width 8 and height 4.
saveas(gcf, 'StockmanSharpe_ShiftTabulatedTest.png', 'png');