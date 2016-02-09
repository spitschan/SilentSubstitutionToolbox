function ShiftNomogramTest
% ShiftNomogramTest
%
% Program to test shifting nomograms and
% related.
%
% 2/3/16    ms      Wrote it.

%% Initialize
close all;

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
xlabel('Wavelength [nm]'); ylabel('Photopigment absorbance');
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
xlabel('Wavelength [nm]'); ylabel('Log photopigment absorbance');
title({'Tabulated vs. nomogram absorbance' 'Logarithmic'})

set(gcf, 'PaperPosition', [0 0 8 4]); % Position plot at left hand corner with width 8 and height 4.
set(gcf, 'PaperSize', [8 4]);         % Set the paper to have width 8 and height 4.
saveas(gcf, 'StockmanSharpe_TabulatedVsNomogramLMSFundamentals.png', 'png');

%% Plot the densities as a function of wavenumber normalized to the maximum wavenumber
%
% Also try log wavenumber normalization.
%
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
%
% Stockman & Sharpe (2000), however prefer the log wavelength
% normalization.

% Wavenumber is given as the reciprocal of wavelength.
wavenumber = 1./wls;

% Find the indices of maximum sensitivity
[~, maxIdx] = max(T_StockmanSharpeAbsorbance, [], 2);

% Normalize the wave number and log wavelength
for ii = 1:3
    wavenumberNorm(ii, :) = wavenumber/wavenumber(maxIdx(ii));
    logWavelengthNorm(ii, :) = log10(wls) - log10(wls(maxIdx(ii)));
end

% Plot normalized pigment absorbances
%
% First wavenumber
figNormalizedWaveNum = figure;
subplot(1, 2, 1);
hold on;
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

% Then log wavelength
subplot(1, 2, 2);
hold on;
for ii = 1:3
    h1(ii) = plot(logWavelengthNorm(ii, :), T_StockmanSharpeAbsorbance(ii, :), 'Color', theLMSCols(ii, :), 'LineWidth', 2);
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
saveas(gcf, 'StockmanSharpe_NormalizedWaveNumberAndLogWavelength.png', 'png');

%% Shift the tabulated spectral absorbances in a log-wavelength, log-sensitivity plane
lambdaMaxShift = 30;
log10_T_StockmanSharpeAbsorbance_Shifted = ShiftPhotopigmentAbsorbance(wls,log10(T_StockmanSharpeAbsorbance),[lambdaMaxShift lambdaMaxShift lambdaMaxShift]);
for ii = 1:3
    [~, maxIdx1] = max(log10(T_StockmanSharpeAbsorbance(ii, :)));
    [~, maxIdx2] = max(log10_T_StockmanSharpeAbsorbance_Shifted(ii, :));
    fprintf('\n');
    fprintf('\t*** lambda-max shift: %.2f nm\n', lambdaMaxShift);
    fprintf('\t>>> Old lambda-max: %.2f nm\n', wls(maxIdx1));
    fprintf('\t>>> New lambda-max: %.2f nm\n', wls(maxIdx2));
end

% Plot to check shifts
figShifted = figure;
for ii = 1:3
    
    subplot(1, 3, ii);
    hold on;
    
    % Plot the unshifted cone absorbance
    h1(ii) = plot(wls, log10(T_StockmanSharpeAbsorbance(ii, :)), '-', 'Color', 'k', 'LineWidth', 2);
    
    % Plot the shifted cone absorbance
    h2(ii) = plot(wls, log10_T_StockmanSharpeAbsorbance_Shifted(ii, :), 'Color', theLMSCols(ii, :), 'LineWidth', 2);
    
    xlabel('Normalized log wavelength [log nm]'); ylabel('Pigment absorbance');
    set(gca, 'TickDir', 'out'); box off;
    xlim([375 785]); ylim([-8 1]);
    xlabel('Wavelength [nm]'); ylabel('log spectral sensitivity');
    title({'Shifted photopigment absorbances' [num2str(lambdaMaxShift) ' nm']});
    legend([h1(ii) h2(ii)], 'Unshifted', 'Shifted', 'Location', 'SouthWest'); legend boxoff;
    pbaspect([1 1 1]);
end
set(gcf, 'PaperPosition', [0 0 9 3]); % Position plot at left hand corner with width 8 and height 4.
set(gcf, 'PaperSize', [9 3]); % Set the paper to have width 8 and height 4.
saveas(figShifted, 'StockmanSharpe_ShiftTabulatedTest.png', 'png');


%% Test if the shifting also works in our wrapper routine GetHumanPhotoreceptorSS
S = [380 2 201]; wls = SToWls(S);
lambdaMaxShift = [-5 10 5]; % Example values
T_energyNormalized1 = GetHumanPhotoreceptorSS(S, {'LConeTabulatedAbsorbance' 'MConeTabulatedAbsorbance', 'SConeTabulatedAbsorbance'}, 30, 32, 6, [0 0 0], [], [], []);
T_energyNormalized2 = GetHumanPhotoreceptorSS(S, {'LConeTabulatedAbsorbance' 'MConeTabulatedAbsorbance', 'SConeTabulatedAbsorbance'}, 30, 32, 6, lambdaMaxShift, [], [], []);

for ii = 1:3
    [~, maxIdx1] = max(log10(T_energyNormalized1(ii, :)));
    [~, maxIdx2] = max(T_energyNormalized2(ii, :));
    fprintf('\n');
    fprintf('\t*** lambda-max shift: %.2f nm\n', lambdaMaxShift(ii));
    fprintf('\t>>> Old lambda-max [after filtering]: %.2f nm\n', wls(maxIdx1));
    fprintf('\t>>> New lambda-max [after filtering]: %.2f nm\n', wls(maxIdx2));
end
fprintf('*** NOTE THAT DUE TO PRE-RECEPTORAL FILTERING, THESE VALUES WILL NOT BE CORRESPOND TO THE SHIFT IN THE PEAK ABSORBANCE. ***\n');

%% Finally, let's see if we get the 10° Stockman-Sharpe fundamentals when we pass in the right parameters
S = [380 2 201]; 
S = WlsToS((390:5:780)');
wls = SToWls(S);
[T_energyNormalized_SST, T_quantalNormalized_SST] = GetHumanPhotoreceptorSS(S, {'LConeTabulatedAbsorbance' 'MConeTabulatedAbsorbance', 'SConeTabulatedAbsorbance'}, 10, 32, 3, [0 0 0], [], [], []);
T_quantal_PTB = ComputeCIEConeFundamentals(S,10,32,3);
T_energy_PTB = EnergyToQuanta(S,T_quantal_PTB')';
load T_cones_ss10
T_cones_ss10_spline = SplineCmf(S_cones_ss10,T_cones_ss10,S,2);
T_cones_ss10_spline_quantal = QuantaToEnergy(S,T_cones_ss10_spline')';

for ii = 1:3
    T_energyNormalized_PTB(ii, :) = T_energy_PTB(ii, :)/max(T_energy_PTB(ii, :));
    T_quantalNormalized_SST(ii, :) = T_quantalNormalized_SST(ii, :)/max(T_quantalNormalized_SST(ii, :));
    T_cones_ss10_spline(ii,:) = T_cones_ss10_spline(ii,:)/max(T_cones_ss10_spline(ii,:));
    T_cones_ss10_spline_quantal(ii,:) = T_cones_ss10_spline_quantal(ii,:)/max(T_cones_ss10_spline_quantal(ii,:));
end
figure;
hold on;
plot(wls, T_energyNormalized_SST, '-k');
plot(wls, T_energyNormalized_PTB, '-b');
plot(SToWls(S_cones_ss10), T_cones_ss10, '-r');
plot(wls, T_cones_ss10_spline, '-g');

% See if we agree with CIEConeFundamentalsTest.m
figure;

% Do exact what is done in CIEConeFundamentalsTest
targetRaw = load('T_cones_ss10');
T_targetEnergy = SplineCmf(targetRaw.S_cones_ss10,targetRaw.T_cones_ss10,S,2);
T_targetQuantal10 = QuantaToEnergy(S,T_targetEnergy')';

for i = 1:3
    T_targetQuantal10(i,:) = T_targetQuantal10(i,:)/max(T_targetQuantal10(i,:));
end
T_predictQuantalCIE10 = ComputeCIEConeFundamentals(S,10,32,3);

subplot(1, 2, 1);
hold on;
for ii = 1:3
plot(wls, (T_targetQuantal10(ii, :)-T_predictQuantalCIE10(ii, :))', 'Color', theLMSCols(ii, :))
end
xlim([380 780]); ylim([-1.5e-3 1.5e-3]);
pbaspect([1 1 1]);
xlabel('Wavelength [nm]');
ylabel('\Delta');
title('CIEConeFundamentalsTest.m');

subplot(1, 2, 2);
hold on;
for ii = 1:3
plot(wls, (T_cones_ss10_spline_quantal(ii, :)-T_quantalNormalized_SST(ii, :))', 'Color', theLMSCols(ii, :))
end
xlim([380 780]); ylim([-1.5e-3 1.5e-3]);
pbaspect([1 1 1]);
xlabel('Wavelength [nm]');
ylabel('\Delta');
title('ShiftNomogramTest.m');

set(gcf, 'PaperPosition', [0 0 8 4]); % Position plot at left hand corner with width 8 and height 4.
set(gcf, 'PaperSize', [8 4]); % Set the paper to have width 8 and height 4.
saveas(gcf, 'StockmanSharpe_ComparisonCIEConeFundamentalsTest.png', 'png');


%% Fit the tabulated absorbances as nomogram mixtures

% Some initialization information
params0.whichNomogram = 'Baylor';
theLambdaMaxes = [558.9 530.3 420.7];

% Set fmincon options
options = optimset('fmincon');
options = optimset(options,'Diagnostics','off','Display','off','LargeScale','off','Algorithm','active-set');

% For each pigment type (L, M, S)
for ii = 1:3
    % Initialize parameters
    params0.lambdaMax1 = theLambdaMaxes(ii);
    params0.lambdaMax2 = theLambdaMaxes(ii)-5;
    params0.weight1 = 0.5;
    
    % Extract absorbance to fit
    theAbsorbance = T_StockmanSharpeAbsorbance(ii,:);
    
    % Initial guess and bounds
    x0 = ParamsToX(params0);
    vlb = [x0(1)-20 x0(1)-20 0];
    vub = [x0(1)+20 x0(1)+20 1];
    
    % Search
    x1 = fmincon(@(x)FitNomogramErrorFunction(x,S,theAbsorbance,params0),x0,[],[],[],[],vlb,vub,[],options)
    params(ii) = XToParams(x1,params0);
    theAbsorbancePred(ii,:) = ComputeNomogramPred(params(ii),S);
end

figFitVsNomogram = figure;
hold on;
for ii = 1:3
    plot(wls, T_StockmanSharpeAbsorbance(ii, :), 'Color', theLMSCols(ii, :), 'LineWidth', 3);
    plot(wls, theAbsorbancePred(ii, :), 'Color', 'k', 'LineWidth', 1.5);
end
pbaspect([1 1 1]);
set(gca, 'TickDir', 'out'); box off;
xlim([375 785]); ylim([-0.01 1.01]);
xlabel('Wavelength [nm]'); ylabel('Photopigment absorbance');
title({'Tabulated vs. fit' 'Linear'});

end

function f = FitNomogramErrorFunction(x,S,theAbsorbance,params0)

params = XToParams(x,params0);
theAbsorbancePred = ComputeNomogramPred(params,S);
theDiff2 = (theAbsorbance-theAbsorbancePred).^2;
f = sqrt(mean(theDiff2));

end

function absorbancePred = ComputeNomogramPred(params,S)

absorbancePred = params.weight1*PhotopigmentNomogram(S,params.lambdaMax1,params.whichNomogram) + ...
    (1-params.weight1)*PhotopigmentNomogram(S,params.lambdaMax2,params.whichNomogram);

end

function x = ParamsToX(params)
x(1) = params.lambdaMax1;
x(2) = params.lambdaMax2;
x(3) = params.weight1;
end

function params = XToParams(x,params0)
params.whichNomogram = params0.whichNomogram;
params.lambdaMax1 = x(1);
params.lambdaMax2 = x(2);
params.weight1 = x(3);
end
