%% HumanEyelidTransmittanceBierman.m
%
% This program calculates the transmittance of the human eyelid. It
% implements the method described in:
%
% Bierman, A., Figueiro, M. G., & Rea, M. S. (2011). Measuring and
% predicting eyelid spectral transmittance. J Biomed Opt, 16(6), 067011.
% doi:10.1117/1.3593151
%
% 11/29/16  spitschan      Wrote it.

%% Dependencies
% `SplineSrf` from PsychToolbox.

%% Source data
% The folder `xSrc` contains files with tabulated data necessary for this
% demo program:
%
% ? Data from Bierman, Figueiro & Rea (2011):
%
%   1) Bierman2011_EyelidTransFromAvgCoef.txt - Predicted eyelid
%   transmittance calculated from the average coefficients. Corresponds to
%   the dashed line in Fig. 5. Sent to me (MS) by Bierman, Figueiro & Rea
%   via personal communication.
%
%   2) Bierman2011_meanTransmittanceN27.txt - Mean measured transmittance
%   calculated from the average coefficients. Corresponds to the solid line
%   in Fig. 5. Sent to me (MS) by Bierman, Figueiro & Rea via personal
%   communication.
%
%   3) Bierman2011_Table_Corrected.txt - Table with the coefficients for
%   the 3-component model (deoxy-hemoglobin, oxy-hemoglobin, melanin,
%   scatter). In the table printed in Bierman, Figueiro & Rea (2011), the
%   scatter coefficients for the subjects 1-17 were wrong and are corrected
%   here. Sent to me (MS) by Bierman, Figueiro & Rea via personal
%   communication.
%
%   4) Bierman2011_Table.txt - Table with the coefficients for
%   the 3-component model (deoxy-hemoglobin, oxy-hemoglobin, melanin,
%   scatter). This is the table printed in the paper with the wrong scatter
%   coefficients described in 3).
%
% - Data from Prahl (2012):
%
%   1) Prahl2012_BilurubinAbsorptionSpectrum.txt - Optical absorption
%   spectrum of bilurubin used in Bierman, Figueiro & Rea (2011). This was
%   obtained from http://omlc.org/spectra/PhotochemCAD/html/119.html (Nov
%   28, 2016).
%
% - Data from Moseley, Bayliss & Fielder (1988):
%
%   1) Moseley1988_EyelidTransmittance.txt ? Eyelid transmittance spectra
%   for 3 adults (Fig. 1A in their paper). This was digitized using
%   WebPlotDigitizer (Nov 28, 2016). Full reference: Moseley, M. J.,
%   Bayliss, S. C., & Fielder, A. R. (1988). Light transmission through the
%   human eyelid: in vivo measurement. Ophthalmic Physiol Opt, 8(2),
%   229-230.

%% Wavelength spacing
% Define the wavelength spacing
S = [380 2 201];
wls = SToWls(S);

%% Melanin absorption
% The skin melanin absorption coefficient can be described by the following
% formula (eq. 2 in Bierman, Figuiero & Rea, 2011; and
% http://omlc.org/spectra/melanin/mua.html):
absorption_melanin = 1.70*10^12*wls.^(-3.48);

% This needs to be scaled by 0.002059. This is not mentioned in the paper
% (personal communication with Bierman, Figueiro & Rea).
absorption_melanin = 0.002059*absorption_melanin;

% Turn into transmittance
trans_melanin = exp(-absorption_melanin);

%% Bilirubin 
% Load bilirubin
tmp = csvread(fullfile('xSrc', 'Prahl2012_BilurubinAbsorptionSpectrum.txt'), 25, 0);
wls_bilirubin = tmp(:, 1);
absorption_bilirubin = tmp(:, 2);

% This needs to be scaled by 2.0*10^-5. This is not mentioned in the paper
% (personal communication with Bierman, Figueiro & Rea).
absorption_bilirubin = 2.0*10^-5*absorption_bilirubin;

% Spline to get the correct wavelength spacing
absorption_bilirubin = SplineSrf(wls_bilirubin, absorption_bilirubin, S);

% Turn into transmittance
trans_bilirubin = exp(-absorption_bilirubin);

%% Hemoglobin 
% We already have the hemoglobin absorption data in the function
% `GetHemoglobin`, so we just use it here.
[S_Prahl, wls_Prahl, absorption_oxyhemoglobin, absorption_deoxyhemoglobin] = GetHemoglobin('Prahl');

% This needs to be scaled by 2.0*10^-5. This is not mentioned in the paper
% (personal communication with Bierman, Figueiro & Rea).
absorption_oxyhemoglobin = 2.0*10^-5*(absorption_oxyhemoglobin);
absorption_deoxyhemoglobin = 2.0*10^-5*(absorption_deoxyhemoglobin);

% Spline to get the correct wavelength spacing
absorption_oxyhemoglobin = SplineSrf(S_Prahl, absorption_oxyhemoglobin, S);
absorption_deoxyhemoglobin = SplineSrf(S_Prahl, absorption_deoxyhemoglobin, S);

% Turn into transmittance
trans_oxyhemoglobin = exp(-absorption_oxyhemoglobin);
trans_deoxyhemoglobin = exp(-absorption_deoxyhemoglobin);

%% Plot the transmittance spectra as in Fig. 4 of Bierman, Figueiro & Rea (2011)
gc1Fig = figure;
h1 = plot(wls, trans_oxyhemoglobin, 'LineWidth', 2); hold on;
h2 = plot(wls, trans_deoxyhemoglobin, 'LineWidth', 2);
h3 = plot(wls, trans_melanin, 'LineWidth', 2);
h4 = plot(wls, trans_bilirubin, 'LineWidth', 2)
ylim([0.1 1.01]);
xlim([450 700]);
xlabel('Wavelength [nm]');
ylabel('Transmittance');
set(gca, 'TickDir', 'out'); box off; pbaspect([1 1 1]);
legend([h1 h2 h3 h4], 'Oxy-hemoglobin', 'Deoxy-hemoglobin', 'Melanin', 'Bilirubin', 'Location', 'SouthEast'); legend boxoff;
title({'Human eyelid transmittance' 'Components'});

% Save the figure
set(gc1Fig, 'PaperPosition', [0 0 4 4]);
set(gc1Fig, 'PaperSize', [4 4]);
saveas(gc1Fig, 'HumanEyelidTransmittanceBierman_Components.pdf', 'pdf');

%% Load in the data from Rea et al.
% Average transmittance from 27 subjects
tmp = csvread(fullfile('xSrc', 'Bierman2011_meanTransmittanceN27.txt'));
wls_eyelidMeas = tmp(:, 1);
trans_eyelidMeas = tmp(:, 2);

% Transmittance from average coeffiecients
tmp = csvread(fullfile('xSrc', 'Bierman2011_EyelidTransFromAvgCoef.txt'));
wls_eyelidPred = tmp(:, 1);
trans_eyelidPredOrig = tmp(:, 2);

% Construct the transmittance from the mean coefficients from the table in Bierman, Figueiro & Rea (2011)
w_deoxyhemoglobin = 0.739;
w_oxyhemoglobin = 1.368;
w_melanin = 2.643;
log_trans_scatter = -0.825;

% Put together
absorption_eyelid = w_deoxyhemoglobin*log(trans_deoxyhemoglobin) + ...
    w_oxyhemoglobin*log(trans_oxyhemoglobin) + ...
    w_melanin*log(trans_melanin) + ...
    log_trans_scatter + 0.1825505; % The constant 0.1825505 ~= log(1.2) corresponds to the best fitting scalar

% Turn into transmittance
trans_eyelidPred = exp(absorption_eyelid);

%% Plot these data
% We're plotting a) our prediction, b) the predicted tabulated
% transmittance (personal communication from Bierman, Figueiro & Rea, 2011)
% - as a check, and c) the average measured transmittance (personal
% communication from Bierman, Figueiro & Rea, 2011).
gc2Fig = figure;
h1 = semilogy(wls, trans_eyelidPred, '-r', 'LineWidth', 2); hold on;
h2 = semilogy(wls_eyelidPred, trans_eyelidPredOrig, '--k', 'LineWidth', 2);
h3 = semilogy(wls_eyelidMeas, trans_eyelidMeas, '-b', 'LineWidth', 2);
xlim([450 700]);
xlabel('Wavelength [nm]');
ylabel('Transmittance');
set(gca, 'TickDir', 'out'); box off; pbaspect([1 1 1]);
legend([h1 h2 h3], 'Predicted transmittance (analytic)', 'Predicted transmittance (tabulated)', ...
    'Measured transmittance', 'Location', 'SouthEast'); legend boxoff;
title({'Human eyelid transmittance'});

% Save the figure
set(gc2Fig, 'PaperPosition', [0 0 4 4]);
set(gc2Fig, 'PaperSize', [4 4]);
saveas(gc2Fig, 'HumanEyelidTransmittanceBierman_Transmittance.pdf', 'pdf');

%% Compare with Moseley et al.'s data
% Load Moseley et al. 1988
gc3Fig = figure;
tmp = csvread(fullfile('xSrc', 'Moseley1988_EyelidTransmittance.txt'), 1, 0);
wls_Moseley = tmp(:, 1);
trans_Moseley = tmp(:, 2)/100; % Divide by 100 to get prop. and not percent
gc3Fig = figure;
h1 = plot(wls, trans_eyelidPred, '-r', 'LineWidth', 2); hold on;
h2 = plot(wls_Moseley, trans_Moseley, '-sk', 'MarkerFaceColor', 'k');
xlim([400 720]); ylim([0 0.2]);
set(gca, 'TickDir', 'out'); box off; pbaspect([1 1 1]);
xlabel('Wavelength [nm]');
ylabel('Transmittance');
legend([h1 h2], 'Bierman, Figueiro & Rea (2011) model', 'Moseley, Bayliss & Fielder (1988) data', 'Location', 'NorthWest'); legend boxoff;
title({'Human eyelid transmittance'});

% Save the figure
set(gc3Fig, 'PaperPosition', [0 0 4 4]);
set(gc3Fig, 'PaperSize', [4 4]);
saveas(gc3Fig, 'HumanEyelidTransmittanceBierman_ComparisonMoseley.pdf', 'pdf');