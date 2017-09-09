%% Do a little bit of house keeping
clear all; close all; clc

%% Set the wavelength spacing
S = [380 2 201];
wls = SToWls(S);

%% Construct the receptor object
receptorObj = SSTReceptorHuman('S', S, 'verbosity', 'high', 'obsAgeInYrs', 32);

% Iterate over the LED bins
binStarts = [441 471 501 531 561 591 621 651];
binEnds = [470 500 530 550 590 620 650 680];
NSteps = 30;
NBins = 8;
for ii = 1:NBins
    eval(['bin' num2str(ii) 'vals = binStarts(' num2str(ii) '):binEnds(' num2str(ii) ');']);
end
% Construct LEDs
%%
currentBestContrast = 0;

c = 1;
for i8 = 1:NSteps
    for i7 = 1:NSteps
        for i6 = 1:NSteps
            for i5 = 1:NSteps
                for i4 = 1:NSteps
                    for i3 = 1:NSteps
                        for i2 = 1:NSteps
                            for i1 = 1:NSteps
                                NLEDs = 8;
                                peakWls = [bin1vals(i1) bin2vals(i2) bin3vals(i3) bin4vals(i4) bin5vals(i5) bin6vals(i6) bin7vals(i7) bin8vals(i8)];
                                fwhm = 10*ones(1, NLEDs);
                                maxPower = ones(1, NLEDs);
                                for i = 1:length(fwhm)
                                    % Figure out the standard deviation.
                                    standardDeviation(i) = FWHMToStd(fwhm(i));
                                end
                                
                                % Make the spectrum.
                                for i = 1:length(fwhm)
                                    spd(:, i) = normpdf(wls, peakWls(i), fwhm(i));
                                    spd(:, i) = spd(:, i)./max(spd(:, i))*maxPower(i);
                                end
                                
                                % Assign the basis matrix
                                B_primary = spd;
                                
                                % Assign the ambient
                                ambientSpd = zeros(size(wls));
                                
                                % Set background to the device midpoint.
                                backgroundPrimary = 0.5*ones(size(B_primary, 2), 1);
                                
                                % Don't pin any primaries.  Do enforce a constraint that we don't
                                % go right to the edge of the gamut.  The head room parameter is
                                % defined in the [0-1] device primary space.  Using a little head
                                % room keeps us a bit away from the hard edge of the device.
                                whichPrimariesToPin = [];
                                primaryHeadRoom = 0.02;
                                
                                % No smoothness constraint envforced here.  It really wouldn't make
                                % to much sense for a three-primary monitor, as the smoothness of a
                                % monitor spectrum is pretty much determined by the spectral shape
                                % of its primarites.
                                maxPowerDiff = 0.005;
                                
                                whichDirection = 'MelDirected';
                                %%
                                whichReceptorsToTarget = {[3] [4]};
                                whichReceptorsToIgnore = {[] []};
                                whichReceptorsToMinimize = {[] []};
                                desiredContrasts = {[] []};
                                directionsYoked = [0 0];
                                directionsYokedAbs = [0 0];
                                
                                % Isolate the receptors by calling the wrapper
                                initialPrimary = backgroundPrimary;
                                unipolarYesNo = false;
                                pegBackground = false;
                                [modulationPrimary backgroundPrimary] = ReceptorIsolateOptim(receptorObj.T.T_energy, whichReceptorsToTarget, ...
                                    whichReceptorsToIgnore,whichReceptorsToMinimize,B_primary,backgroundPrimary,...
                                    initialPrimary,whichPrimariesToPin,primaryHeadRoom,maxPowerDiff,...
                                    desiredContrasts,ambientSpd,directionsYoked,directionsYokedAbs,pegBackground,unipolarYesNo);
                                
                                % Compute the contrasts that we got.
                                for ii = 1:size(modulationPrimary, 2)
                                    bgSpd = B_primary*backgroundPrimary + ambientSpd;
                                    modSpd = B_primary*modulationPrimary{ii} + ambientSpd;
                                    backgroundReceptors = receptorObj.T.T_energy*bgSpd;
                                    modulationReceptors = receptorObj.T.T_energy*modSpd;
                                    contrastReceptors = (modulationReceptors-backgroundReceptors) ./ backgroundReceptors;
                                    %fprintf('\n');
                                    for j = 1:size(receptorObj.T.T_energy,1)
                                        %fprintf('\t%s: contrast = %0.4f\n',receptorObj.labels{j},contrastReceptors(j));
                                    end
                                    %plot(bgSpd, '-k'); hold on;
                                    %pplot(modSpd, '-r');
                                end
                                peakWls0(:, c) = peakWls';
                                contrastReceptors0(:, c) = contrastReceptors;
                                if contrastReceptors(4) > currentBestContrast
                                    currentBestContrast = contrastReceptors(4)
                                end
                                c = c+1;
                            end
                        end
                    end
                end
            end
        end
    end
end