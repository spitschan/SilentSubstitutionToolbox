function theFig = PlotSplatter(theFig, contrastMap, photoreceptorClasses, nominalLambdaMax, targetAge, ageRange, lambdaMaxShiftRange, targetContrast, maxContrast, plotRow, maxRows, SAVEPLOTS, titleSuffix, plotFileName, ageMean)
% theFig = PlotSplatter(theFig, contrastMap, photoreceptorClasses, nominalLambdaMax, targetAge, ageRange, lambdaMaxShiftRange, targetContrast, maxContrast, plotRow, maxRows, SAVEPLOTS, titleSuffix, plotFileName, ageMean)
%
% Produces contrast splatter map.
%
% Input:
%   theFig                          - Figure handle to plot into.
%   contrastMap (cell)              - Cell with map of contrasts.
%   photoreceptorClasses (cell)     - Names of receptors. Should be of same
%                                     dimension as the contrast map cell.
%   nominalLambdaMax (1xK)          - Nominal lambda-max. Should be of same
%                                     dimension as the contrast map cell.
%   targetAge                       - Target age
%   ageRange (1xM)                  - Vector of observer age
%   lambdaMaxShiftRange (1xN)       - Vector of lambda-max shifts
%   targetContrast (cell, 1xK)      - Cell array of the target contrasts
%   maxContrast (cell, 1x2)         - Min. and max. expected contrast.
%                                     Defaults to: [-0.1 0.1] for all
%                                     entries in cell array
%   plotRow                         - Subplot row to plot into.
%   maxRows                         - Maximum number of rows.
%   SAVEPLOTS (bool)                - Flag to save plots.
%   plotFileName (str)              - File name to save out plot.
%   ageMean (1x1)                   - Average age of observers
%
% Output:
%   theFig (1x1)                    - Figure handle
%
% In principle there are two modes of operation, one passes in the target
% contrasts such that one can look at the deviation from the target
% contrasts, or one passes the maximum contrasts to be plotted, then the
% scales are shifted. The behavior changes a little bit, including labeling
% of the color bar (contrast becomes delta-contrast, etc.).
%
%
% 1/21/14   ms      Wrote it.
% 5/26/14   dhb     Fix comment: targetContrast is a cell array.
% 11/21/14  ms      Cleaned up and commented.

% Do a few checks and assign defualt values if not defined otherwise.
if length(contrastMap) > length(photoreceptorClasses)
    error('Contrast map cell array is not same as size as receptor labels passed');
end

if isempty(maxContrast)
    for k = 1:length(contrastMap)
        maxContrast{k} = [-0.1 0.1];
    end
end

if ~isempty(targetContrast)
    targetContrastSupplied = true;
    if length(contrastMap) ~= length(targetContrast)
        error('Contrast map cell array is not same as size as targetcontrast vector passed');
    end
else
    targetContrastSupplied = false;
    for k = 1:length(contrastMap)
        targetContrast{k} = 0;
    end
end

for k = 1:length(contrastMap)
    if size(contrastMap{k}, 1) ~= length(lambdaMaxShiftRange) || size(contrastMap{k}, 2) ~= length(ageRange)
        error('Contrast map for at least one receptor does not correspond to parametric vectors passed');
    end
end

if SAVEPLOTS && isempty(plotFileName)
    plotFileName = 'tmp';
end

%% Open the figure
figure(theFig);

% Set up the color map. In case you are wondering why we are using flipud,
% we want red to be positive contrast and blue to be negative
colormap(flipud(lbmap(400, 'BrownBlue')));

%% Iterate over the k photoreceptor-specific maps in the contrast map array
for k = 1:length(contrastMap)
    % Figure out where to plot this
    plotLocation = k + (length(contrastMap)*(plotRow-1));
    s1 = subplot(maxRows, length(contrastMap), plotLocation);
    imagesc(ageRange, nominalLambdaMax(k)+lambdaMaxShiftRange, contrastMap{k}-targetContrast{k}, maxContrast{k});
    hold on;
    xlabel('Age [yrs]');
    ylabel('\lambda_{max} [nm]')
    
    % Find and plot extreme points
    [r, c] = find(contrastMap{k}==max(contrastMap{k}(:)));
    plot(ageRange(c), nominalLambdaMax(k)+lambdaMaxShiftRange(r), 'ok', 'MarkerFaceColor', 'k');
    
    [r, c] = find(contrastMap{k}==min(contrastMap{k}(:)));
    plot(ageRange(c), nominalLambdaMax(k)+lambdaMaxShiftRange(r), 'ok', 'MarkerFaceColor', 'w');
    
    % Add a cross hair
    plot(repmat(targetAge, 1, length(lambdaMaxShiftRange)), nominalLambdaMax(k)+lambdaMaxShiftRange, '--k');
    plot(ageRange, repmat(nominalLambdaMax(k), 1, length(ageRange)), '--k');
    
    %% Get the confidence ellipses: 95%
    lambdaMaxSD = GetLambdaMaxEstimateSD(photoreceptorClasses{k}, 'M&W');
    ageSD = GetChronicalAgeSDFromLensSD('Xu');
    
    % Construct the 2D Gaussian
    mu = [nominalLambdaMax(k) ageMean];
    Sigma = [lambdaMaxSD ageSD].^2;
    x1 = nominalLambdaMax(k)+lambdaMaxShiftRange;
    x2 = ageRange;
    [X1,X2] = meshgrid(x1,x2);
    F = mvnpdf([X1(:) X2(:)],mu,Sigma);
    
    % Get the stepsize to normalize the pdf to get pmf
    stepSize = median(diff(x1))*median(diff(x2));
    
    %F = reshape(F,length(x2),length(x1))*stepSize;
    %colormap(gray); imagesc(x1,x2,F); hold on;
    [~, x, y] = error_ellipse([lambdaMaxSD.^2 0 ; 0 ageSD.^2], mu, 'conf', 0.9545);
    
    plot(y, x, '-k'); hold on;
    
    % For any given x, y position in contrastMap, test if it's in the
    % error_ellipse. This gives us the 95% confidence interval.
    for i = 1:length(x1)
        for j = 1:length(x2)
            x0 = X1(j, i);
            y0 = X2(j, i);
            xr = x - x0; yr = y - y0; % Translate the vectors' base to (x0,y0)
            xw = xr([2:end,1]); yw = yr([2:end,1]); % Shift ahead by one index
            s = sum(atan2(xr.*yw-yr.*xw,xr.*xw+yr.*yw)); % Sum the projected angles
            in = abs(s) > pi; % 'in' is true if (x0,y0) is inside, otherwise false
            indexMap(i, j) = in;
        end
    end
    
    % Find and plot extreme points
    [r, c] = find(contrastMap{k}==max(contrastMap{k}(indexMap)));
    plot(ageRange(c), nominalLambdaMax(k)+lambdaMaxShiftRange(r), 'ok', 'MarkerFaceColor', 'k');
    
    [r, c] = find(contrastMap{k}==min(contrastMap{k}(indexMap)));
    plot(ageRange(c), nominalLambdaMax(k)+lambdaMaxShiftRange(r), 'ok', 'MarkerFaceColor', 'w');
    
    %% Get the confidence ellipses: 99%
    lambdaMaxSD = GetLambdaMaxEstimateSD(photoreceptorClasses{k}, 'M&W');
    ageSD = GetChronicalAgeSDFromLensSD('Xu');
    
    % Construct the 2D Gaussian
    mu = [nominalLambdaMax(k) ageMean];
    Sigma = [lambdaMaxSD ageSD].^2;
    x1 = nominalLambdaMax(k)+lambdaMaxShiftRange;
    x2 = ageRange;
    [X1,X2] = meshgrid(x1,x2);
    F = mvnpdf([X1(:) X2(:)],mu,Sigma);
    
    % Get the stepsize to normalize the pdf to get pmf
    stepSize = median(diff(x1))*median(diff(x2));
    
    [~, x, y] = error_ellipse([lambdaMaxSD.^2 0 ; 0 ageSD.^2], mu, 'conf', 0.9973);
    
    plot(y, x, '-k'); hold on;
    
    % For any given x, y position in contrastMap, test if it's in the
    % error_ellipse. This gives us the 95% confidence interval.
    for i = 1:length(x1)
        for j = 1:length(x2)
            x0 = X1(j, i);
            y0 = X2(j, i);
            xr = x - x0; yr = y - y0; % Translate the vectors' base to (x0,y0)
            xw = xr([2:end,1]); yw = yr([2:end,1]); % Shift ahead by one index
            s = sum(atan2(xr.*yw-yr.*xw,xr.*xw+yr.*yw)); % Sum the projected angles
            in = abs(s) > pi; % 'in' is true if (x0,y0) is inside, otherwise false
            indexMap(i, j) = in;
        end
    end
    
    % Find and plot extreme points
    [r, c] = find(contrastMap{k}==max(contrastMap{k}(indexMap)));
    plot(ageRange(c), nominalLambdaMax(k)+lambdaMaxShiftRange(r), 'ok', 'MarkerFaceColor', 'k');
    
    [r, c] = find(contrastMap{k}==min(contrastMap{k}(indexMap)));
    plot(ageRange(c), nominalLambdaMax(k)+lambdaMaxShiftRange(r), 'ok', 'MarkerFaceColor', 'w');
    
    % Add a color bar
    s1Pos = get(s1,'position');
    h2 = colorbar('Location','SouthOutside');
    if targetContrastSupplied
        xlabel(h2, {'{\Delta}Contrast' ; ['Deviation from ' num2str(targetContrast{k}, '%.3f')]});
    else
        xlabel(h2, 'Contrast');
    end
    xlim([18 80]);
    set(s1,'position',s1Pos);
    pbaspect([1 1 1]);
    
    % Add a title
    title({[titleSuffix ' modulation: ' photoreceptorClasses{k}] ; ['max: ' num2str(max(contrastMap{k}(:))) ' / min: ' num2str(min(contrastMap{k}(:)))] ; ['at target: ' num2str(contrastMap{k}(find(lambdaMaxShiftRange == 0), find(ageRange == targetAge)))]});
end

if SAVEPLOTS
    %% Save plots
    set(theFig, 'Color', [1 1 1]);
    set(theFig, 'InvertHardCopy', 'off');
    set(theFig, 'PaperPosition', [0 0 20 10]); %Position plot at left hand corner with width 15 and height 6.
    set(theFig, 'PaperSize', [20 10]); %Set the paper to have width 15 and height 6.
    saveas(theFig, plotFileName, 'pdf');
end