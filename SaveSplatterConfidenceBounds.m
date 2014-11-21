function SaveSplatterConfidenceBounds(saveDir, fileNameSuffix, contrastMap, photoreceptorClasses, nominalLambdaMax, ageRange, lambdaMaxShiftRange, targetContrasts, confidenceInterval)
% SaveSplatterConfidenceBounds(saveDir, fileNameSuffix, contrastMap, photoreceptorClasses, nominalLambdaMax, ageRange, lambdaMaxShiftRange, targetContrasts, confidenceInterval)
%
% Saves out the contrast splatter map as well as the statistics on the
% splatter within an error ellipse.
%
% Input:
%   saveDir (str)               - Directory to save map and statistics.
%   fileNameSuffix (str)        - Suffix for file name
%   photoreceptorClasses (cell) - Photoreceptor classes
%   nominalLambdaMax (1xK)      - Nominal lambda-max
%   targetAge                   - Target age
%   ageRange                    - Age range of contrast map
%   lambdaMaxShiftRange         - Shift range of lambda-max
%   targetContrast (1xK)        - Target contrasts
%   confidenceInterval (1x1)    - Confidence interval o finterest

% 11/21/14   ms    Cleaned up and commented

% If previously undefined, save it out in the current directory
if isempty(saveDir)
    saveDir = pwd;
end

%% Compile a few statistics on the splatter and save out as csv
fileID = fopen(fullfile(saveDir, ['Splatter_statistics' fileNameSuffix '.csv']), 'w');
fprintf(fileID,'Class,NominalLambdaMax,LambdaMaxShiftMin,LambdaMaxShiftMax,targetAge,ageRangeMin,ageRangeMax,targetContrast,measuredTargetContrast,meanContrast,stdContrast,meanAbsContrast,stdAbsContrast,minContrast,maxContrast,expectedValue\n');

for o = ageRange
    targetAge = o;
    ageMean = targetAge;
    for k = 1:length(photoreceptorClasses)
        %% Lambda-max
        lambdaMaxSD = GetLambdaMaxEstimateSD(photoreceptorClasses{k}, 'M&W');
        
        %% Age SD from lens
        observerAge = 0;
        ageSD = GetChronicalAgeSDFromLensSD('Xu');
        
        %% Construct the 2D Gaussian
        mu = [nominalLambdaMax(k) ageMean];
        Sigma = [lambdaMaxSD ageSD].^2;
        x1 = nominalLambdaMax(k)+lambdaMaxShiftRange;
        x2 = ageRange;
        [X1,X2] = meshgrid(x1,x2);
        F = mvnpdf([X1(:) X2(:)],mu,Sigma);
        
        % Get the stepsize to normalize the pdf to get pmf
        stepSize = median(diff(x1))*median(diff(x2));
        
        F = reshape(F,length(x2),length(x1))*stepSize;
        
        % Normalize volume to 1 just in case that's not the case.
        F = F/sum(sum(F));
        
        colormap(gray); imagesc(x1,x2,F); hold on;
        [~, x, y] = error_ellipse([lambdaMaxSD.^2 0 ; 0 ageSD.^2], mu, 'conf', confidenceInterval);
        
        plot(x, y); hold on;
        
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
        
        % Save out the confidence bounds
        fprintf(fileID,'%s,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n',photoreceptorClasses{k}, nominalLambdaMax(k), min(lambdaMaxShiftRange), max(lambdaMaxShiftRange), targetAge, min(ageRange), max(ageRange), targetContrasts{k}, contrastMap{k}(find(nominalLambdaMax(k)-lambdaMaxShiftRange == nominalLambdaMax(k)), find(ageRange == targetAge)), mean(contrastMap{k}(:)), std(contrastMap{k}(indexMap)), mean(abs(contrastMap{k}(indexMap))), std(abs(contrastMap{k}(indexMap))), min(contrastMap{k}(indexMap)), max(contrastMap{k}(indexMap)), sum(sum(F .* abs(contrastMap{k})')));
    end
end
fclose(fileID);
fprintf('  - Contrast statistics saved to %s.\n', fullfile(saveDir, ['Splatter_statistics' fileNameSuffix '.csv']));