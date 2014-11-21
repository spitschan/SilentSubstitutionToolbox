function SaveSplatter(saveDir, fileNameSuffix, contrastMap, photoreceptorClasses, nominalLambdaMax, targetAge, ageRange, lambdaMaxShiftRange, targetContrasts)
% SaveSplatter(saveDir, fileNameSuffix, contrastMap, photoreceptorClasses, nominalLambdaMax, ageRange, targetAge, lambdaMaxShiftRange, targetContrasts)
%
% Saves out the contrast splatter map as well as some useful summary
% statistics.
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
%
% 11/21/14   ms    Cleaned up and commented

% If previously undefined, save it out in the current directory
if isempty(saveDir)
    saveDir = pwd;
end

%% Save out the contrast map as .mat
saveFile = fullfile(saveDir, ['Splatter' fileNameSuffix '.mat']);
save(saveFile, 'contrastMap', 'photoreceptorClasses', 'nominalLambdaMax', 'ageRange', 'lambdaMaxShiftRange', 'targetContrasts');
fprintf('  - Contrast map saved to %s.\n', fullfile(saveDir, ['Splatter' fileNameSuffix '.mat']));

%% Compile a few statistics on the splatter and save out as csv
fileID = fopen(fullfile(saveDir, ['Splatter_statistics' fileNameSuffix '.csv']), 'w');
fprintf(fileID,'Class,NominalLambdaMax,LambdaMaxShiftMin,LambdaMaxShiftMax,ageRangeMin,ageRangeMax,targetContrast,measuredTargetContrast,meanContrast,stdContrast,meanAbsContrast,stdAbsContrast,minContrast,maxContrast\n');
for k = 1:length(photoreceptorClasses)
    fprintf(fileID,'%s,%.2f,%.2f,%.2f,%.2f,%.2f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n',photoreceptorClasses{k}, nominalLambdaMax(k), min(lambdaMaxShiftRange), max(lambdaMaxShiftRange), min(ageRange), max(ageRange), targetContrasts{k}, contrastMap{k}(find(nominalLambdaMax(k)-lambdaMaxShiftRange == nominalLambdaMax(k)), find(ageRange == targetAge)), mean(contrastMap{k}(:)), std(contrastMap{k}(:)), mean(abs(contrastMap{k}(:))), std(abs(contrastMap{k}(:))), min(contrastMap{k}(:)), max(contrastMap{k}(:)));
end
fclose(fileID);
fprintf('  - Contrast statistics saved to %s.\n', fullfile(saveDir, ['Splatter_statistics' fileNameSuffix '.csv']));