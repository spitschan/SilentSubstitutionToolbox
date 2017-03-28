% Define the paths
dropboxBasePath = '/Users/spitschan/Dropbox (Aguirre-Brainard Lab)';

theDataPaths = {'MELA_data/MelanopsinMR_fMRI/MaxLMS400pct/HERO_asb1/040716/StimulusFiles/Cache-LMSDirectedSuperMaxLMS/BoxARandomizedLongCableCStubby1_ND00/23-Mar-2016_12_31_27/validation' ...
    'MELA_data/MelanopsinMR_fMRI/MaxLMS400pct/HERO_aso1/033016/StimulusFiles/Cache-LMSDirectedSuperMaxLMS/BoxARandomizedLongCableCStubby1_ND00/23-Mar-2016_12_31_27/validation' ...
    'MELA_data/MelanopsinMR_fMRI/MaxLMS400pct/HERO_gka1/040116/StimulusFiles/Cache-LMSDirectedSuperMaxLMS/BoxARandomizedLongCableCStubby1_ND00/23-Mar-2016_12_31_27/validation' ...
    'MELA_data/MelanopsinMR_fMRI/MaxLMS400pct/HERO_mxs1/040816/StimulusFiles/Cache-LMSDirectedSuperMaxLMS/BoxARandomizedLongCableCStubby1_ND00/23-Mar-2016_12_31_27/validation' ...
    'MELA_data/MelanopsinMR_fMRI/MaxLMSCRF/HERO_asb1/060816/StimulusFiles/Cache-LMSDirectedSuperMaxLMS/BoxARandomizedLongCableBStubby1_ND00/29-May-2016_14_48_42/validation' ...
    'MELA_data/MelanopsinMR_fMRI/MaxLMSCRF/HERO_aso1/060116/StimulusFiles/Cache-LMSDirectedSuperMaxLMS/BoxARandomizedLongCableBStubby1_ND00/29-May-2016_14_48_42/validation' ...
    'MELA_data/MelanopsinMR_fMRI/MaxLMSCRF/HERO_gka1/060616/StimulusFiles/Cache-LMSDirectedSuperMaxLMS/BoxARandomizedLongCableBStubby1_ND00/29-May-2016_14_48_42/validation' ...
    'MELA_data/MelanopsinMR_fMRI/MaxLMSCRF/HERO_mxs1/061016/StimulusFiles/Cache-LMSDirectedSuperMaxLMS/BoxARandomizedLongCableBStubby1_ND00/29-May-2016_14_48_42/validation' ...
    'MELA_data/MelanopsinMR_fMRI/MaxLMSCRF/HERO_mxs1/062816/StimulusFiles/Cache-LMSDirectedSuperMaxLMS/BoxARandomizedLongCableBStubby1_ND03CassetteB/20-Jun-2016_15_17_02/validation' ...
    'MELA_data/MelanopsinMR_fMRI/MaxMel400pct/HERO_asb1/032416/StimulusFiles/Cache-MelanopsinDirectedSuperMaxMel/BoxARandomizedLongCableCStubby1_ND00/23-Mar-2016_12_31_27/validation' ...
    'MELA_data/MelanopsinMR_fMRI/MaxMel400pct/HERO_aso1/032516/StimulusFiles/Cache-MelanopsinDirectedSuperMaxMel/BoxARandomizedLongCableCStubby1_ND00/23-Mar-2016_12_31_27/validation' ...
    'MELA_data/MelanopsinMR_fMRI/MaxMel400pct/HERO_gka1/033116/StimulusFiles/Cache-MelanopsinDirectedSuperMaxMel/BoxARandomizedLongCableCStubby1_ND00/23-Mar-2016_12_31_27/validation' ...
    'MELA_data/MelanopsinMR_fMRI/MaxMel400pct/HERO_mxs1/040616/StimulusFiles/Cache-MelanopsinDirectedSuperMaxMel/BoxARandomizedLongCableCStubby1_ND00/23-Mar-2016_12_31_27/validation' ...
    'MELA_data/MelanopsinMR_fMRI/MaxMelCRF/HERO_asb1/060716/StimulusFiles/Cache-MelanopsinDirectedSuperMaxMel/BoxARandomizedLongCableBStubby1_ND00/29-May-2016_14_48_42/validation' ...
    'MELA_data/MelanopsinMR_fMRI/MaxMelCRF/HERO_aso1/053116/StimulusFiles/Cache-MelanopsinDirectedSuperMaxMel/BoxARandomizedLongCableBStubby1_ND00/29-May-2016_14_48_42/validation' ...
    'MELA_data/MelanopsinMR_fMRI/MaxMelCRF/HERO_gka1/060216/StimulusFiles/Cache-MelanopsinDirectedSuperMaxMel/BoxARandomizedLongCableBStubby1_ND00/29-May-2016_14_48_42/validation' ...
    'MELA_data/MelanopsinMR_fMRI/MaxMelCRF/HERO_mxs1/060916/StimulusFiles/Cache-MelanopsinDirectedSuperMaxMel/BoxARandomizedLongCableBStubby1_ND00/29-May-2016_14_48_42/validation' ...
    'MELA_data/MelanopsinMR_fMRI/MaxMelCRF/HERO_mxs1/061016/StimulusFiles/Cache-MelanopsinDirectedSuperMaxMel/BoxARandomizedLongCableBStubby1_ND00/29-May-2016_14_48_42/validation' ...
    'MELA_data/MelanopsinMR_fMRI/SplatterControlCRF/HERO_asb1/051016/StimulusFiles/Cache-MaxMelPostreceptoralSplatterControl/BoxARandomizedLongCableCStubby1_ND00/18-Apr-2016_10_23_32/validation' ...
    'MELA_data/MelanopsinMR_fMRI/SplatterControlCRF/HERO_aso1/042916/StimulusFiles/Cache-MaxMelPostreceptoralSplatterControl/BoxARandomizedLongCableCStubby1_ND00/18-Apr-2016_10_23_32/validation' ...
    'MELA_data/MelanopsinMR_fMRI/SplatterControlCRF/HERO_gka1/050616/StimulusFiles/Cache-MaxMelPostreceptoralSplatterControl/BoxARandomizedLongCableCStubby1_ND00/18-Apr-2016_10_23_32/validation' ...
    'MELA_data/MelanopsinMR_fMRI/SplatterControlCRF/HERO_mxs1/050916/StimulusFiles/Cache-MaxMelPostreceptoralSplatterControl/BoxARandomizedLongCableCStubby1_ND00/18-Apr-2016_10_23_32/validation'}

theStimuli = {'LMS 400%' 'LMS 400%'  'LMS 400%'  'LMS 400%'  ...
    'LMS CRF'  'LMS CRF'  'LMS CRF'  'LMS CRF' 'LMS CRF'...
    'Mel 400%' 'Mel 400%'  'Mel 400%'  'Mel 400%'  ...
    'Mel CRF'  'Mel CRF'  'Mel CRF'  'Mel CRF' 'Mel CRF' ...
    'Splatter CRF' 'Splatter CRF' 'Splatter CRF'  'Splatter CRF'};
theObservers = {'ASB' 'ASO' 'GKA' 'MXS' ...
    'ASB' 'ASO' 'GKA' 'MXS [1]' 'MXS [2]'...
    'ASB' 'ASO' 'GKA' 'MXS' ...
    'ASB' 'ASO' 'GKA' 'MXS [1]' 'MXS [2]' ...
    'ASB' 'ASO' 'GKA' 'MXS'}

fid = fopen('~/Desktop/test.csv', 'w');
fprintf(fid, 'Stimulus,Observer,Luminance,SD,Chromaticity x,SD,Chromaticity y,SD,L contrast [%s],SD,M contrast [%s],SD,S contrast [%s],SD,LMS contrast [%s],SD,L-M contrast [%s],SD\n', '%', '%', '%', '%', '%');


currDir = pwd;
Mc = [];
for d = 1:length(theDataPaths)
    dataPath = theDataPaths{d};
    
    % Find the folders
    theFolders = dir(fullfile(dropboxBasePath, dataPath));
    
    % Increment the counter
    c = 1;
    clear contrasts;
    clear postRecepContrasts;
    clear luminance;
    clear chromaticity;
    
    % Iterate over the folders
    for f = 1:length(theFolders)
        if ~strcmpi(theFolders(f).name, '.') && ~strcmpi(theFolders(f).name, '..')
            fprintf('>> Validation %s\n', theFolders(f).name);
            
            % Go to the folder
            if isdir(fullfile(dropboxBasePath, dataPath, theFolders(f).name))
                cd(fullfile(dropboxBasePath, dataPath, theFolders(f).name));
            end
            
            % Find the only MAT file there is going to be
            theMATFile = dir([pwd '/*.mat']);
            
            if ~isempty(theMATFile)
                % Load the MAT file
                tmp = load(theMATFile.name);
                
                % Extract the infomation
                S = tmp.cals{1}.describe.cache.data(32).describe.S;
                wls = SToWls(S);
                observerAgeInYrs = tmp.cals{1}.describe.cache.REFERENCE_OBSERVER_AGE;
                fractionBleached = tmp.cals{1}.describe.cache.data(observerAgeInYrs).describe.fractionBleached;
                pupilDiameterMm = tmp.cals{1}.describe.cache.data(observerAgeInYrs).describe.params.pupilDiameterMm;
                fieldSizeDegrees = tmp.cals{1}.describe.cache.data(observerAgeInYrs).describe.params.fieldSizeDegrees;
                
 
                bgSpd = tmp.cals{1}.modulationBGMeas.meas.pr650.spectrum; hold on;
                modSpd = tmp.cals{1}.modulationMaxMeas.meas.pr650.spectrum;

                plot(wls, bgSpd);
                plot(wls, modSpd);
                
                % Calculate luminance and chromaticity
                % Load the CIE functions
                load T_xyz1931
                T_xyz = SplineCmf(S_xyz1931,683*T_xyz1931,WlsToS(wls));
                
                % Calculate luminance and chromaticy
                luminance(c) = T_xyz(2, :)*bgSpd;
                chromaticity(:, c) = (T_xyz([1 2], :)*bgSpd)/sum((T_xyz*bgSpd));
                
                
                
                % Set up the receptor object
                receptorObj = SSTReceptorHuman('obsAgeYrs', observerAgeInYrs, 'fieldSizeDeg', fieldSizeDegrees, 'obsPupilDiameterMm', pupilDiameterMm);
                T_rec = receptorObj.T.T_energyNormalized;
                T_val = tmp.cals{1}.describe.cache.data(observerAgeInYrs).describe.T_receptors;
                
                % Calculate the numerical difference between the assumed and the
                % reconstructed receptor sensitvities
                for ii = 1:size(T_rec, 1)-1
                    tmp1 = sum(T_rec(1, :) - T_val(1, :));
                    if tmp1 > 0
                        error('Error: Couldn''t reconstruct receptor sensitivities...');
                    end
                end
                
                % Calculate the nominal contrast
                for jj = 1:size(T_rec, 1)
                    contrasts(:, c) = (T_rec*(modSpd-bgSpd))./(T_rec*bgSpd);
                end
                postRecepContrasts(:, c) = [1 1 1 ; 1 -1 0]' \ contrasts(:, c);
                
                
                % Increment the counter
                c = c+1;
            end
            
        end
    end
    % Take the average
    lumMean = mean(luminance);
    lumSD = std(luminance);
    chromMean = mean(chromaticity, 2);
    chromSD = std(chromaticity, [], 2);
    
    contrastsMean = mean(contrasts, 2);
    contrastsSD = std(contrasts, [], 2);
    postRecepContrastsMean = mean(postRecepContrasts, 2);
    postRecepContrastsSD = std(postRecepContrasts, [], 2);
    
    M = [];
    M = [M lumMean lumSD chromMean(1) chromSD(1) chromMean(2) chromSD(2)];
    for m = 1:length(contrastsMean)
        M = [M 100*contrastsMean(m) 100*contrastsSD(m)];
    end
    for m = 1:length(postRecepContrastsMean)
        M = [M 100*postRecepContrastsMean(m) 100*postRecepContrastsSD(m)];
    end
    
    fprintf(fid, '%s,%s,', theStimuli{d}, theObservers{d});
    for ii = 1:length(M);
        fprintf(fid, '%.2f,', M(ii));
    end
    fprintf(fid, '\n');
    
end
cd(currDir);
fclose(fid);