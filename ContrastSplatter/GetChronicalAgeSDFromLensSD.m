function ageSD = GetChronicalAgeSDFromLensSD(src)
% ageSD = GetChronicalAgeSDFromLensSD(src)
%
% Returns the variability in age corresponding to a variability in standard
% deviations in lens density at a given wavelength.
%
% We assume standard wavelength spacing and a pupil size of 4.7 mm.
%
% Input:
%   src (str)   - Source of the estimate. 'Xu' and 'M&W' are allowable.
%
% Output:
%   ageSD (1x1) - SD in age units
%
% See also:
%   GetLambdaMaxEstimateSD
%
% References:
%   Webster MA & MacLeod DI (1988) Factors underlying individual
%   differences in the color matches of normal observers. J Opt Soc Am A
%   5(10):1722-1735. ('Xu')
% 
%   Xu J, Pokorny J, & Smith VC (1997) Optical density of the human lens. J
%   Opt Soc Am A Opt Image Sci Vis 14(5):953-960. ('M&W')
%
% 7/15/2014     ms          Wrote it.

switch src
    case 'Xu'
        % Load in the CSV file, which contains lens SD values from Xu,
        % Smith & Pokorny (1997), digitized by MS on July 21, 2014
        M = csvread(fullfile(fileparts(which(mfilename)), 'data', 'XuPokornySmith1997_Fig4A.csv'), 1);
        
        chronologicalAge = M(:, 1);
        predictedAge = M(:, 2);
        ageSD = round(std(chronologicalAge-predictedAge));
    case 'M&W'
        % Set up parameters
        S = [380 2 201];
        pupilSize = 4.7;
        lensSD = 0.18;
        wl= 400;
        
        % The following computations are taken from LensTransmittance.m (in PTB),
        % based on the CIE standard.
        
        % Load CIE age dependent and age independent components
        load den_lens_cie_1
        load den_lens_cie_2
        lensDensity1 = SplineSrf(S_lens_cie_1,den_lens_cie_1,S,2)';
        lensDensity2 = SplineSrf(S_lens_cie_2,den_lens_cie_2,S,2)';
        
        % For the reference density, there will be an associated SD given by the
        % error at 400 nm. This is given in M&W (0.18). Since CIE specifies a lens
        % spectrum and not a single lens density value, we base each calculations
        % on the 400 nm value. So, all lensDensity vectors below are indexed into
        % theWlIndex.
        theWlIndex = find(SToWls(S) == wl);
        
        % Invert the formulae that are below. For a 32 yr old observer, at 400 nm.
        refAge = 32;
        lensDensity = lensDensity1(theWlIndex)*(1+0.02*(refAge-32))+lensDensity2(theWlIndex)
        age = ((lensDensity - lensDensity2(theWlIndex))/lensDensity1(theWlIndex) - 1)./0.02 + 32
        
        % Now, we add the value in lensSD to the lensDensity, and figure
        % out what age that corresponds to.
        lensDensity1SDPlus = lensDensity + lensSD;
        ageSDPlus = ((lensDensity1SDPlus - lensDensity2(theWlIndex))/lensDensity1(theWlIndex) - 1)./0.02 + 32
        
        % Now, we subtract the value in lensSD from the lensDensity spectrum, and figure
        % out what age that corresponds to.
        lensDensity1SDMinus = lensDensity - lensSD;
        ageSDMinus = ((lensDensity1SDMinus - lensDensity2(theWlIndex))/lensDensity1(theWlIndex) - 1)./0.02 + 32
        
        ageSD = ageSDPlus-refAge;
        fprintf(['Lens optical density SD ' num2str(lensSD) ' translates to age SD ' num2str(ageSD)]);
end