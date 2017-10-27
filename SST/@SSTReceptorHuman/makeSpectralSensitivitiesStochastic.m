function makeSpectralSensitivitiesStochastic(obj, varargin)
% makeSpectralSensitivitiesStochastic
%
% Usage:
%     receptorObj.makeSpectralSensitivitiesStochastic
%
% Description:
%     This method allows for stochastic resampling of the individual
%     difference parameters of our photoreceptor model. The standard
%     deviations are taken from Asano et al. (2015) and corresponds to
%     estimates of the natural physiological variability of the
%     photoreceptor model.
%
%     The following individual difference parameters are varied
%       dlens - Deviation in % from CIE computed peak lens density
%       dmac - Deviation in % from CIE peak macular pigment density
%       dphotopigment - Vector of deviations in % from CIE photopigment peak density
%       lambdaMaxShift - Vector of values (in nm) to shift lambda max of each photopigment absorbance by
%
%     Parameters for observer age, field size, and pupil diameter are assumed known and held fixed.
%
%     The lambdaMaxShift is implemented along a linear or log wavelength axis.
%
%     The standard deviations used are given in Table 5 of Asano et al. (2016),
%     doi.org/10.1371/journal.pone.0145671.
%
%     The outputs are returned as cell arrays of spectral sensitivity matrices,
%     to the field "Ts" of the receptor object, in the following formats/units:
%       Ts{}.T_quantalIsomerizations - Quantal isomerizations
%       Ts{}.T_quantalAbsorptions - Quantal absorptions
%       Ts{}.T_quantalAbsorptionsNormalized - Normalized quantal absoprtions
%       Ts{}.T_energy - Energy fundamentals
%       Ts{}.T_energyNormalized - Normalized energy fundamentals
%
%     The spectral sensitivites are specified for light reaching the cornea
%     (i.e., these are receptor fundamentals) and are, in order, L cones, M
%     cones, S cones, Melanopsin, and Rods.  If the object's
%     'doPenumbralConesTrueFalse' field is true, then penumbral L*, M* and
%     S* sensitivities are tacked on at the end.
%
%     Note that the macular pigment parameter often lens to macular
%     transmittances of >1, which leads to out-of-bound sampling. In that case,
%     a note is printed and the sample is rejected.
%
% Input:
%     obj - The receptorObj (e.g. from @SSTReceptor or @SSTReceptorHuman)
%
% Output:
%     None.
%
% Optional key/value pairs:
%     'NSamples' - Number of samples to do (default 1000)
%
%     'RandStream' - String determining which seed should be used to reset
%                    the random number generator (default 'mrg32k3a').
%
% See also:
%     @SSTReceptorHuman, makeSpectralSensitivities,
%     makeSpectralSensitivitiesParametricVariation

% 7/25/17    ms       Commented.
% 9/8/17     dhb      Change reject text to indicate that we expect sometimes to reject a draw.
%            dhb      Comment out warning message - that reduces faith of the user that the code is doing what it should.
%            ms       Updated header comments

% Parse vargin for options passed here
p = inputParser;
p.addParameter('NSamples', 1000, @isnumeric);
p.addParameter('RandStream', 'mrg32k3a', @isstring); % RNG for the resampling
p.KeepUnmatched = true;
p.parse(varargin{:});
NSamples = p.Results.NSamples;
whichRandStream = p.Results.RandStream;

% Note also that this assumes the standard deviation is constant for all
% baseline parameters.
dlensSD = 18.7; % Given in %
dmaculaSD = 36.5; % Given in %
dLConeSD = 9.0; % Given in %
dMConeSD = 9.0; % Given in %
dSConeSD = 7.4; % Given in %
lMaxLConeSD = 2.0; % Given nm
lMaxMConeSD = 1.5; % Given nm
lMaxSConeSD = 1.3; % Given nm

% Create 8 independent random number generators
[s1, s2, s3, s4, s5, s6, s7, s8] = RandStream.create(whichRandStream, 'NumStreams', 8);

% Print out some info if the verbosity level is high
NPrintStep = 200;
if strcmp(obj.verbosity, 'high')
    fprintf('* Generating resampled observers... \n');
end

% Resample!
c = 1;
cr = 1;
while c <= NSamples
    % Print out some info if the verbosity level is high
    if strcmp(obj.verbosity, 'high')
        if mod(c, NPrintStep) == 0
            % Rudimentary status update
            fprintf('  %i/%i <strong>[%.2f%s]</strong>\n', c, NSamples, 100*c/NSamples, '%');
        end
    end
    
    % Sample the lens density
    indDiffParamsLMS.dlens = randn(s1)*dlensSD;
    indDiffParamsMel.dlens = indDiffParamsLMS.dlens;
    indDiffParamsRod.dlens = indDiffParamsMel.dlens;
    
    % Sample the macular pigment density
    indDiffParamsLMS.dmac = randn(s2)*dmaculaSD;
    indDiffParamsMel.dmac = indDiffParamsLMS.dmac;
    indDiffParamsRod.dmac = indDiffParamsMel.dmac;
    
    % Sample the optical pigment density
    indDiffParamsLMS.dphotopigment = [randn(s3)*dLConeSD randn(s4)*dMConeSD randn(s5)*dSConeSD];
    indDiffParamsMel.dphotopigment = [0];
    indDiffParamsRod.dphotopigment = [0];
    
    % Sample the shift in lambda max
    indDiffParamsLMS.lambdaMaxShift = [randn(s6)*lMaxLConeSD randn(s7)*lMaxMConeSD randn(s8)*lMaxSConeSD];
    indDiffParamsMel.lambdaMaxShift = [0];
    indDiffParamsRod.lambdaMaxShift = [0];
    indDiffParamsLMS.shiftType = 'linear';
    indDiffParamsMel.shiftType = 'linear';
    indDiffParamsRod.shiftType = 'linear';
    
    %% Call ComputeCIEConeFundamentals to get the spectral sensitivities
    try
        %% LMS cones
        [T_quantalAbsorptionsNormalizedLMS,T_quantalAbsorptionsLMS,T_quantalIsomerizationsLMS,adjIndDiffParamsLMS] = ComputeCIEConeFundamentals(obj.S,...
            obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,[],[],[], ...
            false,[],[],indDiffParamsLMS);
        
        %% Melanopsin
        [T_quantalAbsorptionsNormalizedMel,T_quantalAbsorptionsMel,T_quantalIsomerizationsMel,adjIndDiffParamsMel] = ComputeCIEMelFundamental(obj.S,...
            obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,indDiffParamsMel);
        
        %% Rods
        [T_quantalAbsorptionsNormalizedRod,T_quantalAbsorptionsRod,T_quantalIsomerizationsRod,adjIndDiffParamsRod] = ComputeCIERodFundamental(obj.S,...
            obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,indDiffParamsRod);
        
        %% L*M*S* (penumbral) cones, if required
        if obj.doPenumbralConesTrueFalse
            % We assume standard parameters here.
            source = 'Prahl';
            vesselOxyFraction = 0.85;
            vesselOverallThicknessUm = 5;
            trans_Hemoglobin = GetHemoglobinTransmittance(obj.S, vesselOxyFraction, vesselOverallThicknessUm, source);
            
            % Expand for the three cones
            trans_Hemoglobin = repmat(trans_Hemoglobin, 1, size(T_quantalAbsorptionsNormalizedLMS, 1));
            
            T_quantalAbsorptionsNormalizedLMSPenumbral = T_quantalAbsorptionsNormalizedLMS .* trans_Hemoglobin';
            T_quantalAbsorptionsNormalizedLMSPenumbral = bsxfun(@rdivide,T_quantalAbsorptionsNormalizedLMSPenumbral,max(T_quantalAbsorptionsNormalizedLMSPenumbral, [], 2));
            T_quantalAbsorptionsLMSPenumbral = T_quantalAbsorptionsLMS .* trans_Hemoglobin';
            T_quantalIsomerizationsLMSPenumbral = T_quantalIsomerizationsLMS .* trans_Hemoglobin';
        end
        
        
        %% Assemble the spectral sensitivities
        T_quantalIsomerizations = [T_quantalIsomerizationsLMS ; T_quantalIsomerizationsMel ; T_quantalIsomerizationsRod];
        T_quantalAbsorptionsNormalized = [T_quantalAbsorptionsNormalizedLMS ; T_quantalAbsorptionsNormalizedMel ; T_quantalAbsorptionsNormalizedRod];
        T_quantalAbsorptions = [T_quantalAbsorptionsLMS ; T_quantalAbsorptionsNormalizedMel ; T_quantalAbsorptionsNormalizedRod];
        
        %% Add the penumbral cones if required
        if obj.doPenumbralConesTrueFalse
            T_quantalAbsorptionsNormalized = [T_quantalAbsorptionsNormalized ; T_quantalAbsorptionsNormalizedLMSPenumbral];
            T_quantalAbsorptions = [T_quantalAbsorptions ; T_quantalAbsorptionsLMSPenumbral];
            T_quantalIsomerizations = [T_quantalIsomerizations ; T_quantalIsomerizationsLMSPenumbral];
        end
        
        %% Convert to energy
        T_energy = EnergyToQuanta(obj.S,T_quantalIsomerizationsNormalized')';
        T_energyNormalized = bsxfun(@rdivide,T_energy,max(T_energy, [], 2));
        
        % Assign the sub-fields in the "Ts" field
        obj.Ts{c}.T_quantalIsomerizations = T_quantalIsomerizations;
        obj.Ts{c}.T_quantalAbsorptions = T_quantalAbsorptions;
        obj.Ts{c}.T_quantalAbsorptionsNormalized = T_quantalAbsorptionsNormalized;
        obj.Ts{c}.T_energy = T_energy;
        obj.Ts{c}.T_energyNormalized = T_energyNormalized;
        obj.Ts{c}.indDiffParams = indDiffParamsLMS;
        adjIndDiffParamsLMS.lambdaMaxShift = indDiffParamsLMS.lambdaMaxShift;
        obj.Ts{c}.adjIndDiffParams = adjIndDiffParamsLMS;
        
        % Increment
        c = c+1;
    catch e
        fprintf('* Sampling not successful for sample %g. Rejecting this sample. It is expected that we will sometimes reject samples, given that we are drawing from normal distributions.\n', c);
        cr = cr + 1; % Add to the counter
    end
    
end
fprintf('* # of rejected samples: %g\n', cr);