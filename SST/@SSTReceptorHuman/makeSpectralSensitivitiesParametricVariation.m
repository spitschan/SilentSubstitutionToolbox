function makeSpectralSensitivitiesParametricVariation(obj, varargin)
% makeSpectralSensitivitiesParametricVariation
%
% Usage:
%     receptorObj.makeSpectralSensitivitiesParametricVariation;
%
% Description:
%     This method of @SSTReceptorHuman creates generates spectral
%     sensitivities which differ in one of the parameters of the
%     8-parameter Asano et al. model.
%
%     The ranges are:
%       ±50% for dlens, dmac, dphotopigment{L|M|S}
%       ±2 nm for lambdaMaxShift{L|M|S}
%       2 to 9 mm for obsPupilDiameterMm.
%    
%     The outputs are returned as sets of receptor matrices in the cell
%     array field "Tp" of the receptor object, in the following
%     formats/units:
%       Tp{}.T_quantalIsomerizations - Quantal isomerizations
%       Tp{}.T_quantalAbsorptions - Quantal absorptions
%       Tp{}.T_quantalAbsorptionsNormalized - Normalized quantal absoprtions
%       Tp{}.T_energy - Energy fundamentals
%       Tp{}.T_energyNormalized - Normalized energy fundamentals
%
%     Each index item in Tp{} also contains the fields indDiffParams and
%     adjIndDiffParams, which correspond to the set of individual
%     difference parameters for that specific set of spectral
%     sensitivities.
%
%     In addition, which parameter values are varied in a given cell are
%     saved in the field Tp{}.parameterVariation, which contains the
%     following information:
%           Tp{}.parameterVariation.labelShort - (short) string telling us which parameter is varied
%           Tp{}.parameterVariation.label - string telling us which parameter is varied
%           Tp{}.parameterVariation.value - value of the parameter (e.g. -20% in dlens)
%           Tp{}.parameterVariation.adjValue - actual lens density value corresponding to the change in Tp{}.parameterVariation.value
%
%     This adds a little bit of redundancy, but helps in interpreting the
%     parameter changes. This is because it is not straightforward to know,
%     from the indDiffParams and adjIndDiffParams which parameters where
%     adjusted.
%
% Input:
%     obj - The receptorObj (e.g. from @SSTReceptor or @SSTReceptorHuman)
%
% Output:
%     None.
%
% Optional key/value pairs:
%     'whichParameter' - Determines which parameter which should be varied.
%                        Possible options are (default, 'dlens'):
%                           'dlens' - lens density
%                           'dmac' - macular pigment density
%                           'dphotopigmentL' - L cone photopigment density
%                           'dphotopigmentM' - M cone photopigment density
%                           'dphotopigmentS' - S cone photopigment density
%                           'lambdaMaxShiftL' - L lambda-max shift
%                           'lambdaMaxShiftM' - M lambda-max shift
%                           'lambdaMaxShiftS' - S lambda-max shift
%                           'obsPupilDiameterMm' - observer's pupil diameter 
%
%     'NTitrations' - The step sizes for each parameter variation, i.e. how
%                     many individual parameter values are calculated.
%                     (Default: 50)
%
% See also:
%     @SSTReceptorHuman, makeSpectralSensitivities,
%     makeSpectralSensitivitiesStochastic

% 7/25/17   ms  Commented.
% 9/7/17    ms  Updated header comments.

% Parse vargin for options passed here
p = inputParser;
p.addParameter('whichParameter', 'dlens', @ischar);
p.addParameter('NTitrations', 50, @isnumeric);
p.KeepUnmatched = true;
p.parse(varargin{:});
NTitrations = p.Results.NTitrations;

% The following individual difference parameters are supported:
%   indDiffParams.dlens - Deviation in % from CIE computed peak lens density
%   indDiffParams.dmac - Deviation in % from CIE peak macular pigment density
%   indDiffParams.dphotopigment - Vector of deviations in % from CIE photopigment peak density
%   indDiffParams.lambdaMaxShift - Vector of values (in nm) to shift lambda max of each photopigment absorbance by
%   indDiffParams.shiftType - 'linear' (default) or 'log'

% Fill the individual differences struct with values.
indDiffParams.dlens = 0;
indDiffParams.dmac = 0;
indDiffParams.dphotopigment = [0 0 0];
indDiffParams.lambdaMaxShift = [0 0 0];
indDiffParams.shiftType = 'linear';

if strcmp(obj.verbosity, 'high')
    % Print out some info if verbosity is high
    fprintf('* Generating parametric variation in <strong>%s</strong>... ', p.Results.whichParameter);
end

% Sample
for ii = 1:NTitrations
    switch p.Results.whichParameter
        case 'dlens' % Lens density
            parv = linspace(-50, 50, NTitrations);
            indDiffParams = DefaultIndDiffParams;
            indDiffParams.dlens = parv(ii);
            parvlabel = '%\DeltaD_{lens}';
            parvlabellong = 'Lens density';
            [T_quantalAbsorptionsNormalizedLMS,T_quantalAbsorptionsLMS,T_quantalIsomerizationsLMS,adjIndDiffParams] = ComputeCIEConeFundamentals(obj.S,...
                obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,[],[],[],[],[],[],indDiffParams);
            
            % Set to 0 for rod and mel
            indDiffParams.lambdaMaxShift = [0]; indDiffParams.dphotopigment = [0];
            [T_quantalAbsorptionsNormalizedMel,T_quantalAbsorptionsMel,T_quantalIsomerizationsMel,adjIndDiffParams] = ComputeCIEMelFundamental(obj.S,...
                obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,indDiffParams);
            [T_quantalAbsorptionsNormalizedRod,T_quantalAbsorptionsRod,T_quantalIsomerizationsRod,adjIndDiffParams] = ComputeCIERodFundamental(obj.S,...
                obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,indDiffParams);
            theIdx = 1;
        case 'dmac' % Macular pigment density
            parv = linspace(-50, 50, NTitrations);
            indDiffParams = DefaultIndDiffParams;
            indDiffParams.dmac = parv(ii);
            parvlabel = '%\DeltaD_{mac}';
            parvlabellong = 'Macular density';
            [T_quantalAbsorptionsNormalizedLMS,T_quantalAbsorptionsLMS,T_quantalIsomerizationsLMS,adjIndDiffParams] = ComputeCIEConeFundamentals(obj.S,...
                obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,[],[],[],[],[],[],indDiffParams);
            
            % Set to 0 for rod and mel
            indDiffParams.lambdaMaxShift = [0]; indDiffParams.dphotopigment = [0];
            [T_quantalAbsorptionsNormalizedMel,T_quantalAbsorptionsMel,T_quantalIsomerizationsMel,adjIndDiffParams] = ComputeCIEMelFundamental(obj.S,...
                obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,indDiffParams);
            [T_quantalAbsorptionsNormalizedRod,T_quantalAbsorptionsRod,T_quantalIsomerizationsRod,adjIndDiffParams] = ComputeCIERodFundamental(obj.S,...
                obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,indDiffParams);
            theIdx = 2;
        case 'dphotopigmentL' % Optical pigment density (L)
            parv = linspace(-50, 50, NTitrations);
            indDiffParams = DefaultIndDiffParams;
            indDiffParams.dphotopigment = [parv(ii) 0 0];
            parvlabel = '%\DeltaD_{pigment} [L]';
            parvlabellong =  {'Photopigment' ; 'optical density [L]'};
            [T_quantalAbsorptionsNormalizedLMS,T_quantalAbsorptionsLMS,T_quantalIsomerizationsLMS,adjIndDiffParams] = ComputeCIEConeFundamentals(obj.S,...
                obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,[],[],[],[],[],[],indDiffParams);
            [T_quantalAbsorptionsNormalizedMel,T_quantalAbsorptionsMel,T_quantalIsomerizationsMel,adjIndDiffParams] = ComputeCIEMelFundamental(obj.S,...
                obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,[]);
            [T_quantalAbsorptionsNormalizedRod,T_quantalAbsorptionsRod,T_quantalIsomerizationsRod,adjIndDiffParams] = ComputeCIERodFundamental(obj.S,...
                obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,[]);
            theIdx = 3;
        case 'dphotopigmentM' % Optical pigment density (M)
            parv = linspace(-50, 50, NTitrations);
            indDiffParams = DefaultIndDiffParams;
            indDiffParams.dphotopigment = [0 parv(ii) 0];
            parvlabel = '%\DeltaD_{pigment} [M]';
            parvlabellong =  {'Photopigment' ; 'optical density [M]'};
            [T_quantalAbsorptionsNormalizedLMS,T_quantalAbsorptionsLMS,T_quantalIsomerizationsLMS,adjIndDiffParams] = ComputeCIEConeFundamentals(obj.S,...
                obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,[],[],[],[],[],[],indDiffParams);
            [T_quantalAbsorptionsNormalizedMel,T_quantalAbsorptionsMel,T_quantalIsomerizationsMel,adjIndDiffParams] = ComputeCIEMelFundamental(obj.S,...
                obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,[]);
            [T_quantalAbsorptionsNormalizedRod,T_quantalAbsorptionsRod,T_quantalIsomerizationsRod,adjIndDiffParams] = ComputeCIERodFundamental(obj.S,...
                obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,[]);
            theIdx = 4;
        case 'dphotopigmentS' % Optical pigment density (S)
            parv = linspace(-50, 50, NTitrations);
            indDiffParams = DefaultIndDiffParams;
            indDiffParams.dphotopigment = [0 0 parv(ii)];
            parvlabel = '%\DeltaD_{pigment} [S]';
            parvlabellong =  {'Photopigment' ; 'optical density [S]'};
            [T_quantalAbsorptionsNormalizedLMS,T_quantalAbsorptionsLMS,T_quantalIsomerizationsLMS,adjIndDiffParams] = ComputeCIEConeFundamentals(obj.S,...
                obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,[],[],[],[],[],[],indDiffParams);
            [T_quantalAbsorptionsNormalizedMel,T_quantalAbsorptionsMel,T_quantalIsomerizationsMel,adjIndDiffParams] = ComputeCIEMelFundamental(obj.S,...
                obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,[]);
            [T_quantalAbsorptionsNormalizedRod,T_quantalAbsorptionsRod,T_quantalIsomerizationsRod,adjIndDiffParams] = ComputeCIERodFundamental(obj.S,...
                obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,[]);
            theIdx = 5;
        case 'lambdaMaxShiftL' % Shift in lambda-max (L)
            parv = linspace(-2, 2, NTitrations);
            indDiffParams = DefaultIndDiffParams;
            indDiffParams.lambdaMaxShift = [parv(ii) 0 0];
            parvlabel = '\Delta\lambda_{max} [L]';
            parvlabellong = {'Peak sensitivity' '\lambda_{max} [L]'};
            [T_quantalAbsorptionsNormalizedLMS,T_quantalAbsorptionsLMS,T_quantalIsomerizationsLMS,adjIndDiffParams] = ComputeCIEConeFundamentals(obj.S,...
                obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,[],[],[],[],[],[],indDiffParams);
            [T_quantalAbsorptionsNormalizedMel,T_quantalAbsorptionsMel,T_quantalIsomerizationsMel,adjIndDiffParams] = ComputeCIEMelFundamental(obj.S,...
                obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,[]);
            [T_quantalAbsorptionsNormalizedRod,T_quantalAbsorptionsRod,T_quantalIsomerizationsRod,adjIndDiffParams] = ComputeCIERodFundamental(obj.S,...
                obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,[]);
            theIdx = 6;
        case 'lambdaMaxShiftM' % Shift in lambda-max (M)
            parv = linspace(-2, 2, NTitrations);
            indDiffParams = DefaultIndDiffParams;
            indDiffParams.lambdaMaxShift = [0 parv(ii) 0];
            parvlabel = '\Delta\lambda_{max} [M]';
            parvlabellong = {'Peak sensitivity' '\lambda_{max} [M]'};
            [T_quantalAbsorptionsNormalizedLMS,T_quantalAbsorptionsLMS,T_quantalIsomerizationsLMS,adjIndDiffParams] = ComputeCIEConeFundamentals(obj.S,...
                obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,[],[],[],[],[],[],indDiffParams);
            [T_quantalAbsorptionsNormalizedMel,T_quantalAbsorptionsMel,T_quantalIsomerizationsMel,adjIndDiffParams] = ComputeCIEMelFundamental(obj.S,...
                obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,[]);
            [T_quantalAbsorptionsNormalizedRod,T_quantalAbsorptionsRod,T_quantalIsomerizationsRod,adjIndDiffParams] = ComputeCIERodFundamental(obj.S,...
                obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,[]);
            theIdx = 7;
        case 'lambdaMaxShiftS' % Shift in lambda-max (S)
            parv = linspace(-2, 2, NTitrations);
            indDiffParams = DefaultIndDiffParams;
            indDiffParams.lambdaMaxShift = [0 0 parv(ii)];
            parvlabel = '\Delta\lambda_{max} [S]';
            parvlabellong = {'Peak sensitivity' '\lambda_{max} [S]'};
            [T_quantalAbsorptionsNormalizedLMS,T_quantalAbsorptionsLMS,T_quantalIsomerizationsLMS,adjIndDiffParams] = ComputeCIEConeFundamentals(obj.S,...
                obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,[],[],[],[],[],[],indDiffParams);
            [T_quantalAbsorptionsNormalizedMel,T_quantalAbsorptionsMel,T_quantalIsomerizationsMel,adjIndDiffParams] = ComputeCIEMelFundamental(obj.S,...
                obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,[]);
            [T_quantalAbsorptionsNormalizedRod,T_quantalAbsorptionsRod,T_quantalIsomerizationsRod,adjIndDiffParams] = ComputeCIERodFundamental(obj.S,...
                obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,[]);
            theIdx = 8;
        case 'obsPupilDiameterMm' % Observer pupil diameter
            parv = linspace(2, 9, NTitrations);
            indDiffParams = DefaultIndDiffParams;
            parvlabel = 'Pupil diameter [mm]';
            parvlabellong = 'Pupil diameter';
            [T_quantalAbsorptionsNormalizedLMS,T_quantalAbsorptionsLMS,T_quantalIsomerizationsLMS,adjIndDiffParams] = ComputeCIEConeFundamentals(obj.S,...
                obj.fieldSizeDeg,obj.obsAgeInYrs,parv(ii),[],[],[],[],[],[],indDiffParams);
            [T_quantalAbsorptionsNormalizedMel,T_quantalAbsorptionsMel,T_quantalIsomerizationsMel,adjIndDiffParams] = ComputeCIEMelFundamental(obj.S,...
                obj.fieldSizeDeg,obj.obsAgeInYrs,parv(ii),[]);
            [T_quantalAbsorptionsNormalizedRod,T_quantalAbsorptionsRod,T_quantalIsomerizationsRod,adjIndDiffParams] = ComputeCIERodFundamental(obj.S,...
                obj.fieldSizeDeg,obj.obsAgeInYrs,parv(ii),[]);
            theIdx = 9;
    end
    
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
    
    %% Assemble the sensitivities
    % Normalized quantal sensitivities
    T_quantalAbsorptionsNormalized = [T_quantalAbsorptionsNormalizedLMS ; T_quantalAbsorptionsNormalizedMel ; T_quantalAbsorptionsNormalizedRod];
    
    % Quantal isomerizations
    T_quantalIsomerizations = [T_quantalIsomerizationsLMS ; T_quantalIsomerizationsMel ; T_quantalIsomerizationsRod];
    
    % Quantal absorption
    T_quantalAbsorptions = [T_quantalAbsorptionsLMS ; T_quantalAbsorptionsMel ; T_quantalAbsorptionsRod];
    
    %% Add the penumbral cones if required
    if obj.doPenumbralConesTrueFalse
        T_quantalAbsorptionsNormalized = [T_quantalAbsorptionsNormalized ; T_quantalAbsorptionsNormalizedLMSPenumbral];
        T_quantalAbsorptions = [T_quantalAbsorptions ; T_quantalAbsorptionsLMSPenumbral];
        T_quantalIsomerizations = [T_quantalIsomerizations ; T_quantalIsomerizationsLMSPenumbral];
    end
    
    % Convert to energy fundamentals
    T_energy = EnergyToQuanta(obj.S,T_quantalAbsorptionsNormalized')';
    
    % And normalize the energy fundamentals
    T_energyNormalized = bsxfun(@rdivide,T_energy,max(T_energy, [], 2));
    
    % Save out the parameter values
    switch p.Results.whichParameter
        case 'dlens'
            parvreal(:, ii) = adjIndDiffParams.lens;
        case 'dmac'
            parvreal(:, ii) = adjIndDiffParams.mac;
        case {'dphotopigmentL' 'dphotopigmentM' 'dphotopigmentS'}
            parvreal(:, ii) = adjIndDiffParams.dphotopigment;
        case {'lambdaMaxShiftL' 'lambdaMaxShiftM' 'lambdaMaxShiftS'}
            parvreal(:, ii) = indDiffParams.lambdaMaxShift;
        case 'obsPupilDiameterMm'
            parvreal(:, ii) = parv(ii);
    end
    
    % Fill in the sub-fields in the "Tp" field
    obj.Tp{theIdx, ii}.T_quantalAbsorptionsNormalized = T_quantalAbsorptionsNormalized;
    obj.Tp{theIdx, ii}.T_quantalAbsorptions = T_quantalAbsorptions;
    obj.Tp{theIdx, ii}.T_quantalIsomerizations = T_quantalIsomerizations;
    obj.Tp{theIdx, ii}.T_energy = T_energy;
    obj.Tp{theIdx, ii}.T_energyNormalized = T_energyNormalized;
    obj.Tp{theIdx, ii}.indDiffParams = indDiffParams;
    obj.Tp{theIdx, ii}.adjIndDiffParams = adjIndDiffParams;
    
    % Assign some additional info
    obj.Tp{theIdx, ii}.parameterVariation.labelShort = parvlabel;
    obj.Tp{theIdx, ii}.parameterVariation.label = parvlabellong;
    obj.Tp{theIdx, ii}.parameterVariation.value = parv;
    obj.Tp{theIdx, ii}.parameterVariation.adjValue = parvreal;
end

if strcmp(obj.verbosity, 'high')
    fprintf('Done.\n');
end