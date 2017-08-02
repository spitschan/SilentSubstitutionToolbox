function [obj, parv, parvlabel, parvlabellong, parvreal] = makeSpectralSensitivitiesParametricVariation(obj, varargin)
% [obj, parv, parvlabel, parvlabellong, parvreal] = makeSpectralSensitivitiesParametricVariation(obj, varargin)
%
% This method calculates spectral sensitivities along variation in the
% individual difference parameters.
%
% Currently, does not support the melanopsin and rod fundamentals.
%
% 7/25/17    ms       Commented.

% Parse vargin for options passed here
p = inputParser;
p.addParameter('WhichParameter', 'dlens', @ischar);
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
    fprintf('* Generating parametric variation in <strong>%s</strong>... ', p.Results.WhichParameter);
end

% Sample
for ii = 1:NTitrations
    switch p.Results.WhichParameter
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
    T_quantalAbsorptionsNormalized = [T_quantalAbsorptionsNormalizedLMS ; T_quantalAbsorptionsNormalizedMel ; T_quantalAbsorptionsNormalizedRod];
    T_quantalAbsorptions = [T_quantalAbsorptionsLMS ; T_quantalAbsorptionsMel ; T_quantalAbsorptionsRod];
    T_quantalIsomerizations = [T_quantalIsomerizationsLMS ; T_quantalIsomerizationsMel ; T_quantalIsomerizationsRod];
    
    % Get the cone fundamentals
    T_energy = EnergyToQuanta(obj.S,T_quantalAbsorptionsNormalized')';
    T_energyNormalized = bsxfun(@rdivide,T_energy,max(T_energy, [], 2));
    
    % Save out the parameter values
    switch p.Results.WhichParameter
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
end

if strcmp(obj.verbosity, 'high')
    fprintf('Done.\n');
end