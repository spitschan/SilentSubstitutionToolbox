function [obj parv, parvlabel, parvlabellong, parvreal] = makeSpectralSensitivitiesParametricVariation(obj, varargin)

% Parse vargin for options passed here
p = inputParser;
p.addParameter('WhichParameter', 'dlens', @ischar);
p.addParameter('NTitrations', 50, @isnumeric);
p.KeepUnmatched = true;
p.parse(varargin{:});
NTitrations = p.Results.NTitrations;

% We use a few individual difference parameters
% indDiffParams.dlens - deviation in % from CIE computed peak lens density
% indDiffParams.dmac - deviation in % from CIE peak macular pigment density
% indDiffParams.dphotopigment - vector of deviations in % from CIE photopigment peak density.
% indDiffParams.lambdaMaxShift - vector of values (in nm) to shift lambda max of each photopigment absorbance by.
% indDiffParams.shiftType - 'linear' (default) or 'log'.

% Fill the individual differences struct with values.
indDiffParams.dlens = 0;
indDiffParams.dmac = 0;
indDiffParams.dphotopigment = [0 0 0];
indDiffParams.lambdaMaxShift = [0 0 0];
indDiffParams.shiftType = 'linear';

% Sample
for ii = 1:NTitrations
    switch p.Results.WhichParameter
        case 'dlens'
            parv = linspace(-50, 50, NTitrations);
            indDiffParams = DefaultIndDiffParams;
            indDiffParams.dlens = parv(ii);
            parvlabel = '%\DeltaD_{lens}';
            parvlabellong = 'Lens density';
            [T_quantalAbsorptionsNormalized,T_quantalAbsorptions,T_quantalIsomerizations,adjIndDiffParams] = ComputeCIEConeFundamentals(obj.S,...
                obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,[],[],[], ...
                false,[],[],indDiffParams);
            theIdx = 1;
        case 'dmac'
            parv = linspace(-50, 50, NTitrations);
            indDiffParams = DefaultIndDiffParams;
            indDiffParams.dmac = parv(ii);
            parvlabel = '%\DeltaD_{mac}';
            parvlabellong = 'Macular density';
            [T_quantalAbsorptionsNormalized,T_quantalAbsorptions,T_quantalIsomerizations,adjIndDiffParams] = ComputeCIEConeFundamentals(obj.S,...
                obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,[],[],[], ...
                false,[],[],indDiffParams);
            theIdx = 2;
        case 'dphotopigmentL'
            parv = linspace(-50, 50, NTitrations);
            indDiffParams = DefaultIndDiffParams;
            indDiffParams.dphotopigment = [parv(ii) 0 0];
            parvlabel = '%\DeltaD_{pigment}';
            parvlabellong =  {'Photopigment' ; 'optical density'};
            [T_quantalAbsorptionsNormalized,T_quantalAbsorptions,T_quantalIsomerizations,adjIndDiffParams] = ComputeCIEConeFundamentals(obj.S,...
                obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,[],[],[], ...
                false,[],[],indDiffParams);
            theIdx = 3;
        case 'dphotopigmentM'
            parv = linspace(-50, 50, NTitrations);
            indDiffParams = DefaultIndDiffParams;
            indDiffParams.dphotopigment = [0 parv(ii) 0];
            parvlabel = '%\DeltaD_{pigment}';
            parvlabellong =  {'Photopigment' ; 'optical density'};
            [T_quantalAbsorptionsNormalized,T_quantalAbsorptions,T_quantalIsomerizations,adjIndDiffParams] = ComputeCIEConeFundamentals(obj.S,...
                obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,[],[],[], ...
                false,[],[],indDiffParams);
            theIdx = 4;
        case 'dphotopigmentS'
            parv = linspace(-50, 50, NTitrations);
            indDiffParams = DefaultIndDiffParams;
            indDiffParams.dphotopigment = [0 0 parv(ii)];
            parvlabel = '%\DeltaD_{pigment}';
            parvlabellong =  {'Photopigment' ; 'optical density'};
            [T_quantalAbsorptionsNormalized,T_quantalAbsorptions,T_quantalIsomerizations,adjIndDiffParams] = ComputeCIEConeFundamentals(obj.S,...
                obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,[],[],[], ...
                false,[],[],indDiffParams);
            theIdx = 5;
        case 'lambdaMaxShiftL'
            parv = linspace(-2, 2, NTitrations);
            indDiffParams = DefaultIndDiffParams;
            indDiffParams.lambdaMaxShift = [parv(ii) 0 0];
            parvlabel = '\Delta\lambda_{max}';
            parvlabellong = {'Peak sensitivity' '\lambda_{max}'};
            [T_quantalAbsorptionsNormalized,T_quantalAbsorptions,T_quantalIsomerizations,adjIndDiffParams] = ComputeCIEConeFundamentals(obj.S,...
                obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,[],[],[], ...
                false,[],[],indDiffParams);
            theIdx = 6;
        case 'lambdaMaxShiftM'
            parv = linspace(-2, 2, NTitrations);
            indDiffParams = DefaultIndDiffParams;
            indDiffParams.lambdaMaxShift = [0 parv(ii) 0];
            parvlabel = '\Delta\lambda_{max}';
            parvlabellong = {'Peak sensitivity' '\lambda_{max}'};
            [T_quantalAbsorptionsNormalized,T_quantalAbsorptions,T_quantalIsomerizations,adjIndDiffParams] = ComputeCIEConeFundamentals(obj.S,...
                obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,[],[],[], ...
                false,[],[],indDiffParams);
            theIdx = 7;
        case 'lambdaMaxShiftS'
            parv = linspace(-2, 2, NTitrations);
            indDiffParams = DefaultIndDiffParams;
            indDiffParams.lambdaMaxShift = [0 0 parv(ii)];
            parvlabel = '\Delta\lambda_{max}';
            parvlabellong = {'Peak sensitivity' '\lambda_{max}'};
            [T_quantalAbsorptionsNormalized,T_quantalAbsorptions,T_quantalIsomerizations,adjIndDiffParams] = ComputeCIEConeFundamentals(obj.S,...
                obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,[],[],[], ...
                false,[],[],indDiffParams);
            theIdx = 8;
        case 'obsPupilDiameterMm'
            parv = linspace(3, 9, NTitrations);
            indDiffParams = DefaultIndDiffParams;
            parvlabel = 'Pupil diameter [mm]';
            parvlabellong = 'Pupil diameter';
            [T_quantalAbsorptionsNormalized,T_quantalAbsorptions,T_quantalIsomerizations,adjIndDiffParams] = ComputeCIEConeFundamentals(obj.S,...
                obj.fieldSizeDeg,obj.obsAgeInYrs,parv(ii),[],[],[], ...
                false,[],[],indDiffParams);
            theIdx = 9;
    end
    
    % Get the cone fundamentals
    T_energy = EnergyToQuanta(obj.S,T_quantalAbsorptionsNormalized')';
    T_energyNormalized = bsxfun(@rdivide,T_energy,max(T_energy, [], 2));
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
    
    % obj = makeSpectralSensitivities(obj)
    obj.Tp{theIdx, ii}.T_quantalAbsorptionsNormalized = T_quantalAbsorptionsNormalized;
    obj.Tp{theIdx, ii}.T_quantalAbsorptions = T_quantalAbsorptions;
    obj.Tp{theIdx, ii}.T_quantalIsomerizations = T_quantalIsomerizations;
    obj.Tp{theIdx, ii}.T_energy = T_energy;
    obj.Tp{theIdx, ii}.T_energyNormalized = T_energyNormalized;
    obj.Tp{theIdx, ii}.indDiffParams = indDiffParams;
    obj.Tp{theIdx, ii}.adjIndDiffParams = adjIndDiffParams;
end