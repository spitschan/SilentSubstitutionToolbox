function obj = makeSpectralSensitivitiesParametricVariation(obj, varargin)

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

switch p.Results.WhichParameter
    case 'dlens'
        theRange = linspace(-50, 50, NTitrations);
        
        % Sample
        for ii = 1:NTitrations
            indDiffParams = DefaultIndDiffParams;
            indDiffParams.dlens = theRange(ii);

            % Get the cone fundamentals
            [T_quantalAbsorptionsNormalized,T_quantalAbsorptions,T_quantalIsomerizations] = ComputeCIEConeFundamentals(obj.S,...
                obj.fieldSizeDeg,obj.obsAgeInYrs,obj.obsPupilDiameterMm,[],[],[], ...
                false,[],[],indDiffParams);
            T_energy = EnergyToQuanta(obj.S,T_quantalAbsorptionsNormalized')';
            T_energyNormalized = bsxfun(@rdivide,T_energy,max(T_energy, [], 2));
            
            % obj = makeSpectralSensitivities(obj)
            obj.Ts{ii}.T_quantalAbsorptionsNormalized = T_quantalAbsorptionsNormalized;
            obj.Ts{ii}.T_quantalAbsorptions = T_quantalAbsorptions;
            obj.Ts{ii}.T_quantalIsomerizations = T_quantalIsomerizations;
            obj.Ts{ii}.T_energy = T_energy;
            obj.Ts{ii}.T_energyNormalized = T_energyNormalized;
        end
end

% % The standard deviations are given in Table 5 in Asano et al. 2015
% dlensSD = 18.7;
% dmaculaSD = 10;%36.5;
% dLConeSD = 9.0;
% dMConeSD = 9.0;
% dSConeSD = 7.4;
% lMaxLConeSD = 2.0;
% lMaxMConeSD = 1.5;
% lMaxSConeSD = 1.3;
%
% % Create 8 random number generators
% [s1, s2, s3, s4, s5, s6, s7, s8] = RandStream.create('mrg32k3a','NumStreams',8);
%
