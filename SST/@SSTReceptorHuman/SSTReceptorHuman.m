classdef SSTReceptorHuman < SSTReceptor
    % SSTReceptorHuman
    %
    %   receptor = SSTReceptorHuman(varargin)
    %
    % Inherits optional key/value pairs from parent class SSTReceptor.
    
    % Public, read-only properties.
    properties (SetAccess = private, GetAccess = public)
        obsAgeInYrs; % Age of the observer in years
        obsPupilDiameterMm; % Pupil diameter
        fieldSizeDeg; % Field size in degrees
        T;  % Fundamentals created using makeSpectralSensitivities
        Tp; % Fundamentals created using makeSpectralSensitivitiesParametricVariation
        Ts; % Fundamentals created using makeSpectralSensitivitiesStochastic
        labels; % Labels for the spectral sensitivities
        MD5Hash; % MD5 hash of the receptor object
        doPenumbralConesTrueFalse;
    end
    
    % Private properties. Only methods of the parent class can set these
    properties(Access = private)
    end
    
    % Public methods
    methods
    end
    
    properties (Dependent)
    end
    
    % Methods.  Most public methods are implemented in a separate
    % function, but we put the class constructor here.
    methods (Access=public)
        % Constructor
        function obj = SSTReceptorHuman(varargin)
            % Base class constructor
            obj = obj@SSTReceptor(varargin{:});
            
            % Parse vargin for options passed here
            p = inputParser;
            p.addParameter('obsAgeInYrs', 32, @isnumeric);
            p.addParameter('obsPupilDiameterMm', 3, @isnumeric);
            p.addParameter('fieldSizeDeg', 10, @isnumeric);
            p.addParameter('doPenumbralConesTrueFalse', false, @islogical);
            p.addParameter('verbosity', 'high', @ischar);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            % Display the parameters and list any ones which didn't get matched
            theUnmatchedFields = fields(p.Unmatched);
            if ~isempty(theUnmatchedFields)
                warning('* There are unmatched parameters.');
                % Print the valid parameters
                fprintf('* Valid parameters for this model are:\n');
                for ii = 1:length(p.Parameters)
                    fprintf('  <strong>%s</strong>\n', p.Parameters{ii});
                end
                fprintf('\n');
                
                fprintf('* Unmatched input parameters for this model are:\n');
                for ii = 1:length(theUnmatchedFields)
                    fprintf('  <strong>%s</strong>\n', theUnmatchedFields{ii});
                end
                fprintf('\n');
            end
            
            % Assign the parameter values from the parser
            obj.obsAgeInYrs = p.Results.obsAgeInYrs;
            obj.obsPupilDiameterMm = p.Results.obsPupilDiameterMm;
            obj.fieldSizeDeg = p.Results.fieldSizeDeg;
            obj.doPenumbralConesTrueFalse = p.Results.doPenumbralConesTrueFalse;
            obj = makeSpectralSensitivities(obj);

            % Print out some output
            if strcmp(obj.verbosity, 'high')
                fprintf('* Setting up receptor object with parameters:\n');
                fprintf('  <strong>Age [yrs]</strong>:\t\t\t%i\n', obj.obsAgeInYrs);
                fprintf('  <strong>Pupil diameter [mm]</strong>:\t\t%.2f\n', obj.obsPupilDiameterMm);
                fprintf('  <strong>Field size [deg]</strong>:\t\t%.2f\n', obj.fieldSizeDeg);
                if obj.doPenumbralConesTrueFalse == 1
                    fprintf('  <strong>Including penumbral cones?</strong>:\t%s\n\n', 'True');
                elseif obj.doPenumbralConesTrueFalse == 0
                    fprintf('  <strong>Including penumbral cones?</strong>:\t%s\n\n', 'False');
                end
            end
        end
        makeSpectralSensitivitiesStochastic(obj, varargin);
        [obj, parv, parvlabel, parvlabellong, parvreal] = makeSpectralSensitivitiesParametricVariation(obj, varargin);
        setMD5Hash(obj);
    end
    
    % Get methods for dependent properties
    methods
    end
    
    % Methods may be called by the subclasses, but are otherwise private
    methods (Access = protected)
    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
    end
end