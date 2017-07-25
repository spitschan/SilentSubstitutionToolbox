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
        T;  % Cone fundamentals created using makeSpectralSensitivities
        Tp; % Cone fundamentals created using makeSpectralSensitivitiesParametricVariation
        Ts; % Cone fundamentals created using makeSpectralSensitivitiesStochastic
        MD5Hash; % MD5 hash of the receptor object
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
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            obj.obsAgeInYrs = p.Results.obsAgeInYrs;
            obj.obsPupilDiameterMm = p.Results.obsPupilDiameterMm;
            obj.fieldSizeDeg = p.Results.fieldSizeDeg;
            obj = makeSpectralSensitivities(obj);
            
            if strcmp(obj.verbosity, 'high')
                fprintf('* Setting up receptor object with parameters:\n');
                fprintf('  <strong>Age [yrs]</strong>:\t\t%i\n', obj.obsAgeInYrs);
                fprintf('  <strong>Pupil diameter [mm]</strong>:\t%.2f\n', obj.obsPupilDiameterMm);
                fprintf('  <strong>Field size [deg]</strong>:\t%.2f\n\n', obj.fieldSizeDeg);
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