classdef SSTReceptorHuman < SSTReceptor;
    % SSTReceptorHuman
    %
    %   receptor = SSTReceptorHuman(varargin)
    %
    % Inherits optional key/value pairs from parent class SSTReceptor.
    
    % Public, read-only properties.
    properties (SetAccess = private, GetAccess = public)
        obsAgeInYrs;
        obsPupilDiameterMm;
        fieldSizeDeg;
        T;
        Ts;
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
            p.addParameter('obsPupilDiameterMm', 4, @isnumeric);
            p.addParameter('fieldSizeDeg', 10, @isnumeric);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            obj.obsAgeInYrs = p.Results.obsAgeInYrs;
            obj.obsPupilDiameterMm = p.Results.obsPupilDiameterMm;
            obj.fieldSizeDeg = p.Results.fieldSizeDeg;
            obj = makeSpectralSensitivities(obj);
        end
        makeSpectralSensitivitiesStochastic(obj, varargin);
        makeSpectralSensitivitiesParametricVariation(obj, varargin);
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