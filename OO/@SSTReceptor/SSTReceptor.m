classdef SSTReceptor < handle
% The @SSTReceptor parent class for modeling spectral sensitivities.  This implements
% methods that are model independent and defines the API for subclasses
% that implement particular models.
%
% Methods implemented in this class:
%   plotSpectralSensitivities - Plot response and fit
%
% Methods that must be implemented in subclasses
%   makeSpectralSensitivities
%
% The class constructor can take optional key/value pairs as follows.
%      'verbosity' - string (default 'none').  Verbsoity level for method
%      diagnostic printout
%      'none' - No diagnostic printout
%      'high' - Print everything we can think might be useful.

    % Public read/write properties
    properties              
    end
    
    % Dependent properties, computed from other parameters
    properties
    end
        
    % Public, read-only properties.  These can be set by methods of the
    % parent class (that is, this class) but not by methods of subclasses.
    properties (SetAccess = private, GetAccess = public)
        % Verbosity level, set at constructor time.
        verbosity
        S
    end
    
    % Protected properties; Methods of the parent class and all of its
    % subclasses can set these.
    properties (SetAccess = protected, GetAccess = public)
    end
    
    % Private properties. Only methods of the parent class can set or read these
    properties (Access = private)      
    end
    
    % Public methods
    %
    % Methods defined in separate files are public by default, so we don't
    % explicitly decare them.
    %
    % We do stick the constructor method here.
    methods (Access=public)
        % Constructor
        function obj = SSTReceptor(varargin)
                     
            %% Parse vargin for options passed here
            p = inputParser;
            p.addParameter('verbosity','high',@ischar);
            p.addParameter('S',[380 2 201],@isnumeric);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            obj.verbosity = p.Results.verbosity;
            obj.S = p.Results.S;
        end
    end
    
    % Methods that must only be implemented in the subclasses.
    % If a subclass does not implement each and every of these methods
    % it cannot instantiate objects.
    methods (Abstract, Access=public)
        obj = makeSpectralSensitivities(obj);
    end
    
    % Methods may be called by the subclasses, but are otherwise private 
    methods (Access = protected)

    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)

    end
    
end
