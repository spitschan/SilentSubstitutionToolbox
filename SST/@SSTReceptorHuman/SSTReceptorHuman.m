classdef SSTReceptorHuman < SSTReceptor
    % SSTReceptorHuman
    %
    % Usage:
    %     receptorObj = SSTReceptorHuman;
    %
    % Description:
    %     @SSTReceptor is the class for the human photoreceptor spectral
    %     sensitivities. Upon instantiating an object of this class, a
    %     standard set of nominal spectral sensitivities will be generated
    %     using this class's makeSpectralSensitivities function.
    %
    %     The spectral sensitivites are specified for light reaching the
    %     cornea (i.e., these are receptor fundamentals) and are, in order,
    %     L cones, M cones, S cones, Melanopsin, and Rods.  If the object's
    %     'doPenumbralConesTrueFalse' field is true, then penumbral L*, M*
    %     and S* sensitivities are tacked on at the end.  These in up as
    %     rows in the separate fields of the various T... (e.g. T, Tp,
    %     Tp_i, Ts) properites of the object, where each T... is a struct.
    %     The labels property of the object provides string labels for each
    %     receptor type.
    %   
    %     When fundamentals are computed, we get them in various units as
    %     fields of the corresponding struct.  The individual T_ matrices
    %     follow Psychtoolbox conventions.  Each row is the sensitivity for
    %     one class, and each column represents one wavelength band.  In
    %     this form, you can multiply T_ into a column vector that
    %     describes the radiance and get back a column vector that contains
    %     the receptor responses.
    %
    %        T_quantalIsomerizations - light specified as radiance with
    %        power in quanta/sec per wavelength band and result is the
    %        number of isomerizations per second.  The geometric units of
    %        radiance are [m2-sr].
    %
    %        T_quantalAbsorptions - light specific as radiance in quantal
    %        units as just above, but result is the number of quanta
    %        absorbed per second.  Not all absorptions lead to an isomerization, so
    %        this is different from isomerizations in our nomenclature.
    %        The two differ by a factor for each receptor type which is called
    %        the quantal efficiency.  [NOTE: DHB - I can't imagine why
    %        anyone would want to actually use this, if I correctly
    %        understand what it is.  It may be a mistake to expose it,
    %        because the most likely effect is to cause confusion.  Some
    %        people use "absorptions" as a synonym for "isomerizations",
    %        which is one reason I wrote this long comment but as still
    %        worried.]
    %
    %        T_quantalAbsorptionsNormalized: These are like the quantal
    %        absorptions but with each sensitivity normlizeed to a peak of
    %        1. Since they are normalized, what you get out doesn't really
    %        correspond to anything physical, and we might call the results
    %        cone coordinates.  The coordinates will predict color matches
    %        as long as the input was in quantal units.  [NOTE: DHB - I am
    %        not sure we want to encourage use of this either, and might
    %        delete it as a propery.]
    %
    %        T_energy: These are fundamentals that expect radiance to be
    %        specified in terms of Watts, rather than quanta/sec, but
    %        none-the-less return the number of isomerizations per second.
    %
    %        T_energyNormalized: These are fundamentals each normalized to
    %        a max of 1.  This is often how tabulated fundamentals are
    %        provided.  If you use these, you get cone coordinates where
    %        the numbers are on a human scale and predict color matches,
    %        but there is no direct relation to the biopysically
    %        meaningful quantity of isomerizations/sec.  
    %
    % Input:
    %     None.
    %
    % Output:
    %     receptorObj - The receptor object
    % 
    % The receptor object has the following output fields at instantiation:
    %
    %     receptorObj.obsAgeInYrs - Age of the observer in years
    %     receptorObj.obsPupilDiameterMm - Observer pupil diameter
    %     receptorObj.T - Fundamentals created using the makeSpectralSensitivities method
    %     receptorObj.labels - Labels for the spectral sensitvities contained in receptorObj.T
    %     receptorObj.doPenumbralConesTrueFalse - Compute penumbral cone sensitivities?
    % 
    % Additional fields are created using methods of this class. These are:
    %     
    %     receptorObj.Tp - Fundamentals created using makeSpectralSensitivitiesParametricVariation.
    %                      These describe variation in steps along the parameters of the Asano et al.
    %                      model.
    %
    %     receptorObj.Tp_i - Information about the parametric variations
    %
    %     receptorObj.Ts - Fundamentals created using makeSpectralSensitivitiesStochastic.  These
    %                      represent draws from the statistical model of Asano et al.
    %
    %     receptorObj.MD5Hash - MD5 checksum of the receptor object. Useful to check for integrity of specific
    %                           resampled receptor sets.
    %
    % Optional key/value pairs:
    %     'obsAgeInYrs' - Observer age in years (Default: 32 years).
    %
    %     'obsPupilDiameter' - Assumed pupil diameter of the observer
    %                          (Default: 3 mm)
    %
    %     'fieldSizeDeg' - Assume field size (Default: 10 deg)
    %
    %     'doPenumbralConesTrueFalse' - Logical determining whether the
    %                                   penumbral cone spectral
    %                                   sensitivities should also be
    %                                   generated.
    %
    %     'verbosity' - Verbosity level for printing diagnostics. Possible
    %                   options:
    %                       'none' - No diagnostic print-out
    %                       'high' - Print everything that might be useful
    %
    %     'S' - The wavelength sampling specification in standard PTB
    %           notation. (Default: [380 2 201]).
    %
    % Methods that are implemented in this subclasses:
    %     makeSpectralSensitivities
    %     makeSpectralSensitivitiesParametricVariation
    %     makeSpectralSensitivitiesStochastic
    %     setMD5Hash
    %       
    % Methods implemented in the parent class SSTReceptor:
    %     plotSpectralSensitivities - Plots the spectral sensitivities
    %                                 generated by any subclass object
    %
    % See also:
    %     @SSTReceptor, makeSpectralSensitivities,
    %     makeSpectralSensitivitiesStochastic,
    %     makeSpectralSensitivitiesParametricVariation, setMD5Hash
    
    % 9/8/17  ms  Added header comments.
    
    % Public, read-only properties.
    properties (SetAccess = private, GetAccess = public)
        obsAgeInYrs; % Age of the observer in years
        obsPupilDiameterMm; % Pupil diameter
        fieldSizeDeg; % Field size in degrees
        T;  % Fundamentals created using makeSpectralSensitivities
        Tp; % Fundamentals created using makeSpectralSensitivitiesParametricVariation
        Tp_i; % Information about the parametric variation from makeSpectralSensitivitiesParametricVariation 
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
            p.addParameter('S',[380 2 201],@isnumeric);
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
            obj.makeSpectralSensitivities;
            
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
        
        % Declare the rest of the methods
        makeSpectralSensitivitiesStochastic(obj, varargin);
        [parv, parvlabel, parvlabellong, parvreal] = makeSpectralSensitivitiesParametricVariation(obj, varargin); 
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