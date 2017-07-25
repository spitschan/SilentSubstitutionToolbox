function h = plotSpectralSensitivities(obj, varargin)
% plotSpectralSensitivities(obj, varargin)
%
% Plot the spectral sensitivities in rudimentary form.
%
% Key/value pairs
%   'NewWindow' - true/false (default true).  Create new window?
%
% 7/25/17   ms  Commented.

%% Parse vargin for options passed here
%
% Setting 'KeepUmatched' to true means that we can pass the varargin{:})
% along from a calling routine without an error here, if the key/value
% pairs recognized by the calling routine are not needed here.
p = inputParser; p.KeepUnmatched = true;
p.addParameter('NewWindow', true, @islogical);
p.parse(varargin{:});

%% Make the figure
if (p.Results.NewWindow)
    h = figure; clf; hold on;
else
    h = gcf;
end

% Plot
theRGB = DefaultReceptorColors(obj.labels);
wls = SToWls(obj.S);
NReceptorsToPlot = length(obj.labels);
for ii = 1:NReceptorsToPlot
    % Determine the peak wavelength
    [~, idx] = max(obj.T.T_energyNormalized(ii, :));
    
    plot(wls, obj.T.T_energyNormalized(ii, :), '-', 'LineWidth', 2, 'Color', theRGB(ii, :)); hold on;
end