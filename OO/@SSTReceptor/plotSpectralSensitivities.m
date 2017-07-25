function h = plotSpectralSensitivities(obj, varargin)
% plotSpectralSensitivities(obj, varargin)
%
% Plot the spectral sensitivities in rudimentary form.
%
% Key/value pairs
%   'NewWindow' - Logical determing whether to create a new window or to
%                 keep on plotting in gcf. Can be true or false.
%                 Default: true
%
%   'whichFormat' - String determining which types of fundamentals to plot.
%                   Possible ones are:
%                       T_quantalIsomerizations
%                       T_quantalAbsorptions
%                       T_quantalAbsorptionsNormalized
%                       T_energy
%                       T_energyNormalized
%                   Default: T_energyNormalized
%
%  'logUnits' - Logical determining whether plot should be in log units.
%
% 7/25/17   ms  Commented.

%% Parse vargin for options passed here
%
% Setting 'KeepUmatched' to true means that we can pass the varargin{:})
% along from a calling routine without an error here, if the key/value
% pairs recognized by the calling routine are not needed here.
p = inputParser; p.KeepUnmatched = true;
p.addParameter('NewWindow', true, @islogical);
p.addParameter('whichFormat', 'T_energyNormalized', @ischar);
p.addParameter('logUnits', false, @islogical);
p.parse(varargin{:});

% Determine what we want to plot based on the varargin
whichFormat = p.Results.whichFormat;
switch whichFormat
    case 'T_quantalIsomerizations'
        T_receptors = obj.T.T_quantalIsomerizations;
        sublabel = 'Quantal isomerizations';
        yLabel = 'Quantal sensitivity';
    case 'T_quantalAbsorptions'
        T_receptors = obj.T.T_quantalAbsorptions;
        sublabel = 'Quantal absorptions';
        yLabel = 'Quantal sensitivity';
    case 'T_quantalAbsorptionsNormalized'
        T_receptors = obj.T.T_quantalAbsorptionsNormalized;
        sublabel = 'Quantal absorptions [normalized]';
        yLabel = 'Relative sensitivity';
    case 'T_energy'
        T_receptors = obj.T.T_energy;
        sublabel = 'Energy fundamentals';
        yLabel = 'Energy sensitivity';
    case 'T_energyNormalized'
        T_receptors = obj.T.T_energyNormalized;
        sublabel = 'Energy fundamentals [normalized]';
        yLabel = 'Relative sensitivity';
    otherwise
        error(['Unknown type ' whichFormat']);
end

if p.Results.logUnits
   T_receptors = log10(T_receptors);
   yLabel = ['log_{10} ' lower(yLabel)];
end

%% Make the figure
if (p.Results.NewWindow)
    h = figure; clf; hold on;
else
    h = gcf;
end

% Plot
% Determine the wavelength spacing
wls = SToWls(obj.S);

% Get the colors
theRGB = DefaultReceptorColors(obj.labels);

% Find out xy maxima
xLims = [wls(1) wls(end)];
yLims = [min(min(T_receptors))*0.9  max(max(T_receptors))*1.1];

% Determine how many receptors we need to plot
NReceptorsToPlot = length(obj.labels);
for ii = 1:NReceptorsToPlot
    % Plot the fundamentals
    h(ii) = plot(wls, T_receptors(ii, :), '-', 'LineWidth', 2, 'Color', theRGB(ii, :)); hold on;
end
% Tune the plot properties
title(sublabel);
legend(h, obj.labels); legend boxoff;
xlim(xLims); ylim(yLims);
ylabel(yLabel);
pbaspect([1 1 1]); set(gca, 'TickDir', 'out');