function f1 = plotSpectralSensitivities(obj, varargin)
% plotSpectralSensitivities(obj, varargin)
%
% Usage:
%     f1 = plotSpectralSensitivities(obj, varargin)
%
% Description:
%     This method plots the spectral sensitivities of a receptorObj. The
%     color for each spectral sensitivity is determined using the function
%     SSTDefaultReceptorColors(), which reads the labels in the receptor
%     object and returns predefined RGB triplets.
%
% Input:
%     obj - The receptorObj (e.g. from @SSTReceptor or @SSTReceptorHuman)
%
% Output:
%     f1 - Figure handle to the created plot.
%
% Optional key/value pairs:
%     'ax' - Axes object to plot in, rather than gca. If none is provided,
%            creates a new figure window and new axes object to plot in.
%
%     'NewWindow' - LEGACY, replaced by 'ax'. If passed instead of 'ax',
%                   functions as expected. Ignored if both are passed.
%                   Logical determing whether to  create a new window or to
%                   keep on plotting in gcf. Can be true or false. 
%                   (Default: true)
%
%     'whichFormat' - String determining which field of the receptorObj
%                     gets plotted. Possible options are:
%                           T_quantalIsomerizations
%                           T_quantalAbsorptions
%                           T_quantalAbsorptionsNormalized
%                           T_energy
%                           T_energyNormalized
%                     (Default: T_energyNormalized)
%
%     'logUnits' - Logical determining whether the quantity specified in
%                  'whichFormat' should be plotted on a log scale.
%                  (Default: false)
%
%     'saveFig' - Logical determining whether the figure should be saved.
%                 (Default: true)
%
%     'saveFigPath' - String determining the path where the figure should
%                     be saved. The default is empty ('') so that it will
%                     just get saved in the current directory.
%
% See also:
%     @SSTReceptorHuman, SSTDefaultReceptorColors

% 9/8/17  ms  Added header comments.

%% Parse vargin for options passed here
%
% Setting 'KeepUmatched' to true means that we can pass the varargin{:})
% along from a calling routine without an error here, if the key/value
% pairs recognized by the calling routine are not needed here.
p = inputParser; p.KeepUnmatched = true;
p.addParameter('ax',[], @ishandle);
p.addParameter('NewWindow', true, @islogical);
p.addParameter('whichFormat', 'T_energyNormalized', @ischar);
p.addParameter('logUnits', false, @islogical);
p.addParameter('saveFig', true, @islogical);
p.addParameter('saveFigPath', '', @ischar);
p.parse(varargin{:});

% Extract some parameters
whichFormat = p.Results.whichFormat;
saveFigYesNo = p.Results.saveFig;
saveFigPath = p.Results.saveFigPath;
if isempty(saveFigPath)
    saveFigPath = fullfile(['SpectralSensitivities_' whichFormat '.png']);
else
    saveFigPath = fullfile(saveFigPath, ['SpectralSensitivities_' whichFormat '.png']);
end

% Determine what we want to plot based on the varargin
switch whichFormat
    case 'T_quantalIsomerizations'
        T_receptors = obj.T.T_quantalIsomerizations;
        sublabel = 'Quantal isomerizations';
        yLabel = 'Quantal sensitivity';
    case 'T_quantalAbsorptions'
        T_receptors = obj.T.T_quantalAbsorptions;
        sublabel = 'Quantal absorptions';
        yLabel = 'Quantal absorptions';
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
if isempty(p.Results.ax)
    if p.Results.NewWindow
        f1 = figure(); clf; hold on;
    else
        f1 = gcf;
    end
    ax = gca;
else
    ax = p.Results.ax;
    axes(ax);
    f1 = gcf;
end

% Plot
% Determine the wavelength spacing
wls = SToWls(obj.S);

% Get the colors
theRGB = SSTDefaultReceptorColors(obj.labels);

% Find out xy limit
xLims = [wls(1) wls(end)];
yLims = [min(min(T_receptors))*0.9  max(max(T_receptors))*1.1];

% Determine how many receptors we need to plot
NReceptorsToPlot = length(obj.labels);
for ii = 1:NReceptorsToPlot
    % Plot the fundamentals
    h(ii) = plot(ax,wls, T_receptors(ii, :), '-', 'LineWidth', 2, 'Color', theRGB(ii, :)); hold on;
end

% Tune the plot properties
title(sublabel);
legend(h, obj.labels); legend boxoff;
xlim(xLims); ylim(yLims);
ylabel(yLabel);
pbaspect([1 1 1]); set(gca, 'TickDir', 'out');

if saveFigYesNo
    % Save out the figures
    set(f1, 'PaperPosition', [0 0 15 15]);
    set(f1, 'PaperSize', [15 15]);
    set(f1, 'Color', 'w');
    set(f1, 'InvertHardcopy', 'off');
    saveas(f1, saveFigPath, 'png');
end