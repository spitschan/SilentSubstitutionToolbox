function ScatterPlotNoHistogram(x, y, varargin)
% ScatterPlotWithHistogram(x, y, varargin)
%
% Maks a scatter lot with histograms in the margins.

% Parse vargin for options passed here
p = inputParser;
p.addParameter('XLim', [min(x) max(x)], @isnumeric);
p.addParameter('YLim', [min(y) max(y)], @isnumeric);
p.addParameter('XBinWidth', 0.001, @isnumeric);
p.addParameter('YBinWidth', 0.001, @isnumeric);
p.addParameter('XLabel', 'X', @ischar);
p.addParameter('YLabel', 'Y', @ischar);
p.addParameter('XRefLines', [], @isnumeric);
p.addParameter('YRefLines', [], @isnumeric);
p.addParameter('XNominalContrast', [], @isnumeric);
p.addParameter('YNominalContrast', [], @isnumeric)
p.addParameter('Color', [1 0 0 ; 0 0 1], @isnumeric);

p.KeepUnmatched = true;
p.parse(varargin{:});

% Set axis limits
minx = p.Results.XLim(1);
maxx = p.Results.XLim(2);
miny = p.Results.YLim(1);
maxy = p.Results.YLim(2);

%% Scatterplot
% Add reference lines
hold on;
plot(p.Results.XRefLines(1, :), p.Results.XRefLines(2, :), ':k');
plot(p.Results.YRefLines(1, :), p.Results.YRefLines(2, :), ':k');

% Add scatter plot
scatter(x, y, 8, p.Results.Color(1, :), '.');
scatter(x, y, 5, p.Results.Color(2, :), '.');

plot(mean(x), mean(y), '+r');
plot(median(x), median(y), '+k');

xlabel(p.Results.XLabel);
ylabel(p.Results.YLabel);
set(gca, 'TickDir', 'out');
box off; pbaspect([1 1 1]);

xlim([minx maxx]);
ylim([miny maxy]);