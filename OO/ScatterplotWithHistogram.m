function ScatterPlotWithHistogram(x, y, varargin)
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
p.addParameter('PlotMean', false, @islogical);
p.addParameter('PlotMedian', false, @islogical);
p.addParameter('MaxP', [], @isnumeric);
p.addParameter('PlotMarginals', true, @islogical);

p.KeepUnmatched = true;
p.parse(varargin{:});

% Set axis limits
minx = p.Results.XLim(1);
maxx = p.Results.XLim(2);
miny = p.Results.YLim(1);
maxy = p.Results.YLim(2);

%% Y histogram
if p.Results.PlotMarginals
    ah1 = subplot(2, 2, 4);
    hold on;
    histy = histogram(y, 'Orientation','horizontal', 'Normalization', 'probability', 'BinWidth', p.Results.YBinWidth);
    histy.FaceColor = p.Results.Color(2, :);
    plot([0 0], [miny maxy], '-k');
    
    box off; pbaspect([1 1 1]);
    set(gca, 'TickDir', 'out');
    set(gca, 'YAxisLocation', 'right');
    xlabel('Probability');
    set(gca, 'Color', [0.9 0.9 0.9]);
    
    %% X histogram
    ah2 = subplot(2, 2, 1);
    hold on;
    histx = histogram(x, 'Normalization', 'probability', 'BinWidth', p.Results.XBinWidth);
    histx.FaceColor = p.Results.Color(1, :);
    plot([minx maxx], [0 0], '-k');
    
    box off; pbaspect([1 1 1]);
    set(gca, 'TickDir', 'out');
    set(gca, 'XAxisLocation', 'top');
    ylabel('Probability');
    set(gca, 'Color', [0.9 0.9 0.9]);
end

%% Scatterplot
if p.Results.PlotMarginals
    ah3 = subplot(2, 2, 3);
end
hold on;

% Add reference lines
if ~isempty(p.Results.XRefLines)
    plot(p.Results.XRefLines(1, :), p.Results.XRefLines(2, :), ':k');
end

if ~isempty(p.Results.YRefLines)
    plot(p.Results.YRefLines(1, :), p.Results.YRefLines(2, :), ':k');
end

% Add scatter plot
scatter(x, y, 8, p.Results.Color(1, :), '.');
scatter(x, y, 5, p.Results.Color(2, :), '.');
xlabel(p.Results.XLabel);
ylabel(p.Results.YLabel);
set(gca, 'TickDir', 'out');
box off; pbaspect([1 1 1]);

%% Fix the histogram
% Figure out the maximum robability and add a little headroom
if isempty(p.Results.MaxP)
    maxp =  max([histx.Values histy.Values])*1.2;
else
    maxp = p.Results.MaxP;
end

if p.Results.PlotMarginals
    % Add mean, median and nominal symbols
    subplot(2, 2, 1);
    if p.Results.PlotMean
        plot([mean(histx.Data) mean(histx.Data)], [0 1000], '-r');
    end
    if p.Results.PlotMedian
        plot([median(histx.Data) median(histx.Data)], [0 1000], '-k');
    end
    if ~isempty(p.Results.XNominalContrast)
        plot([p.Results.XNominalContrast p.Results.XNominalContrast], [0 1000], ':k');
    end
    
    subplot(2, 2, 4);
    if p.Results.PlotMean
        plot([0 1000], [mean(histy.Data) mean(histy.Data)], '-r');
    end
    if p.Results.PlotMedian
        plot([0 1000], [median(histy.Data) median(histy.Data)], '-k');
    end
    if ~isempty(p.Results.YNominalContrast)
        plot([0 1000], [p.Results.YNominalContrast p.Results.YNominalContrast], ':k');
    end
end

if p.Results.PlotMarginals
    subplot(2, 2, 3);
end
if p.Results.PlotMean
    plot(mean(histx.Data), mean(histy.Data), '+r');
end
if p.Results.PlotMedian
    plot(median(histx.Data), median(histy.Data), '+k');
end

if p.Results.PlotMarginals
    % Adjust the axis limits
    ah1.XLim = [0 maxp];
    ah1.YLim = [miny maxy];
    
    ah2.XLim = [minx maxx];
    ah2.YLim = [0 maxp];
    
    ah3.XLim = [minx maxx];
    ah3.YLim = [miny maxy];
else
    xlim([minx maxx]);
    ylim([miny maxy]);
end