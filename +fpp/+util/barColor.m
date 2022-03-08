
% [b,s,e] = fpp.util.barColor(X[,colorArray,errorVec,gapInd,dataPoints])
%
% Function to plot a bar graph with specific bar colors, and overlaid
% circles for individual data points.
%
% Arguments: 
% - X (numeric 2D matrix): samples by condition array of values to plot.
%       May contain NaN values, which are ignored.
% 
% Optional arguments:
% - colorArray (cell array): array of colors by condition
% - errorVec (scalar or vector): whether to use error bars, or vector value
%       for error. 0 = no error bars, 1 = standard error, vector = use
%       specified values.
% - gapInd (vector): gaps in the graph will be placed after these bars
% - dataPoints (boolean): whether to plot individual data points
% - dotArea (scalar): area of data points (size input to scatter)
%
% Outputs:
% - b (Bar object): bar graph handle
% - s (Scatter object): scatter plot handle
% - e (ErrorBar object): error bar graph handle, if errors were used
% - barInd (vector): y-position of bars in graph

function [b,s,e,barInd] = barColor(X,colorArray,errorVec,gapInd,dataPoints,dotArea)

e = [];
s = [];
nSamples = size(X,1);
nConds = size(X,2);

% Check inputs variables
if ~exist('colorArray','var') || isempty(colorArray)
   colorArray = {'r','b','g','m','c','k'};
elseif size(colorArray,1)>size(colorArray,2)
    colorArray = colorArray';
end
while length(colorArray)<nConds
    colorArray = [colorArray colorArray];
end
colorArray = colorArray(1:nConds);
if ~exist('errorVec','var') || isempty(errorVec)
    errorVec = 0;
elseif size(errorVec,1)>size(errorVec,2)
    errorVec = errorVec';
end
if ~exist('gapInd','var') || isempty(gapInd)
    gapInd = [];
end
if ~exist('dataPoints','var') || isempty(dataPoints)
    dataPoints = 1;
end
if ~exist('dotArea','var') || isempty(dotArea)
    dotArea = 50;
end

% Add gaps between bars
nGaps = length(gapInd);
nSlots = nConds+nGaps;
allInd = 1:nSlots;
gapIndNew = [];
for i=1:length(gapInd)
    gapIndNew(i) = gapInd(i)+i;
end
barInd = setdiff(allInd,gapIndNew);
XNew = zeros(nSamples,nSlots);
XNew(:,barInd) = X;
X = XNew;
colorArrayNew = repmat({'k'},[1 nSlots]);
colorArrayNew(barInd) = colorArray;
colorArray = colorArrayNew;

% Compute mean values, ignoring NaN vals
for i=1:size(X,2)
    meanVals(i) = mean(X(~isnan(X(:,i)),i),1);
end

% Determine error bar values
errors = [];
if ~isscalar(errorVec)
    errors = zeros(1,nSlots);
    errors(barInd) = errorVec;
    errorVec = 1;
elseif errorVec==1
    for i=1:size(X,2)
        errors(i) = std(X(~isnan(X(:,i)),i),1)/sqrt(sum(~isnan(X(:,1))));
    end
end

% Plot bar graph
b = bar(meanVals);
b.FaceColor = 'flat';

% Change bar colors
for c=1:size(b.CData,1)
    thisColor = colorArray{c};
    if ischar(thisColor)
        thisColor = rem(floor((strfind('kbgcrmyw', thisColor) - 1) * [0.25 0.5 1]), 2);
    elseif isscalar(thisColor)
        thisColor = repmat(thisColor,[1 3]);
    elseif size(thisColor,1)==3
        thisColor = thisColor';
    end
    b.CData(c,:) = thisColor;
end

% Plot individual data points
if dataPoints
    scatterX = []; scatterY = [];
    for i=barInd
        prevLength = length(scatterX);
        scatterX = [scatterX; i*ones(sum(~isnan(X(:,i))),1)];
        scatterY = [scatterY; X(~isnan(X(:,i)),i)];
        
        % Deal with multiple dots positioned at the same y-value. Stagger
        % x-position of dots to display multiple dots.
        allYVals = unique(X(:,i));
        if length(allYVals)<length(X(:,i))
            for yInd=1:length(allYVals)
                if sum(X(:,i)==allYVals(yInd))>1
                    ind = find(X(:,i)==allYVals(yInd));
                    xPos = linspace(i-b.BarWidth/2,i+b.BarWidth/2,length(ind)+2);
                    xPos = xPos(2:end-1);
                    for d=1:length(ind)
                        scatterX(prevLength+ind(d)) = xPos(d);
                    end
                end
            end
        end
    end
    hold on;
    s = scatter(scatterX,scatterY,dotArea,'k','filled');
end

% Plot errors
if errorVec
    hold on;
    e = errorbar(meanVals,errors,'.k');
    set(e,'LineWidth',2);
end

% Set figure/axis properties
box off;
set(gcf,'Color',[1 1 1]);
set(gca,'LineWidth',2,'FontSize',24,'XTick',[]);
set(b,'LineWidth',2,'BarWidth',.75);
set(b.BaseLine,'LineWidth',2);
ylim = get(gca,'YLim');
if ylim(1)<0    % Remove horizontal line beneath 0
    set(gca,'XColor',[1 1 1]);
end
 
end