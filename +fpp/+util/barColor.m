
% [b,s,e] = barColor(X,[colorArray,errorVec,gapInd])
%
% Function to plot a bar graph with specific bar colors, and overlaid
% circles for individual data points.
%
% Arguments: 
% - X (numeric 2D matrix): samples by condition array of values to plot
% 
% Optional arguments:
% - colorArray (cell array): array of colors by condition
% - errorVec (scalar or vector): whether to use error bars, or vector value
%       for error. 0 = no error bars, 1 = standard error, vector = use
%       specified values.
% - gapInd (vector): gaps in the graph will be placed after these bars
%
% Outputs:
% - b (Bar object): bar graph handle
% - s (Scatter object): scatter plot handle
% - e (ErrorBar object): error bar graph handle, if errors were used
% - barInd (vector): y-position of bars in graph

function [b,s,e,barInd] = barColor(X,colorArray,errorVec,gapInd)

e = [];
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

% Determine error bar values
errors = [];
if ~isscalar(errorVec)
    errors = zeros(1,nSlots);
    errors(barInd) = errorVec;
    errorVec = 1;
elseif errorVec==1
    errors = std(X,1)/sqrt(nSamples);
end

% Plot bar graph
b = bar(mean(X,1));
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
scatterX = []; scatterY = [];
for i=barInd
    scatterX = [scatterX; i*ones(nSamples,1)];
    scatterY = [scatterY; X(:,i)];
end
hold on;
s = scatter(scatterX,scatterY,50,'k','filled');

% Plot errors
if errorVec
    e = errorbar(mean(X,1),errors,'.k');
    set(e,'LineWidth',2);
end

% Set figure/axis properties
box off;
set(gcf,'Color',[1 1 1]);
set(gca,'LineWidth',2,'FontSize',30,'XTick',[]);
set(b,'LineWidth',2,'BarWidth',.75);
set(b.BaseLine,'LineWidth',2);
ylim = get(gca,'YLim');
if ylim(1)<0    % Remove horizontal line beneath 0
    set(gca,'XColor',[1 1 1]);
end
 
end