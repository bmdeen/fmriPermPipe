
% Function to generate carpet plot: plot of time series from all voxels in
% an image, separated by image segment, and measures of data quality.
%
% fpp.util.carpetPlot(inputPath,maskPath,segmentMaskPaths,nuisanceSeries,...
%   nuisanceNames[,segmentColors,nuisanceColors])
%
% Arguments:
% - inputPath (string): path to input functional data
% - maskPath (string): path to input brain mask. Should exclude any voxels
%   with zero signal in functional data
% - segmentMaskPaths (cell array of strings): paths to segment masks.
%   Segments should be non-overlapping. First mask is assumed to be GM or 
%   cortical GM, for the sake of computing % signal change. Default color
%   scheme is intended for: {GM, WM, CSF}
% - nuisanceSeries (numeric matrix): D x T matrix of D nuisance time series
%   to plot, of length T. Default color scheme is intended for: {Global
%   mean, WM mean, CSF mean, Framewise Displacement, DVARS_std}
% - nuisanceNames (cell array of strings): name of each nuisance signal
% - outputPath (string): path to output image
%
% Optional arguments:
% - segmentColors (cell array of color vectors): colors to label each
%   segment
% - nuisanceColors (cell array of color vectors): colors to plot each
%   nuisance time series
%
% For further information about the plot, see
% https://www.jonathanpower.net/2017-ni-the-plot.html


function carpetPlot(inputPath,maskPath,segmentMaskPaths,nuisanceSeries,...
    nuisanceNames,outputPath,segmentColors,nuisanceColors)

% Define parameters
hLims = [-5 5];   % Grayscale heatmap limits, in % signal change relative to GM mean
titleFontSize = 14;
labelFontSize = 14;
tickFontSize = 10;
figSizeDefault = 1200;  % Figure height/widght in pixels

if ~exist('segmentColors','var') || isempty(segmentColors)
    segmentColors = {[0 .5 1],[0 1 0],[1 1 0]};
end
nSegs = length(segmentMaskPaths);
while length(segmentColors)<nSegs
    segmentColors = [segmentColors segmentColors];
end

if ~exist('nuisanceColors','var') || isempty(nuisanceColors)
    nuisanceColors = {[61 165 193]/255,[45 167 111]/255,[164 146 43]/255,...
        [246 102 126]/255,[197 111 242]/255};
end
if size(nuisanceSeries,1) > size(nuisanceSeries,2)
    nuisanceSeries = nuisanceSeries';
end
nNuis = size(nuisanceSeries,1);
while length(nuisanceColors)<nNuis
    nuisanceColors = [nuisanceColors nuisanceColors];
end
nSubPlotHeatMap = 5;    % How many subplots the heat map should occupy
nSubPlot = nNuis + nSubPlotHeatMap;   % One subplot for each nuisance signal, nSubPlotHeatMap for heat map

% Load brain mask and functional data
maskData = fpp.util.mriRead(maskPath);
dataMat = fpp.util.readDataMatrix(inputPath,maskData.vol);
vols = size(dataMat,2); xLims = [1 vols];
if vols~=size(nuisanceSeries,2)
    error('Size of functional data and nuisance time series do not match.');
end

% Load image segment masks, within brain mask
segLengths = [];
for s=1:nSegs
    segMat{s} = fpp.util.readDataMatrix(segmentMaskPaths{s},maskData.vol);
    segLengths(s) = sum(segMat{s}==1);
end
segLengthsCumulative = [0 cumsum(segLengths)];

% Demean data, convert to % signal change relative to gray matter mean
meanMat = repmat(mean(dataMat,2),[1 vols]);
gmMean = mean(mean(dataMat(segMat{1}==1,:)));   % Assuming first segment is GM or cortical GM
dataMat = (dataMat-meanMat)/gmMean*100;

% Arrange data by image segment
dataMatOrdered = [];
for s=1:nSegs
    dataMatOrdered = [dataMatOrdered; dataMat(segMat{s}==1,:)];
end

% Check screen size, reduce figure size if necessary
set(0,'units','pixels')
screenRect = get(0,'screensize');
figSize = min([figSizeDefault screenRect(3:4)]);

% Initiate plot
h=figure;
set(h,'position',[0 0 figSize figSize]);
set(h,'visible','on');
set(h,'Color',[1 1 1]);     % white bg

% Plot nuisance time series
for n=1:nNuis
    subplot(nSubPlot,1,n);
    plot(nuisanceSeries(n,:),'Color',nuisanceColors{n},'LineWidth',2);
    set(gca,'FontSize',tickFontSize);
    xlim(xLims);
    title(nuisanceNames{n},'FontSize',titleFontSize,'Interpreter','none');
    box off;
end

% Plot heat map
subplot(nSubPlot,1,nSubPlot-nSubPlotHeatMap+1:nSubPlot);
imagesc(dataMatOrdered,hLims);
set(gca,'ytick',[],'FontSize',tickFontSize);
colormap(gray);
xlim(xLims);
ylabel('Voxels','FontSize',labelFontSize);
xlabel('Volume #','FontSize',labelFontSize);

% Draw colored rectangles on heat map, defining each segment
xStart = 0; xEnd = ceil(vols/100); f = 1:4;
for s=1:nSegs
    yStart = segLengthsCumulative(s)+.5; yEnd = segLengthsCumulative(s+1)+.5;
    v = [xStart yStart; xStart yEnd; xEnd yEnd; xEnd yStart];
    patch('Faces',f,'Vertices',v,'FaceColor',segmentColors{s},'LineStyle','none');
end

% Save and close plot
saveas(h,outputPath);
close(h);

end