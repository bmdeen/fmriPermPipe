
% Function to plot GLM regressor information: regressors, regressor
% correlations, and variance of task regressors explained by nuisance
% regressors.
%
% fpp.func.analysis.plotRegressors(taskRegrMat,nuisRegrMat,regrNames,outputDir,outputNameGeneric,closeFigs)

function plotRegressors(taskRegrMat,nuisRegrMat,regrNames,outputDir,outputNameGeneric,closeFigs)

titleSize = 8;     % Font size of plot titles

% Determine description field
descExp = '_desc-[0-9a-zA-Z]+_';
[descInd,descIndEnd] = regexp(outputNameGeneric,descExp);
if ~isempty(descInd)
    descValue = outputNameGeneric(descInd(1)+6:descIndEnd(1)-1);
else
    descValue = '';
end

X = [taskRegrMat nuisRegrMat];
nConds = size(taskRegrMat,2);

% 1) Plot all regressors
regrFig = figure('units','normalized','position',[.1 .1 .8 .8]);
if closeFigs, set(regrFig,'visible','off'); end
numRegrs = length(regrNames);
if numRegrs>15
    maxRows = 10;
else
    maxRows = 5;
end
figRows = min(numRegrs,maxRows);
figCols = ceil(numRegrs/maxRows);
for r=1:numRegrs
    subplot(figRows,figCols,r)
    plot(X(:,r));
    xlim([1 size(X,1)]);
    if r==1
        title({outputNameGeneric,regrNames{r}},'interpreter','none','FontSize',titleSize);
    else
        title(regrNames{r},'interpreter','none');
    end
end
set(gcf,'Color',[1 1 1]);
saveas(regrFig,[outputDir '/' fpp.bids.changeName(outputNameGeneric,'desc',...
    [descValue 'Regressors']) '_image.png']);
if closeFigs, close(regrFig); end

% 2) Plot regressor correlation matrix, scale [-1 1]
corrFig = figure;
if closeFigs, set(corrFig,'visible','off'); end
imagesc(corr(X),[-1 1]); colorbar;
set(gca,'XTick',1:numRegrs,'YTick',1:numRegrs,'XTickLabel',...
    [],'YTickLabel',regrNames,'TickLabelInterpreter','none');
set(gcf,'Color',[1 1 1]);
title({outputNameGeneric,'Regressor Correlations'},'interpreter','none','FontSize',titleSize);
saveas(corrFig,[outputDir '/' fpp.bids.changeName(outputNameGeneric,'desc',...
    [descValue 'RegressorCorrelations']) '_image.png']);
if closeFigs, close(corrFig); end

% 3) Plot proportion of task regressor variance removed by nuisance regressors
if ~isempty(nuisRegrMat)
    varianceFig = figure;
    if closeFigs, set(varianceFig,'visible','off'); end
    taskRegrMatProj = taskRegrMat-nuisRegrMat*inv(nuisRegrMat'*nuisRegrMat)*nuisRegrMat'*taskRegrMat;
    varRemoved = ones(1,nConds)-sum(taskRegrMatProj.^2)./sum(taskRegrMat.^2);
    bar(varRemoved);
    ylim([0 1]);
    set(gca,'XTick',1:nConds,'XTickLabel',regrNames(1:nConds),'TickLabelInterpreter','none');
    set(gcf,'Color',[1 1 1]);
    title({outputNameGeneric,'Variance Explained by Nuisance Regrs'},'interpreter','none','FontSize',titleSize);
	saveas(varianceFig,[outputDir '/' fpp.bids.changeName(outputNameGeneric,'desc',...
        [descValue 'RegressorVarianceRemoved']) '_image.png']);
    if closeFigs, close(varianceFig); end
end



end