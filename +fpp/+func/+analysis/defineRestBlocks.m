
% Function to define 3-column rest block regressor, based on 3-column task
% regressor input.
%
% restRegr = fpp.func.analysis.defineRestBlocks(taskRegr,exptDuration,regrWeight)

function restRegr = defineRestBlocks(taskRegr,exptDuration,regrWeight)

restRegr = [];
[~,indSort] = sort(taskRegr(:,1));
taskRegr = taskRegr(indSort,:);
stimEnd = taskRegr(end,1)+taskRegr(end,2);  % End of experimental stims
if size(taskRegr,1)==1
    if taskRegr(1)>0
        restRegr = [0 taskRegr(1,1) regrWeight];
    end
else
    blockEnd = taskRegr(:,1)+taskRegr(:,2); % Block end times
    blockBeforeRestInd = find(taskRegr(2:end,1)~=blockEnd(1:end-1));  % Indices of blocks preceding rest periods
    if ~isempty(blockBeforeRestInd)
        restRegr = [0 taskRegr(1,1) regrWeight; [blockEnd(blockBeforeRestInd) taskRegr(blockBeforeRestInd+1,1)-...
            blockEnd(blockBeforeRestInd) regrWeight*ones(length(blockBeforeRestInd),1)]];
    else
        restRegr = [0 taskRegr(1,1) regrWeight];
    end
    restRegr = restRegr(restRegr(:,2)>0,:);
end
if exptDuration>stimEnd
    restRegr = [restRegr; stimEnd exptDuration-stimEnd regrWeight];
end

end