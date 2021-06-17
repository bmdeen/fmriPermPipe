
% Function to convert BIDS-format events.tsv file to FSL-format 3-column
% regressor files.
%
% fpp.util.convertEventsTSV2FSL(inputTSVPath,outputStem)

function convertEventsTSV2FSL(inputTSVPath,outputStem)

tsvData = bids.util.tsvread(inputTSVPath);

conds = unique(tsvData.trial_type);

for c=1:length(conds)
    conditionInd = find(ismember(tsvData.trial_type,conds{c}));
    regr3Col = [tsvData.onset(conditionInd) tsvData.duration(conditionInd) ones(length(conditionInd),1)];
    outPath = [outputStem '_desc-' conds{c} '_regr3col'];
    
    fid = fopen(outPath,'w+');
    fprintf(fid,'%6.2f\t%6.2f\t%6.2f\n',regr3Col');
    fclose(fid);
end


end