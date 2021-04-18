
% Function to read a set of volume indices from a BIDS-style outliers.tsv
% file, with delta function regressors at each outlier time point.
%
% [badVols,nVols] = fpp.util.readOutlierTSV(outlierPath)

function [badVols,nVols] = readOutlierTSV(outlierPath)

badVols = [];

tsv = bids.util.tsvread(outlierPath);
fields = fieldnames(tsv);
for f=1:length(fields)
    thisData = tsv.(fields{f});
    if isnumeric(thisData) && isvector(thisData) && sum(thisData==1)==1 % If field has delta function value
        badVols = [badVols find(thisData==1)];
        nVols = length(thisData);
    end
end

end