
% Function to check the value specified for a given key in a BIDS filename.
%
% For instance, fpp.bids.checkNameValue('sub-01_desc-stuff.nii.gz','desc')
% returns 'stuff'
%
% value = fpp.bids.checkNameValue(inputPath,key)
%
% Arguments:
% - inputPath (string): path or filename to evaluate
% - key (string): key to check value of

function value = checkNameValue(inputPath,key)

value = '';

[~,inputName,~] = fpp.util.fileParts(inputPath);
keyLength = length(key);

if strcmp(inputName(1:keyLength),key)
    [startInd,finishInd] = regexp(inputName,[key '-[a-zA-Z0-9]+_']);
    if sum(startInd>0)
        value = inputName(startInd(1)+keyLength+1:finishInd(1)-1);
    end
else
    [startInd,finishInd] = regexp(inputName,['_' key '-[a-zA-Z0-9]+_']);
    if sum(startInd>0)
        value = inputName(startInd(1)+keyLength+2:finishInd(1)-1);
    end
end

end