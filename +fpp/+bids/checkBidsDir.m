
% Function to decompose a BIDS-style path into a BIDS base directory -
% defined as one level above the sub directory - and everything after
% the base directory.
%
% [bidsBaseDir,bidsInternalDirs] = fpp.bids.checkBidsDir(inputPath)

function [bidsBaseDir,bidsInternalDirs] = checkBidsDir(inputPath)

% Check for BIDS directory structure
if strcmp(inputPath(end),'/'), inputPath = inputPath(1:end-1); end
inputPathTerms = split(inputPath,'/');
if contains(inputPathTerms{end},'.') || exist(inputPath)==2
    startInd = 1;   % If last term in path is a file, don't consider it a potential directory
else
    startInd = 0;
end
for i=startInd:length(inputPathTerms)-1
    if i<length(inputPathTerms)-1 && ~isempty(regexp(inputPathTerms{end-i},'ses-.*'))...
            && ~isempty(regexp(inputPathTerms{end-i-1},'sub-.*'))
        % Input directory contains a BIDS session directory
        bidsBaseDir = join(inputPathTerms(1:end-i-2),'/');
        bidsBaseDir = bidsBaseDir{1};
        bidsInternalDirs = join(inputPathTerms(end-i-1:end),'/');
        bidsInternalDirs = bidsInternalDirs{1};
        return;
    elseif ~isempty(regexp(inputPathTerms{end-i},'sub-.*'))
        % Input directory contains a BIDS subject directory, no session dir
        bidsBaseDir = join(inputPathTerms(1:end-i-1),'/');
        bidsBaseDir = bidsBaseDir{1};
        bidsInternalDirs = join(inputPathTerms(end-i:end),'/');
        bidsInternalDirs = bidsInternalDirs{1};
        return;
    end
end

% Not BIDS formatted
bidsBaseDir = '';
bidsInternalDirs = inputPath;

end