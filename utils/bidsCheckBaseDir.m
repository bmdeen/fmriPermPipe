
% Function to decompose a BIDS-style path into a BIDS base directory,
% and everything after the BIDS base directory.

function [bidsBaseDir,bidsInternalDirs] = bidsCheckBaseDir(inputPath)

% Check for BIDS directory structure
if strcmp(inputPath(end),'/'), inputPath = inputPath(1:end-1); end
inputPathTerms = split(inputPath,'/');
for i=0:length(inputPathTerms)-1
    if i<length(inputPathTerms)-1 && ~isempty(regexp(inputPathTerms{end-i},'ses-.*'))...
            && ~isempty(regexp(inputPathTerms{end-i-1},'sub-.*'))
        % Input directory contains a BIDS session directory
        disp(['hi' int2str(i)]);
        bidsBaseDir = join(inputPathTerms(1:end-i-2),'/');
        bidsInternalDirs = join(inputPathTerms(end-i-1:end),'/');
        return;
    elseif ~isempty(regexp(inputPathTerms{end-i},'sub-.*'))
        % Input directory contains a BIDS subject directory, no session dir
        disp(i);
        bidsBaseDir = join(inputPathTerms(1:end-i-1),'/');
        bidsInternalDirs = join(inputPathTerms(end-i:end),'/');
        return;
    end
end
% Not BIDS formatted
bidsBaseDir = '';
bidsInternalDirs = inputPath;

end