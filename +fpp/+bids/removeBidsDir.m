
% Function to remove BIDS base directory (defined as above sub directory)
% from path.

function outputPath = removeBidsDir(inputPath)

bidsBaseDir = fpp.bids.checkBidsDir(inputPath);  % BIDS derivative root directory
if isempty(bidsBaseDir), bidsBaseDir = '123059815310419841xyz'; end     % Hack so that strrep commands below work for non-BIDS input
outputPath = strrep(inputPath,[bidsBaseDir '/'],'');

end