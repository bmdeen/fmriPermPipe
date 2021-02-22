
% Wrapper for bids-matlab function get_metadata
%
% meta = fpp.bids.getMetadata(filename)

function meta = getMetadata(filename)

% Replace initial tilda with home directory, otherwise get_metadata will
% fail for tilda input
if strcmp(filename(1:2),'~/')
    filename = [getenv('HOME') filename(2:end)];
end

% If file is specified in current dir, add directory string
[fPath,~,~] = fpp.util.fileParts(filename);
if isempty(fPath)
    filename = [pwd '/' filename];
end

% Move to bids-matlab private dir, run get_metadata, move back
currentDir = pwd;
bidsMatlabDir = fileparts(which('bids.query'));
if isempty(bidsMatlabDir)
    error('bids-matlab is not sourced.');
end
cd([bidsMatlabDir '/private']);
meta = get_metadata(filename);
cd(currentDir);

return