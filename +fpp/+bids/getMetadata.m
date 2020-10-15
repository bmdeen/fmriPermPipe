function meta = getMetadata(filename, pattern)

% Wrapper for bids-matlab function get_metadata

% Replace initial tilda with home directory, otherwise get_metadata will
% fail for tilda input
if strcmp(filename(1:2),'~/')
    filename = [getenv('HOME') filename(2:end)];
end

currentDir = pwd;
bidsMatlabDir = fileparts(which('bids.query'));
if isempty(bidsMatlabDir)
    error('Error: bids-matlab is not sourced.');
end
cd([bidsMatlabDir '/private']);
if nargin == 1
    meta = get_metadata(filename);
else
    meta = get_metadata(filename,pattern);
end
cd(currentDir);

return