function meta = bidsGetMetadata(filename, pattern)

% Wrapper for bids-matlab function get_metadata

currentDir = pwd;
bidsMatlabDir = fileparts(which('bids.query'));
cd([bidsMatlabDir '/private']);
if nargin == 1
    meta = get_metadata(filename);
else
    meta = get_metadata(filename,pattern);
end
cd(currentDir);

return