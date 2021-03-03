
% Wrapper for bids-matlab function get_metadata
%
% meta = fpp.bids.getMetadata(filename)

function meta = getMetadata(filename)

meta = struct();

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

% Don't check metadata for BIDS-invalid filenames
if ~fpp.bids.validFilename(filename), return; end

meta = bids.internal.get_metadata(filename);

return