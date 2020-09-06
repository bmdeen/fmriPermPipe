
% Wrapper for MATLAB's fileparts, deal with .gz extension, e.g. in .nii.gz

function [filepath,name,ext] = fileParts(filename)

if strcmpi(filename(end-2:end),'.gz')
    extGZ = filename(end-2:end);
    filename = filename(1:end-3);
    [filepath,name,ext] = fileparts(filename);
    ext = [ext extGZ];
else
    [filepath,name,ext] = fileparts(filename);
end


end