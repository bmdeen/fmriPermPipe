
% Wrapper for MATLAB's fileparts
% - Accommodates .gz extension, e.g. in .nii.gz
% - Accomodates GIFTI/CIFTI extensions, e.g. .surf.gii or .dscalar.nii

function [filepath,name,ext] = fileParts(filename)

if length(filename)>=3 && strcmpi(filename(end-2:end),'.gz')
    extGZ = filename(end-2:end);
    filename = filename(1:end-3);
    [filepath,name,ext] = fileparts(filename);
    ext = [ext extGZ];
else
    [filepath,name,ext] = fileparts(filename);
end

giiCiiLabels = {'surf','shape','label','dconn','dscalar','dtseries','dlabel',...
    'dpconn','pconn','pdconn','pscalar','ptseries','plabel','sdseries','fiberTemp'};

for i=1:length(giiCiiLabels)
    if contains(name,['.' giiCiiLabels{i}])
        strrep(name,['.' giiCiiLabels{i}],'');
        ext = ['.' giiCiiLabels{i} ext];
        break;
    end
end

end