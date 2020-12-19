
% Script to copy MRI image file (functional, structural, or field map), as
% well as JSON sidecar (if it exists), maintaining fields used by
% fmriPermPipe, based on BIDS 1.4.1 specification.
%
% fpp.util.copyImageAndJson
%
% Arguments:
% - inputImagePath (string)
% - outputImagePath (string)
% - fieldsToKeep (optional, cell array of strings): JSON fields to copy
%    OR (string): label for image type, determining which fields to keep.
%    Options: midprepFMRI, fMRI, MRI, Mask, Seg, Surf, Cifti, keepAll
%
% Dependencies: bids-matlab (required), bids-matlab-tools (recommended for
% JSONio)

function copyImageAndJson(inputImagePath,outputImagePath,fieldsToKeep)

[~,~,inputExt] = fpp.util.fileParts(inputImagePath);
[~,~,outputExt] = fpp.util.fileParts(outputImagePath);

if strcmpi(inputExt,outputExt)
    fpp.util.system(['cp ' inputImagePath ' ' outputImagePath]);
else
    fpp.util.system(['mri_convert ' inputImagePath ' ' outputImagePath]);
end

if exist('fieldsToKeep','var') && ~isempty(fieldsToKeep)
    fpp.bids.jsonReconstruct(inputImagePath,outputImagePath,fieldsToKeep);
else
    fpp.bids.jsonReconstruct(inputImagePath,outputImagePath);
end


end