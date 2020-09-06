
% Script to copy MRI image file (functional, structural, or field map), as
% well as JSON sidecar (if it exists), maintaining fields used by
% fmriPermPipe, based on BIDS 1.4.0 specification.
%
% Arguments:
% - inputImagePath (string)
% - outputImagePath (string)
%
% Dependencies: bids-matlab (required), bids-matlab-tools (recommended for
% JSONio)

function copyImageAndJson(inputImagePath,outputImagePath)


[~,~,inputExt] = fpp.util.fileParts(inputImagePath);
[~,~,outputExt] = fpp.util.fileParts(outputImagePath);

inputJsonPath = strrep(inputImagePath,inputExt,'.json');
outputJsonPath = strrep(outputImagePath,outputExt,'.json');

if strcmpi(inputExt,outputExt)
    system(['cp ' inputImagePath ' ' outputImagePath]);
else
    system(['mri_convert ' inputImagePath ' ' outputImagePath]);
end

fpp.bids.jsonReconstruct(inputJsonPath,outputJsonPath);


end