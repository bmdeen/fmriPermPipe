
% Function to reorient a 3D/4D MRI dataset to LAS orientation (radiological
% convention).
%
% CAUTION: SHOULD ONLY BE USED WHEN YOU'RE CONFIDENT THAT AN IMAGE HAS
% ACCURATE ORIENTATION INFO IN NIFTI HEADER.

function reorientToLAS(inputPath,outputPath)

warning(sprintf('%s\n',['CAUTION: reorientToLAS should only be used when you are confident '...
    'that the image''s NIFTI header contains accurate orientation info.']));

[~,~,inputExt] = fpp.util.fileparts(inputPath);
[outputDir,~,outputExt] = fpp.util.fileparts(outputPath);
if isempty(outputDir), outputDir = pwd; end
tmpDir = [outputDir '/tmp_reorientToLASIntermediateDir'];
mkdir(tmpDir);

tmpPath = [tmpDir '/newImage.nii.gz'];
tmpPath2 = [tmpDir '/newImageSwap.nii.gz'];
if strcmp(inputExt,'.nii.gz')
    fpp.util.system(['cp ' inputPath ' ' tmpPath]);
else
    fpp.util.system(['mri_convert ' inputPath ' ' tmpPath]);
end

inputOrientation = fpp.util.getImageOrientation(tmpPath);
if ~strcmp(inputOrientation,'LAS')
    
    % If image isn't in LAS/RAS orientation, rotate to this orientation
    if ~strcmp(inputOrientation,'RAS')
        fpp.util.system(['fslreorient2std ' tmpPath ' ' tmpPath]);
    end
    
    % If fslreorient2std converted image to RAS, a left-right reflection is
    % required to get to LAS.
    imageOrientation = fpp.util.getImageOrientation(tmpPath);
    if ~strcmp(imageOrientation,'LAS')
        fpp.util.system(['cp ' tmpPath ' ' tmpPath2]);
        fpp.util.system(['fslorient -swaporient ' tmpPath2]);
        fpp.util.system(['echo $''-1 0 0 0\n0 1 0 0\n0 0 1 0\n0 0 0 1'' > ' tmpDir '/flipLR.mat']);
        fpp.util.system(['flirt -in ' tmpPath ' -ref ' tmpPath2 ' -applyxfm -init ' ...
            tmpDir '/flipLR.mat -out ' tmpPath]);
    end
end

if strcmp(outputExt,'.nii.gz')
    fpp.util.system(['cp ' tmpPath ' ' outputPath]);
else
    fpp.util.system(['mri_convert ' tmpPath ' ' outputPath]);
end
fpp.util.system(['rm -rf ' tmpDir]);

% Convert orientation of accompanying JSON file if it exist
inputJsonPath = strrep(inputPath,inputExt,'.json');
if exist(inputJsonPath,'file')
    outputJsonPath = strrep(outputPath,outputExt,'.json');
    fpp.util.system(['cp ' inputJsonPath ' ' outputJsonPath]);
    fpp.bids.jsonReorient(outputJsonPath,inputOrientation,'LAS')
end

end