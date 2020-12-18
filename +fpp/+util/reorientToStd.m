
% Function to reorient a 3D/4D MRI dataset to LAS/RAS orientation via a
% rotation. Also modifies relevant fields in JSON metadata.

function reorientToStd(inputPath,outputPath)

[~,~,inputExt] = fpp.util.fileParts(inputPath);
[outputDir,~,outputExt] = fpp.util.fileParts(outputPath);
if isempty(outputDir), outputDir = pwd; end
tmpDir = [outputDir '/tmp_reorientToStdIntermediateDir'];
mkdir(tmpDir);

tmpPath = [tmpDir '/newImage.nii.gz'];
if strcmp(inputExt,'.nii.gz')
    fpp.util.system(['cp ' inputPath ' ' tmpPath]);
else
    fpp.util.system(['mri_convert ' inputPath ' ' tmpPath]);
end

inputOrientation = fpp.util.getImageOrientation(tmpPath);
if sum(ismember(inputOrientation,{'LAS','RAS'}))==0
    
    % If image isn't in LAS/RAS orientation, rotate to this orientation
    if ~strcmp(inputOrientation,'RAS')
        fpp.util.system(['fslreorient2std ' tmpPath ' ' tmpPath]);
    end
    outputOrientation = fpp.util.getImageOrientation(tmpPath);
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
    if exist('outputOrientation','var') && ~strcmp(inputOrientation,outputOrientation)
        fpp.bids.jsonReorient(outputJsonPath,inputOrientation,outputOrientation);
    end
end

end