
% Function to compute tSNR map from fMRI data.

function tsnrMap(inputPath,outputPath,outputDescription)

% Define i/o directories/names
[inputDir,inputName,~] = fpp.util.fileParts(inputPath);
if isempty(inputDir), inputDir = pwd; end
if ~exist('outputPath','var')
    outputPath = [inputDir '/' strrep(inputName,'_bold','') '_tsnr.nii.gz'];
end
[outputDir,~,~] = fpp.util.fileParts(outputPath);
if isempty(outputDir), outputDir = pwd; end

% Compute tSNR
fpp.util.system(['fslmaths ' inputPath ' -Tmean ' inputDir '/input_mean2350982452.nii.gz']);
fpp.util.system(['fslmaths ' inputPath ' -Tstd ' inputDir '/input_std2350982452.nii.gz']);
fpp.util.system(['fslmaths ' inputDir '/input_mean2350982452.nii.gz -div ' inputDir ...
    '/input_std2350982452.nii.gz ' outputPath]);
fpp.util.system(['rm -rf ' inputDir '/input_mean2350982452.nii.gz ' inputDir '/input_std2350982452.nii.gz']);

% Generate tSNR map json file
if ~isempty(fpp.bids.getMetadata(inputPath))
    fpp.bids.jsonReconstruct(inputPath,outputPath,'mri');
    if exist('outputDescription','var') && ischar(outputDescription)
        fpp.bids.jsonChangeValue(outputPath,'Description',outputDesciption);
    else
        fpp.bids.jsonChangeValue(outputPath,'Description','- tSNR map',1);
    end
    fpp.bids.jsonChangeValue(outputPath,'Sources',fpp.bids.removeBidsDir(inputPath));
end

end