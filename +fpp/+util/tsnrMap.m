
% fpp.util.tsnrMap(inputPath[,outputPath,outputDescription])
%
% Function to compute tSNR map from fMRI data (NIFTI or CIFTI dtseries)

function tsnrMap(inputPath,outputPath,outputDescription)

% Define input file parts
[inputDir,inputName,inputExt] = fpp.util.fileParts(inputPath);
if isempty(inputDir), inputDir = pwd; end

% Check whether input is NIFTI or CIFTI
if ismember(lower(inputExt),{'.nii.gz','.nii'})
    isCifti = 0;
    imageType = 'volume';
    outputExt = '.nii.gz';
elseif ismember(lower(inputExt),{'.dtseries.nii'})
    isCifti = 1;
    imageType = 'cifti';
    outputExt = '.dscalar.nii';
else
    error('fpp.util.tsnrMap is intended for NIFTI or dtseries CIFTI files.');
end

% Define output path if it doesn't exist
if ~exist('outputPath','var') || isempty(outputPath)
    outputPath = fpp.bids.changeName(inputPath,[],[],'tsnr',outputExt);
end
[outputDir,~,~] = fpp.util.fileParts(outputPath);
if isempty(outputDir), outputDir = pwd; end

% Compute tSNR
if ~isCifti
    fpp.util.system(['fslmaths ' inputPath ' -Tmean ' inputDir '/' inputName '_mean2350982452.nii.gz']);
    fpp.util.system(['fslmaths ' inputPath ' -Tstd ' inputDir '/' inputName '_std2350982452.nii.gz']);
    fpp.util.system(['fslmaths ' inputDir '/' inputName '_mean2350982452.nii.gz -div ' inputDir ...
        '/' inputName '_std2350982452.nii.gz ' outputPath]);
    fpp.util.system(['rm -rf ' inputDir '/' inputName '_mean2350982452.nii.gz ' inputDir '/' inputName '_std2350982452.nii.gz']);
else
    [funcMat,hdr] = fpp.util.readDataMatrix(inputPath);
    tsnrVec = mean(funcMat,2)./std(funcMat,0,2);
    hdr.diminfo{2}.maps(1).name = 'tSNR map';
    hdr.diminfo{2}.maps(1).metadata.key = '';
    hdr.diminfo{2}.maps(1).metadata.value = '';
    fpp.util.writeDataMatrix(tsnrVec,hdr,outputPath);
end

% Generate tSNR map json file
if ~isempty(fpp.bids.getMetadata(inputPath))
    fpp.bids.jsonReconstruct(inputPath,outputPath,'mri');
    if exist('outputDescription','var') && ischar(outputDescription)
        fpp.bids.jsonChangeValue(outputPath,'Description',outputDescription);
    else
        fpp.bids.jsonChangeValue(outputPath,'Description',' - tSNR map',1);
    end
    fpp.bids.jsonChangeValue(outputPath,'Sources',fpp.bids.removeBidsDir(inputPath));
end

end