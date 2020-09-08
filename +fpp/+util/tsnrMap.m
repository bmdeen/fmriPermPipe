
% Function to compute tSNR map from fMRI data.
%
% TODO: Add json file

function tsnrMap(inputPath,outputPath)

[inputDir,inputName,~] = fpp.util.fileParts(inputPath);
if isempty(inputDir), inputDir = pwd; end
if ~exist('outputPath','var')
    outputPath = [inputDir '/' strrep(inputName,'_bold','') '_tsnr.nii.gz'];
end

system(['fslmaths ' inputPath ' -Tmean ' inputDir '/input_mean2350982452.nii.gz']);
system(['fslmaths ' inputPath ' -Tstd ' inputDir '/input_std2350982452.nii.gz']);
system(['fslmaths ' inputDir '/input_mean2350982452.nii.gz -div ' inputDir ...
    '/input_std2350982452.nii.gz ' outputPath]);
system(['rm -rf ' inputDir '/input_mean2350982452.nii.gz ' inputDir '/input_std2350982452.nii.gz']);

end