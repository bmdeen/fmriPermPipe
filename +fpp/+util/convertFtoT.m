
% Function to take Freesurfer-style F.nii and sig.nii volumes, and generate
% a t-statistic map. NOTE: F-statistic must have come from a single linear
% contrast.
%
% Arguments: full paths to F.nii and sig.nii volumes.

function convertFtoT(fMapPath,sigMapPath,outputPath)

if ~exist(fMapPath,'file') || ~exist(sigMapPath,'file')
    fprintf('%s\n','ERROR: input images could not be found.');
    return;
end

[mapDir,~,~] = fileparts(fMapPath);
fileRoot = [mapDir '/tmp240592460928512089'];

extensions = {'.mgz','.MGZ','.mgh','.MGH'};
for e=1:length(extensions)
    if contains(fMapPath,extensions{e})
        fpp.util.system(['mri_convert ' fMapPath ' ' strrep(fMapPath,extensions{e},'.nii.gz')]);
        fMapPath = strrep(fMapPath,extensions{e},'.nii.gz');
    end
    if contains(sigMapPath,extensions{e})
        fpp.util.system(['mri_convert ' sigMapPath ' ' strrep(sigMapPath,extensions{e},'.nii.gz')]);
        sigMapPath = strrep(sigMapPath,extensions{e},'.nii.gz');
    end
end

fpp.util.system(['fslmaths ' sigMapPath ' -div 1000 -bin ' fileRoot '-posVals.nii.gz']);
fpp.util.system(['fslmaths ' sigMapPath ' -div 1000 -mul -1 -bin ' fileRoot '-negVals.nii.gz']);
fpp.util.system(['fslmaths ' fMapPath ' -sqrt -mul ' fileRoot '-posVals.nii.gz ' fileRoot '-posTVals.nii.gz']);
fpp.util.system(['fslmaths ' fMapPath ' -sqrt -mul ' fileRoot '-negVals.nii.gz ' fileRoot '-negTVals.nii.gz']);
fpp.util.system(['fslmaths ' fileRoot '-posTVals.nii.gz -sub ' fileRoot '-negTVals.nii.gz ' outputPath]);
fpp.util.system(['rm -rf ' fileRoot '-posVals.nii.gz ' fileRoot '-negVals.nii.gz ' ...
    fileRoot '-posTVals.nii.gz ' fileRoot '-negTVals.nii.gz']);

end