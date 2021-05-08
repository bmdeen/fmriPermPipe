
% Function to identify good/bad voxels in fMRMI data. Remove voxels with a 
% temporal coefficient of variance greater than 0.5 standard deviations of 
% their local neighborhood.
% 
% fpp.wb.defineGoodVoxels(inputPath,gmPath,goodVoxPath,badVoxPath)
%
% Arguments:
% - inputPath (string): path to input fMRI data
% - gmPath (string): path to cortical gray matter (ribbon) mask image
% - goodVoxPath (string): output good voxel mask path
% - badVoxPath (string, optional): output bad voxel mask path

function defineGoodVoxels(inputPath,gmPath,goodVoxPath,badVoxPath)

% Constants
neighborhoodSmoothing = 5;  % Gaussian sigma defining extent of local neighborhood
sdFactor = .5;              % Remove voxels with CoV this many SDs above local mean

[outputDir,~,~] = fpp.util.fileParts(goodVoxPath);
if isempty(outputDir), outputDir = pwd; end
[~,inputName,inputExt] = fpp.util.fileParts(inputPath);
outputStem = [outputDir '/' fpp.bids.changeName(inputName,'desc',...
    'tmpDefineGoodVoxels12093518033','','')];

if strcmp(inputExt,'.nii.gz')
    fpp.util.system(['cp ' inputPath ' ' outputStem 'input.nii.gz']);
else
    fpp.util.system(['mri_convert ' inputPath ' ' outputStem 'input.nii.gz']);
end

% Compute coefficient of variation
fpp.fsl.maths([outputStem 'input'],'-Tmean',[outputStem 'mean']);
fpp.fsl.maths([outputStem 'input'],'-Tstd',[outputStem 'std']);
fpp.fsl.maths([outputStem 'std'],['-div ' outputStem 'mean'],[outputStem 'cov']);

% Mask CoV by cortical gm
fpp.fsl.maths([outputStem 'std'],['-mas ' gmPath],[outputStem 'cov_ribbon']);

% Normalize CoV (divide by mean)
[~,result] = fpp.util.system(['fslstats ' outputStem 'cov_ribbon -M']);
covMean = str2num(strtrim(result));
fpp.fsl.maths([outputStem 'cov_ribbon'],['-div ' num2str(covMean)],[outputStem 'cov_ribbon_norm']);

% Smooth normalized CoV within ribbon
fpp.fsl.maths([outputStem 'cov_ribbon_norm'],['-bin -s ' num2str(neighborhoodSmoothing)],[outputStem 'SmoothNorm']);
fpp.fsl.maths([outputStem 'cov_ribbon_norm'],['-s '  num2str(neighborhoodSmoothing) ' -div '...
    outputStem 'SmoothNorm -dilD'],[outputStem 'cov_ribbon_norm_s' num2str(neighborhoodSmoothing)]);

% Define CoV modulation: CoV unsmoothed / CoV smoothed.
fpp.fsl.maths([outputStem 'cov'],['-div ' num2str(covMean) ' -div ' outputStem...
    'cov_ribbon_norm_s' num2str(neighborhoodSmoothing)],[outputStem 'cov_norm_modulate']);
fpp.fsl.maths([outputStem 'cov_norm_modulate'],['-mas ' gmPath],[outputStem 'cov_norm_modulate_ribbon']);

% Compute mean and standard deviation of CoV modulation
[~,result] = fpp.util.system(['fslstats ' outputStem 'cov_norm_modulate_ribbon -S']);
stdVal = str2num(strtrim(result));
[~,result] = fpp.util.system(['fslstats ' outputStem 'cov_norm_modulate_ribbon -M']);
meanVal = str2num(strtrim(result));
covThresh = meanVal+stdVal*sdFactor;

% Define good and bad voxels based on CoV modulation
fpp.fsl.maths([outputStem 'mean'],'-bin',[outputStem 'mask']);
fpp.fsl.maths([outputStem 'cov_norm_modulate'],['-thr ' num2str(covThresh) ' -bin -sub '...
    outputStem 'mask -mul -1'],goodVoxelPath);
if exist('badVoxelPath','var') && ~isempty(badVoxelPath)
    fpp.fsl.maths([outputStem 'mask'],['-sub ' goodVoxelPath],badVoxelPath);
end

%fpp.util.system(['rm -rf ' outputStem '*']);

end