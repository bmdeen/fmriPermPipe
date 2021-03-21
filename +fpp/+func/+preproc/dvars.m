
% Function to compute DVARS for an fMRI dataset - IE, the root-mean-square
% (over space) of the temporal difference image, as a function of time.
%
% [output,outputStd] = dvars(inputPath,maskPath)
%
% Arguments:
% - inputPath (string): path to input functional data
% - maskPath (string): path to brain mask
%
% Standardized DVARS computed based on Nichols (2013), "Notes on Creating a 
% Standardized Version of DVARS"

function [output,outputStd] = dvars(inputPath,maskPath)

[inputDir,inputName,~] = fpp.util.fileParts(inputPath);
if isempty(inputDir), inputDir = pwd; end
tmpDir = [inputDir '/' inputName '_dvars1305981354134'];
mkdir(tmpDir);

% Define intermediate image paths
tmeanPath = [tmpDir '/' fpp.bids.changeName(inputName,'desc','tmean') '.nii.gz'];
lqPath = [tmpDir '/' fpp.bids.changeName(inputName,'desc','lowerquartile') '.nii.gz'];
uqPath = [tmpDir '/' fpp.bids.changeName(inputName,'desc','upperquartile') '.nii.gz'];
sdPath = [tmpDir '/' fpp.bids.changeName(inputName,'desc','stddev') '.nii.gz'];
tdiffSDPath = [tmpDir '/' fpp.bids.changeName(inputName,'desc','tdiffStdDevEstimate') '.nii.gz'];
ar1Path = [tmpDir '/' fpp.bids.changeName(inputName,'desc','ar1') '.nii.gz'];

% Compute temporal mean
fpp.fsl.maths(inputPath,'-Tmean',tmeanPath);

% Compute robust standard deviation estimate
fpp.fsl.maths(inputPath,'-Tperc 25',lqPath);
fpp.fsl.maths(inputPath,'-Tperc 75',uqPath);
fpp.fsl.maths(uqPath,['-sub ' lqPath ' -div 1.349'],sdPath);

% Compute (non-robust) estimate of lag-1 autocorrelation
fpp.fsl.maths(inputPath,['-sub ' tmeanPath ' -Tar1'],ar1Path);

% Compute (predicted) standard deviation of temporal difference time series
fpp.fsl.maths(ar1Path,['-mul -1 -add 1 -mul 2 -sqrt -mul ' sdPath],tdiffSDPath);

% Load brain mask, data, and temporal difference standard deviation
mask = fpp.util.mriRead(maskPath);
dataMat = fpp.util.readDataMatrix(inputPath,mask.vol);
tdiffSDMat = fpp.util.readDataMatrix(tdiffSDPath,mask.vol);

% Compute DVARS
dataDiff = zeros(size(dataMat));
dataDiff(:,2:end) = diff(dataMat,1,2);
output = rms(dataDiff)';

% Compute standardized DVARS
outputStd = output/mean(tdiffSDMat);

fpp.util.system(['rm -rf ' tmpDir]);

end