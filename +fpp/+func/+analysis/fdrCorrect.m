
% Function to FDR-correct z-statistic image, output thresholded image
% 
% [critZ,critP] = fpp.func.analysis.fdrCorrect(inputPath,outputPath[,maskPath,qThresh,tails])
% 
% Arguments:
%   - inputPath (string): path to input z-statistic NIFTI/CIFTI file
%   - outputPath (string): path to output (corrected) NIFTI/CIFTI file
% 
% Optional arguments:
%   - maskPath (string optional): path to 3D mask volume (NIFTI only)
%   - qThresh (scalar): FDR q-threshold
%   - tails (scalar): 1 or 2 tailed test
% 
% Outputs:
%   - critZ: z-value cutoff
%   - critP: p-value cutoff
% 
% Note: script will treat Z=0 values as not tested.

function [critZ,critP] = fdrCorrect(inputPath,outputPath,maskPath,qThresh,tails)

method = 'pdep';    % pdep = Benjamini & Hochberg, dep = Benjamini & Yekutieli
[~,~,inputExt] = fpp.util.fileParts(inputPath);
isCifti = 0;
if strcmpi(inputExt,'.dscalar.nii')
    isCifti = 1;
end

if ~exist('maskPath','var') || isempty(maskPath)
    maskPath = '';
    maskData = [];
end
if ~exist('qThresh','var') || ~(isnumeric(qThresh) && isscalar(qThresh) && ...
        qThresh<1 && qThresh>0)
    qThresh = .05;
end

if ~exist('tails','var') || ~(isnumeric(tails) && isscalar(tails) && ismember(tails,[1 2]))
    tails = 2;
end

% Load mask
if ~isempty(maskPath)
    if isCifti
        maskData = fpp.util.readDataMatrix(maskPath);
    else
        maskImage = fpp.util.mriRead(maskPath);
        maskData = maskImage.vol;
    end
end

% Load data
[zVec,hdr] = fpp.util.readDataMatrix(inputPath,maskData);

% Convert z-values to p-values
pVec = zeros(size(zVec));
if tails==1
    pVec = normcdf(-zVec);
else
    pVec(zVec>0) = 2*normcdf(-zVec(zVec>0));
    pVec(zVec<0) = 2*normcdf(zVec(zVec<0));
end

% Compute FDR threshold
[h1,critP,~,~] = fpp.func.analysis.fdrBH(pVec(zVec~=0),qThresh,method);
h = zeros(size(zVec));
h(zVec~=0) = h1;
critZ = abs(icdf('norm',critP,0,1));
zVec(h==0) = 0;

% Write output thresholded image
fpp.util.writeDataMatrix(zVec,hdr,outputPath,maskData);

end