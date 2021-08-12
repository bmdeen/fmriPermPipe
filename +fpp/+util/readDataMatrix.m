
% Wrapper function to read data from a NIFTI 3D/4D volume or CIFTI dscalar
% or dtseries file, and convert it to a coordinate by 1 vector or
% coordinate by time matrix.
%
% [dataMat,hdr] = fpp.util.readDataMatrix(inputPath,maskData)
%
% Arguments:
% - inputPath (string): path to NIFTI/CIFTI input file
%
% Optional arguments:
% - maskData (vector or 3D matrix): mask data, either 3D volume (NIFTI)
%       or 1D vector (CIFTI) of 1s and 0s. dataMat will only contain values
%       where maskData==1.
%
% Dependencies: Freesurfer matlab toolbox, CIFTI matlab toolbox

function [dataMat,hdr] = readDataMatrix(inputPath,maskData)

[~,~,inputExt] = fpp.util.fileParts(inputPath);
if ismember(lower(inputExt),{'.nii.gz','.nii'})
    isCifti = 0;
elseif ismember(lower(inputExt),{'.dlabel.nii','.dscalar.nii','.dtseries.nii'})
    isCifti = 1;
else
    error('fpp.util.readDataMatrix is intended for NIFTI or dlabel/dscalar/dtseries CIFTI files.');
end

% Check maskData input
useMask = 0;
if exist('maskData','var') && ~isempty(maskData)
    useMask = 1;
    if isCifti && ~isvector(maskData)
        error('For CIFTI input, maskData must be a vector.');
    elseif ~isCifti && numel(size(maskData))~=3
        error('For NIFTI input, maskData must be a 3D matrix.');
    end
end

if isCifti
    hdr = cifti_read(inputPath);
    dataMat = hdr.cdata;
    hdr.cdata = [];
    if useMask
        dataMat = dataMat(maskData==1,:);
    end
else
    hdr = fpp.util.mriRead(inputPath);
    dims = size(hdr.vol);
    if length(dims)<=3
        dataMat = hdr.vol(:);
    else
        dataMat = reshape(hdr.vol,[prod(dims(1:3)) dims(4)]);
    end
    if useMask
        dataMat = dataMat(maskData(:)==1,:);
    end
    hdr.vol = [];
end

end