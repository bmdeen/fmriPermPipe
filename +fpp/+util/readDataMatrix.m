
% Wrapper function to read data from a NIFTI 3D/4D volume or CIFTI dscalar
% or dtseries file, and convert it to a coordinate by 1 vector or
% coordinate by time matrix.
%
% [dataMat,hdr] = fpp.util.readDataMatrix(inputPath,maskVol)
%
% Arguments:
%   - inputPath (string): path to NIFTI/CIFTI input file
%   - maskVol (3D matrix, optional): 3D mask volume (NIFTI only)
%
% Dependencies: Freesurfer matlab toolbox, CIFTI matlab toolbox

function [dataMat,hdr] = readDataMatrix(inputPath,maskVol)

[~,~,inputExt] = fpp.util.fileParts(inputPath);
if ismember(lower(inputExt),{'.nii.gz','.nii'})
    isCifti = 0;
elseif ismember(lower(inputExt),{'.dlabel.nii','.dscalar.nii','.dtseries.nii'})
    isCifti = 1;
else
    error('fpp.util.readDataMatrix is intended for NIFTI or dlabel/dscalar/dtseries CIFTI files.');
end

if isCifti
    hdr = cifti_read(inputPath);
    dataMat = hdr.cdata;
    hdr.cdata = [];
else
    hdr = fpp.util.mriRead(inputPath);
    dims = size(hdr.vol);
    if length(dims)<=3
        dataMat = hdr.vol(:);
    else
        dataMat = reshape(hdr.vol,[prod(dims(1:3)) dims(4)]);
    end
    if exist('maskVol','var') && ~isempty(maskVol)
        dataMat = dataMat(maskVol(:)==1,:);
    end
    hdr.vol = [];
end

end