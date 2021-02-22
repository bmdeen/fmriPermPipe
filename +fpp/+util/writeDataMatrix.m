
% Wrapper function to write data stored as a coordinate by 1 vector or
% coordinate by time matrix to a NIFTI 3D/4D volume or CIFTI discalar or
% dtseries file.
%
% writeDataMatrix(dataMat,hdr,outputPath,maskPath)
%
% Arguments:
%   - dataMat (2D numeric matrix): coordinate by time data matrix, or
%       coordinate by 1 data vector
%   - hdr (struct): hdr structure read by fpp.util.readDataMatrix,
%       fpp.util.mriRead, or cifti_read
%   - outputPath (string): path to NIFTI/CIFTI output file
%   - maskVol (3D matrix, optional): 3D mask volume (NIFTI only)
%
% Dependencies: Freesurfer matlab toolbox, CIFTI matlab toolbox

function writeDataMatrix(dataMat,hdr,outputPath,maskVol)

[~,~,outputExt] = fpp.util.fileParts(outputPath);
if ismember(lower(outputExt),{'.nii.gz','.nii'})
    isCifti = 0;
elseif ismember(lower(outputExt),{'.dscalar.nii','.dtseries.nii'})
    isCifti = 1;
else
    error('fpp.util.readDataMatrix is intended for NIFTI or dscalar/dtseries CIFTI files.');
end

dtseriesFields = {'seriesStart','seriesStep','seriesUnit'};
dscalarFields = {'maps'};

if isCifti
    hdr.cdata = dataMat;
    hdr.diminfo{1}.length = size(dataMat,1);
    hdr.diminfo{2}.length = size(dataMat,2);
    if strcmpi(outputExt,'.dscalar.nii')
        hdr.diminfo{2}.type = 'scalars';
        for f=1:length(dtseriesFields)
            if isfield(hdr.diminfo{2},dtseriesFields{f})
                rmfield(hdr.diminfo{2},dtseriesFields{f});
            end
        end
    elseif strcmpi(outputExt,'.dtseries.nii')
        hdr.diminfo{2}.type = 'series';
        for f=1:length(dscalarFields)
            if isfield(hdr.diminfo{2},dscalarFields{f})
                rmfield(hdr.diminfo{2},dscalarFields{f});
            end
        end
    end
    cifti_write(hdr,outputPath);
else
    dims = hdr.volsize;
    if size(dataMat,2)>1, dims(4) = size(dataMat,2); end
    hdr.vol = zeros(dims);
    if exist('maskVol','var') && ~isempty(maskVol)
        if length(dims)==4
            for t=1:dims(4)
                tmp = zeros(dims(1:3));
                tmp(maskVol(:)==1) = dataMat(:,t);
                hdr.vol(:,:,:,t) = tmp;
            end
        else
            hdr.vol(maskVol(:)==1) = dataMat;
        end
    else
        hdr.vol = reshape(dataMat,dims);
    end
    fpp.util.mriWrite(hdr,outputPath);
end

end