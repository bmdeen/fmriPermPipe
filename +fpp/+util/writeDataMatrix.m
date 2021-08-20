
% Wrapper function to write data stored as a coordinate by 1 vector or
% coordinate by time matrix to a NIFTI 3D/4D volume or CIFTI discalar, 
% dlabel, or dtseries file.
%
% writeDataMatrix(dataMat,hdr,outputPath,maskData)
%
% Arguments:
%   - dataMat (2D numeric matrix): coordinate by time data matrix, or
%       coordinate by 1 data vector
%   - hdr (struct): hdr structure read by fpp.util.readDataMatrix,
%       fpp.util.mriRead, or cifti_read
%   - outputPath (string): path to NIFTI/CIFTI output file
%
% Optional arguments:
% - maskData (vector or 3D matrix): mask data, either 3D volume (NIFTI)
%       or 1D vector (CIFTI) of 1s and 0s. dataMat will only contain values
%       where maskData==1.
%
% Dependencies: Freesurfer matlab toolbox, CIFTI matlab toolbox

function writeDataMatrix(dataMat,hdr,outputPath,maskData)

[~,~,outputExt] = fpp.util.fileParts(outputPath);
if ismember(lower(outputExt),{'.nii.gz','.nii'})
    isCifti = 0;
elseif ismember(lower(outputExt),{'.dlabel.nii','.dscalar.nii','.dtseries.nii'})
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

dtseriesFields = {'seriesStart','seriesStep','seriesUnit'};
dscalarFields = {'maps'};

if isCifti
    if useMask
        tmpMat = zeros(hdr.diminfo{1}.length,size(dataMat,2));
        tmpMat(maskData(1:hdr.diminfo{1}.length)==1,:) = dataMat;
        hdr.cdata = tmpMat;
    else
        hdr.cdata = dataMat;
        hdr.diminfo{1}.length = size(dataMat,1);    % Can't use if mask was applied
    end
    hdr.diminfo{2}.length = size(dataMat,2);
    if strcmpi(outputExt,'.dscalar.nii')
        hdr.diminfo{2}.type = 'scalars';
        for f=1:length(dtseriesFields)
            if isfield(hdr.diminfo{2},dtseriesFields{f})
                rmfield(hdr.diminfo{2},dtseriesFields{f});
            end
        end
        if ~isfield(hdr.diminfo{2},'maps')    % Add map names, if not specified
            for i=1:hdr.diminfo{2}.length
                hdr.diminfo{2}.maps(i).name = int2str(i);
                hdr.diminfo{2}.maps(i).metadata = struct('key',{},'value',{});
            end
        end
    elseif strcmpi(outputExt,'.dtseries.nii')
        hdr.diminfo{2}.type = 'series';
        for f=1:length(dscalarFields)
            if isfield(hdr.diminfo{2},dscalarFields{f})
                rmfield(hdr.diminfo{2},dscalarFields{f});
            end
        end
    elseif strcmpi(outputExt,'.dlabel.nii')
        hdr.diminfo{2}.type = 'labels';
        for f=1:length(dtseriesFields)
            if isfield(hdr.diminfo{2},dtseriesFields{f})
                rmfield(hdr.diminfo{2},dtseriesFields{f});
            end
        end
        if ~isfield(hdr.diminfo{2},'maps')    % Add map names, if not specified
            for i=1:hdr.diminfo{2}.length
                hdr.diminfo{2}.maps(i).name = int2str(i);
                hdr.diminfo{2}.maps(i).metadata = struct('key',{},'value',{});
            end
        end
    end
    cifti_write(hdr,outputPath);
else
    dims = hdr.volsize([2 1 3]);
    if size(dataMat,2)>1, dims(4) = size(dataMat,2); end
    hdr.vol = zeros(dims);
    if useMask
        if length(dims)==4
            for t=1:dims(4)
                tmp = zeros(dims(1:3));
                tmp(maskData(:)==1) = dataMat(:,t);
                hdr.vol(:,:,:,t) = tmp;
            end
        else
            hdr.vol(maskData(:)==1) = dataMat;
        end
    else
        hdr.vol = reshape(dataMat,dims);
    end
    fpp.util.mriWrite(hdr,outputPath);
end

end