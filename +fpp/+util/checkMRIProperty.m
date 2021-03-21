%
% propertyValue = fpp.util.checkMRIProperty(propertyName,inputPath)
%
% Function to check property of input MR image, using sidecar JSON file if
% it exists, or alternative method if available. Returns null if it can't
% find a way to check the property.
%
% Arguments:
% - propertyName (string): property to check
%   Options:
%   + TR - Repetition time (s)
%   + TE - Echo time (ms), in vector form for multi-echo data
%   + Vols - # of volumes in 4D dataset
%   + Dims - 3D or 4D image dimensions
%   + VoxelSize - voxel sizes in each dimension (mm)
%   + PEDir - Phase-encode direction, BIDS format (e.g. "j-")
%   + PEDirStr - Phase-encode direction, orientation string format (e.g.
%       "AP")
%   + SEDir - Slice encoding direction, BIDS format
%   + SEDirStr - Slice-encode direction, orientation string format (e.g.
%       "AP")
%   + Topup - Spin-echo EPI properties for FSL's topup (phase dir/timing)
%   + ST/SliceTiming - vector of slice acquisition times relative to start
%       of volume acquisition
% - inputPath (string): path to input image
%
% Dependencies: bids-matlab

function propertyValue = checkMRIProperty(propertyName,inputPath)

propertyValue = [];
if exist(fpp.bids.jsonPath(inputPath),'file')
    jsonData = fpp.bids.getMetadata(inputPath);
end
[inputDir,inputName,inputExt] = fpp.util.fileParts(inputPath);
if isempty(inputDir), inputDir = pwd; end

switch lower(propertyName)
    case 'tr'
        if exist('jsonData','var') && isfield(jsonData,'RepetitionTime')
            propertyValue = jsonData.RepetitionTime;
        else
            [~, tr] = fpp.util.system(['fslval ' inputPath ' pixdim4']);
            [~, tu] = fpp.util.system(['fslval ' inputPath ' time_units']);
            tu = strtrim(tu);
            if strcmp(tu,'s')
                propertyValue = str2num(tr);
            elseif strcmp(tu,'ms')
                propertyValue = str2num(tr)/1000;
            elseif strcmp(tu,'us')
                propertyValue = str2num(tr)/1000000;
            elseif strcmp(tu,'Unknown')
                propertyValue = [];
            end
        end
    case 'te'
        if exist('jsonData','var') && isfield(jsonData,'EchoTime')
            % Check for BIDS-formatted multi-echo data, return multiple TE values if so
            if any(regexp(inputPath,'_echo-[0-9]+_'))
                % Get list of input paths
                inputNames = struct2cell(fpp.util.regExpDir(inputDir,regexprep([inputName inputExt],'_echo-[0-9]+_','_echo-[0-9]+_')));
                inputNames = inputNames(1,:)';
                for e=1:length(inputNames)
                    inputPaths{e} = [inputDir '/' inputNames{e}];
                    jsonData = fpp.bids.getMetadata(inputPaths{e});
                    propertyValue(e) = jsonData.EchoTime*1000;
                end
            else
                propertyValue = jsonData.EchoTime*1000;
            end
        end
    case 'pedir'
        if exist('jsonData','var') && isfield(jsonData,'PhaseEncodingDirection')
            propertyValue = jsonData.PhaseEncodingDirection;
        end
    case 'pedirstr'
        if exist('jsonData','var') && isfield(jsonData,'PhaseEncodingDirection')
            propertyValue = fpp.util.convertBidsPEDirToStr(inputPath,jsonData.PhaseEncodingDirection);
        end
    case 'sedir'
        if exist('jsonData','var') && isfield(jsonData,'SliceEncodingDirection')
            propertyValue = jsonData.SliceEncodingDirection;
        end
    case 'sedirstr'
        if exist('jsonData','var') && isfield(jsonData,'SliceEncodingDirection')
            imageOrientation = fpp.util.getImageOrientation(inputPath);
            orientationLabels = {'L','R','A','P','S','I'};
            orientationLabelsInverted = {'R','L','P','A','I','S'};
            switch jsonData.SliceEncodingDirection(1)
                case 'i'
                    seDir = 1;
                case 'j'
                    seDir = 2;
                case 'k'
                    seDir = 3;
            end
            propertyValue = [orientationLabelsInverted{strcmp(imageOrientation(seDir),orientationLabels)} ...
                imageOrientation(seDir)];
            if length(jsonData.SliceEncodingDirection)>1
                propertyValue = fliplr(propertyValue);
            end
        end
    case 'topup'
        if exist('jsonData','var')  && isfield(jsonData,'PhaseEncodingDirection') ...
        	&& isfield(jsonData,'TotalReadoutTime')
            switch jsonData.PhaseEncodingDirection(1)
                case 'i'
                    phaseEncodeVector = [1 0 0];
                case 'j'
                    phaseEncodeVector = [0 1 0];
                case 'k'
                    phaseEncodeVector = [0 0 1];
            end
            if length(jsonData.PhaseEncodingDirection)>1
                phaseEncodeVector = -phaseEncodeVector;
            end
            propertyValue = [phaseEncodeVector jsonData.TotalReadoutTime];
        end
    case {'st','slicetiming'}
        if exist('jsonData','var') && isfield(jsonData,'SliceTiming')
            propertyValue = jsonData.SliceTiming;
        end
    case 'vols'
        [~,vols] = fpp.util.system(['fslval ' inputPath ' dim4']);
        propertyValue = str2num(strtrim(vols));
    case 'dims'
        [~,dim1] = fpp.util.system(['fslval ' inputPath ' dim1']);
        propertyValue(1) = str2num(strtrim(dim1));
        [~,dim2] = fpp.util.system(['fslval ' inputPath ' dim2']);
        propertyValue(2) = str2num(strtrim(dim2));
        [~,dim3] = fpp.util.system(['fslval ' inputPath ' dim3']);
        propertyValue(3) = str2num(strtrim(dim3));
        [~,vols] = fpp.util.system(['fslval ' inputPath ' dim4']);
        vols = str2num(strtrim(vols));
        if ~isempty(vols) && vols>1, propertyValue(4) = vols; end
    case {'voxelsize','voxsize'}
        [~,pixdim1] = fpp.util.system(['fslval ' inputPath ' pixdim1']);
        propertyValue(1) = str2num(strtrim(pixdim1));
        [~,pixdim2] = fpp.util.system(['fslval ' inputPath ' pixdim2']);
        propertyValue(2) = str2num(strtrim(pixdim2));
        [~,pixdim3] = fpp.util.system(['fslval ' inputPath ' pixdim3']);
        propertyValue(3) = str2num(strtrim(pixdim3));
        [~, vu] = fpp.util.system(['fslval ' inputPath ' vox_units']);
        vu = strtrim(vu);
        if strcmp(vu,'m')
            propertyValue = propertyValue*1000;
        elseif strcmp(vu,'um')
            propertyValue = propertyValue/1000;
        elseif strcmp(vu,'Unknown')
            propertyValue = [];
        end
end