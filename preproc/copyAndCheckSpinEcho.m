
% Function to copy spin echo field maps from raw to derivative directory,
% check phase encode information, and verify that all information required
% for topup-based distortion correction can be identified.
%
% If data are BIDS-formatted, and functional data have defined IntendedFor 
% metadata field, only the first two arguments are needed.

function [errorMsg,spinEchoPaths,fieldMapParamPath] = copyAndCheckSpinEcho(inputFuncPath,fmapPreprocDir,spinEchoPaths,...
    spinEchoPhaseEncodeDirections,funcDataPhaseEncodeDirection,fieldMapParamPath)

errorMsg = [];

[inputDir,inputName,~] = filepartsGZ(inputFuncPath);
if isempty(inputDir), inputDir = pwd; end
fmapRawDir = [inputDir '/../fmap'];

% Find spin echo field map files, if not specified, based on IntendedFor
% field in JSON sidecar.
if ~exist('spinEchoPaths','var') || isempty(spinEchoPaths)
    spinEchoPaths = {};
    fmapRawNames = dir([fmapRawDir '/*_epi.nii.gz']);
    for f=1:length(fmapRawNames)
        fmapRawPath = [fmapRawDir '/' fmapRawNames(f).name];
        fmapJsonPath = strrep(fmapRawPath,'.nii.gz','.json');
        if exist(fmapJsonPath,'file')
            jsonData = bids.util.jsondecode(fmapJsonPath);
            if isfield(jsonData,'IntendedFor')
                intendedFor = jsonData.IntendedFor;
                if ~iscell(intendedFor), intendedFor = {intendedFor}; end
                for i=1:length(intendedFor)
                    if any(strfind(intendedFor{i},inputName))
                        spinEchoPaths{end+1} = fmapRawPath;
                    end
                end
            end
        end
    end
elseif size(spinEchoPaths,1)>size(spinEchoPaths,2)
    spinEchoPaths = spinEchoPaths'; % Ensure row array
end

% If intended field maps couldn't be found, return
if length(spinEchoPaths)~=2
    errorMsg = 'ERROR: Two spin echo field maps intended for the functional data could not be found.';
    return;
end

% Check spin echo field map PE directions, if not specified
if ~exist('spinEchoPhaseEncodeDirections','var') || isempty(spinEchoPhaseEncodeDirections)
    for f=1:length(spinEchoPaths)
        tmp = checkMRIProperty('PEDir',spinEchoPaths{f});
        if ~isempty(tmp)
            spinEchoPhaseEncodeDirections{f} = tmp;
        end
    end
elseif size(spinEchoPhaseEncodeDirections,1)>size(spinEchoPhaseEncodeDirections,2)
    spinEchoPhaseEncodeDirections = spinEchoPhaseEncodeDirections'; % Ensure row array
end
% Convert PE directions to string output format (e.g. 'AP')
for f=1:length(spinEchoPaths)
    tmp = checkMRIProperty('PEDirStr',spinEchoPaths{f});
    if ~isempty(tmp)
        spinEchoPhaseEncodeDirectionsStr{f} = tmp;
    end
end

% Check functional data PE direction
if ~exist('funcDataPhaseEncodeDirection','var') || isempty(funcDataPhaseEncodeDirection)
    tmp = checkMRIProperty('PEDir',inputFuncPath);
    if ~isempty(tmp)
        funcDataPhaseEncodeDirection = tmp;
    end
end

% Check for appropriate PE directions
if length(spinEchoPhaseEncodeDirections)~=2
    errorMsg = 'ERROR: Spin echo phase encode direction information could not be found.';
    return;
end
if isempty(funcDataPhaseEncodeDirection)
    errorMsg = 'ERROR: Could not determine functional data phase encode direction.';
    return;
end
desiredPhaseEncodeDirections = {funcDataPhaseEncodeDirection(1),[funcDataPhaseEncodeDirection(1) '-']};
if ~(all(ismember(spinEchoPhaseEncodeDirections,desiredPhaseEncodeDirections)) && ...
        length(unique(spinEchoPhaseEncodeDirections))==2)
    errorMsg = 'ERROR: Spin echo acquisitions do not have approppriate PE direction.';
    return;
end

% Reorder spin echo files to make first one match PE dir of data
if ~strcmp(spinEchoPhaseEncodeDirections{1},funcDataPhaseEncodeDirection)
    spinEchoPhaseEncodeDirections = fliplr(spinEchoPhaseEncodeDirections);
    spinEchoPhaseEncodeDirectionsStr = fliplr(spinEchoPhaseEncodeDirectionsStr);
    spinEchoPaths = fliplr(spinEchoPaths);
end

% Copy spin echo images to fmapPreprocDir, and concatenate
for f=1:length(spinEchoPaths)
    [~,spinEchoName,~] = filepartsGZ(spinEchoPaths{f});
    spinEchoPathsNew{f} = [fmapPreprocDir '/' bidsChangeEntity(spinEchoName,'dir',spinEchoPhaseEncodeDirectionsStr{f},'epi') '.nii.gz'];
    copyImageAndJson(spinEchoPaths{f},spinEchoPathsNew{f});
end
spinEchoPaths = spinEchoPathsNew;
spinEchoPaths{3} = bidsChangeEntity(spinEchoPaths{1},'dir','Both');
system(['fslmerge -t ' spinEchoPaths{3} ' ' spinEchoPaths{1} ' ' spinEchoPaths{2}]);
inputJsonPath = strrep(spinEchoPaths{1},'.nii.gz','.json');
outputJsonPath = strrep(spinEchoPaths{3},'.nii.gz','.json');
jsonReconstruct(inputJsonPath,outputJsonPath);
jsonChangeValue(outputJsonPath,{'PhaseEncodingDirection'},{[]});

% Check topup properties, write to field map param file
if ~exist('fieldMapParamPath','var') || isempty(fieldMapParamPath) || ~exist(fieldMapParamPath,'file')
    fmapProperties = [];
    for f=1:2
        fmapProperties(end+1,:) = checkMRIProperty('Topup',spinEchoPaths{f});
    end
    if isempty(fmapProperties)
        errorMsg = 'ERROR: Could not determine spin echo topup parameters.';
        return;
    end
    fieldMapParamPath = strrep(spinEchoPaths{3},'_epi.nii.gz','_topupparams.txt');
    fid = fopen(fieldMapParamPath,'w');
    fprintf(fid,'%d %d %d %f\n',fmapProperties');
    fclose(fid);
end

end
    
    