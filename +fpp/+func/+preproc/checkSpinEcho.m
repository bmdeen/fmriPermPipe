
% Function to determine which spin echo and undistortion warp/jacobian 
% files to use for a given raw fMRI dataset, based on IntendedFor info in 
% spin-echo image json metadata. Should be run after fpp.fmap.preproc.
%
% [topupWarpPath,topupJacobianPath] = fpp.func.preproc.checkSpinEcho(inputFuncPath,fmapPreprocDir)

function [spinEchoPath,topupWarpPath,topupJacobianPath] = checkSpinEcho(inputFuncPath,fmapPreprocDir)

topupWarpPath = '';
topupJacobianPath = '';
spinEchoPath = '';
spinEchoPaths = {};

[inputDir,inputName,~] = fpp.util.fileParts(inputFuncPath);
if isempty(inputDir), inputDir = pwd; end

% Find spin echo field map files based on IntendedFor field in JSON sidecar
spinEchoPaths = {};
fmapNames = dir([fmapPreprocDir '/*_epi.nii.gz']);
for f=1:length(fmapNames)
    fmapPath = [fmapPreprocDir '/' fmapNames(f).name];
    if strfind(fmapNames(f).name,'dir-Both'), continue; end  % Don't include combined field map files
    fmapMetadata = fpp.bids.getMetadata(fmapPath);
    if isfield(fmapMetadata,'IntendedFor')
        intendedFor = fmapMetadata.IntendedFor;
        if ~iscell(intendedFor), intendedFor = {intendedFor}; end
        for i=1:length(intendedFor)
            if any(strfind(intendedFor{i},inputName))
                spinEchoPaths{end+1} = fmapPath;
            end
        end
    end
end

% If intended field maps couldn't be found, return
if isempty(spinEchoPaths)
    error('Spin echo field maps intended for the functional data could not be found.');
end
[~,spinEchoName,~] = fpp.util.fileParts(spinEchoPaths{1});

% Check functional data phase encode direction
peDirStr = fpp.util.checkMRIProperty('pedirstr',inputFuncPath);

% Define warp/jacobian paths
if ~isempty(peDirStr)
    spinEchoPath = [fmapPreprocDir '/' fpp.bids.changeName(spinEchoName,'dir',peDirStr) '.nii.gz'];
    topupWarpPath = [fmapPreprocDir '/' fpp.bids.changeName(spinEchoName,{'from','to','mode','desc','dir'},...
        {'native','Undistorted','image',['UndistortionWarp' peDirStr],''},'xfm','.nii.gz')];
    topupJacobianPath = fpp.bids.changeName(topupWarpPath,{},{},'jacobian');
else
    warning(['Phase-encode direction could not be determined for functional data ' ...
        inputFuncPath ' - assuming that it is matched to undistortion warp direction. Check this!!']);
    topupWarpPathTemplate = [fmapPreprocDir '/' fpp.bids.changeName(spinEchoName,{'from','to','mode','desc','dir'},...
        {'native','Undistorted','image','UndistortionWarp*',''},'xfm','.nii.gz')];
    warpFiles = dir(topupWarpPathTemplate);
    if isempty(warpFiles)
        error('Undistortion warp could not be found. Run fpp.fmap.preproc.');
    elseif length(warpFiles)>1
        error(['Multiple undistortion warp files exist for a given run; '...
            'functional data phase encode direction must be specified.']);
    end
    topupWarpPath = [fmapPreprocDir '/' warpFiles(1).name];
    topupJacobianPath = fpp.bids.changeName(topupWarpPath,{},{},'jacobian');
    [~,ind] = regexpi(topupWarpPath,'desc-UndistortionWarp');
    peDirStr = topupWarpPath(ind+1:ind+2);
    spinEchoPath = [fmapPreprocDir '/' fpp.bids.changeName(spinEchoName,'dir',peDirStr) '.nii.gz'];
end

end
