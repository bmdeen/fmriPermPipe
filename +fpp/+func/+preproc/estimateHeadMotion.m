
% Function to estimate head motion using FSL's MCFLIRT
%
% ADD: computation of total translation/rotation, and DVARS. Definition of
% json file sidecar for confounds.tsv.

function motionParams = estimateHeadMotion(inputPath,outputDir,moCorrTargetVolNum)

if ~exist('moCorrTargetVolNum','var')
    vols = fpp.util.checkMRIProperty('vols',inputPath);
    moCorrTargetVolNum = ceil(vols/2);
end
[~,inputName,inputExt] = fpp.util.fileParts(inputPath);
mcBase = [outputDir '/' strrep(inputName,'_bold','') '_motion'];
confoundFile = [outputDir '/' fpp.bids.changeName([inputName inputExt],{},{},'confounds','.tsv')];
if exist(confoundFile,'file')
    tsvData = bids.util.tsvread(confoundFile);
end

% Run FSL's MCFLIRT
fpp.util.system(['mcflirt -in ' inputPath ' -out ' mcBase ...
    ' -mats -plots -refvol ' int2str(moCorrTargetVolNum)]);

% Load motion parameters
motionParams = load([mcBase '.par']);

% Write to BIDS confounds.tsv file
tsvData.trans_x = motionParams(:,1);
tsvData.trans_y = motionParams(:,2);
tsvData.trans_z = motionParams(:,3);
tsvData.rot_x = motionParams(:,4);
tsvData.rot_y = motionParams(:,5);
tsvData.rot_z = motionParams(:,6);
fpp.bids.tsvWrite(confoundFile,tsvData);

% Rename xfm outputs
matFiles = dir([mcBase '.mat/MAT_*']);
for f=1:length(matFiles)
    fpp.util.system(['mv ' mcBase '.mat/' matFiles(f).name ' ' outputDir '/' ...
        fpp.bids.changeName(inputName,{'from','to','mode'},{['orig' ...
        matFiles(f).name(5:end)],'orig','image'},'xfm','.mat')]);
end
fpp.util.system(['rm -rf ' mcBase '.nii.gz ' mcBase '.mat ' mcBase '.par']);


end