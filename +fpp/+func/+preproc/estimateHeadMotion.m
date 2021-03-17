
% Function to estimate head motion using FSL's MCFLIRT
%
% motionParams = fpp.func.preproc.estimateHeadMotion(inputPath,outputDir,moCorrTargetVolNum)

function motionParams = estimateHeadMotion(inputPath,outputDir,moCorrTargetVolNum)

% TODO:
% - Add definition of json file sidecar for confounds.tsv.

if ~exist('moCorrTargetVolNum','var')
    vols = fpp.util.checkMRIProperty('vols',inputPath);
    moCorrTargetVolNum = ceil(vols/2);
end
[~,inputName,inputExt] = fpp.util.fileParts(inputPath);
mcBase = [outputDir '/' strrep(inputName,'_bold','') '_motion'];
confoundFile = [outputDir '/../' fpp.bids.changeName([inputName inputExt],{'desc','echo'},{[],[]},'confounds','.tsv')];
if exist(confoundFile,'file')
    tsvData = bids.util.tsvread(confoundFile);
end

% Run FSL's MCFLIRT
fpp.util.system(['mcflirt -in ' inputPath ' -out ' mcBase ...
    ' -mats -plots -refvol ' int2str(moCorrTargetVolNum)]);

% Load motion parameters
motionParams = load([mcBase '.par']);

% Add motion params to confounds.tsv data (rotation in degrees)
tsvData.rot_x = motionParams(:,1)*180/pi;
tsvData.rot_y = motionParams(:,2)*180/pi;
tsvData.rot_z = motionParams(:,3)*180/pi;
tsvData.trans_x = motionParams(:,4);
tsvData.trans_y = motionParams(:,5);
tsvData.trans_z = motionParams(:,6);

% Add mean translation/rotation in mm/degrees
moDiff = zeros(size(motionParams));
moDiff(2:end,:) = diff(motionParams);
tsvData.trans_total = sqrt(sum(moDiff(:,4:6).^2,2));
tsvData.rot_total = acos((cos(moDiff(:,1)).*cos(moDiff(:,2)) + cos(moDiff(:,1)).*cos(moDiff(:,3)) + ...
    cos(moDiff(:,2)).*cos(moDiff(:,3)) + sin(moDiff(:,1)).*sin(moDiff(:,2)).*sin(moDiff(:,3)) - 1)/2)*180/pi;

% Add FramewiseDisplacement
moDiffMM = moDiff;
moDiffMM(:,1:3) = 50*moDiffMM(:,1:3);   % Distance on 50mm-radius sphere, roughly the radius of cortex
tsvData.framewise_displacement = sum(abs(moDiffMM),2);

% Add DVARS
[tsvData.dvars,tsvData.dvars_std] = fpp.func.preproc.dvars(inputPath);

% Write to BIDS confounds.tsv file
fpp.bids.tsvWrite(confoundFile,tsvData);

% Rename xfm outputs
matFiles = dir([mcBase '.mat/MAT_*']);
for f=1:length(matFiles)
    fpp.util.system(['mv ' mcBase '.mat/' matFiles(f).name ' ' outputDir '/' ...
        fpp.bids.changeName(inputName,{'desc','from','to','mode'},{[],['native' ...
        matFiles(f).name(5:end)],'native','image'},'xfm','.mat')]);
end
fpp.util.system(['rm -rf ' mcBase '.nii.gz ' mcBase '.mat ' mcBase '.par']);


end