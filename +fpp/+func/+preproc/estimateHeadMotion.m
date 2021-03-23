
% Function to estimate head motion using FSL's MCFLIRT. Generates
% transformation matrices to be used for later one-shot motion and
% distortion correction.
%
% confoundTSV = fpp.func.preproc.estimateHeadMotion(inputPath,outputDir,moCorrTargetVolNum)

function confoundTSV = estimateHeadMotion(inputPath,outputDir,moCorrTargetVolNum)

% TODO:
% - Add definition of json file sidecar for confounds.tsv.

if ~exist('moCorrTargetVolNum','var')
    vols = fpp.util.checkMRIProperty('vols',inputPath);
    moCorrTargetVolNum = ceil(vols/2);
end
[~,inputName,inputExt] = fpp.util.fileParts(inputPath);
mcBase = [outputDir '/' strrep(inputName,'_bold','') '_motion'];

% Run FSL's MCFLIRT
fpp.util.system(['mcflirt -in ' inputPath ' -out ' mcBase ...
    ' -mats -plots -refvol ' int2str(moCorrTargetVolNum)]);

% Load motion parameters
motionParams = load([mcBase '.par']);

% Add motion params to confound TSV data (mm/degrees)
confoundTSV.rot_x = motionParams(:,1)*180/pi;
confoundTSV.rot_y = motionParams(:,2)*180/pi;
confoundTSV.rot_z = motionParams(:,3)*180/pi;
confoundTSV.trans_x = motionParams(:,4);
confoundTSV.trans_y = motionParams(:,5);
confoundTSV.trans_z = motionParams(:,6);

% Compute framewise translation/rotation (mm/degrees)
moDiff = zeros(size(motionParams));
moDiff(2:end,:) = diff(motionParams);
confoundTSV.framewise_translation = sqrt(sum(moDiff(:,4:6).^2,2));
confoundTSV.framewise_rotation = acos((cos(moDiff(:,1)).*cos(moDiff(:,2)) + cos(moDiff(:,1)).*cos(moDiff(:,3)) + ...
    cos(moDiff(:,2)).*cos(moDiff(:,3)) + sin(moDiff(:,1)).*sin(moDiff(:,2)).*sin(moDiff(:,3)) - 1)/2)*180/pi;

% Compute framewise displacement (mm)
moDiffMM = moDiff;
moDiffMM(:,1:3) = 50*moDiffMM(:,1:3);   % Rotation-driven distance on 50mm-radius sphere, roughly the radius of cortex
confoundTSV.framewise_displacement = sum(abs(moDiffMM),2);

% Rename xfm outputs
matFiles = dir([mcBase '.mat/MAT_*']);
for f=1:length(matFiles)
    fpp.util.system(['mv ' mcBase '.mat/' matFiles(f).name ' ' outputDir '/' ...
        fpp.bids.changeName(inputName,{'desc','from','to','mode'},{[],['native' ...
        matFiles(f).name(5:end)],'native','image'},'xfm','.mat')]);
end
fpp.util.system(['rm -rf ' mcBase '.nii.gz ' mcBase '.mat ' mcBase '.par']);


end