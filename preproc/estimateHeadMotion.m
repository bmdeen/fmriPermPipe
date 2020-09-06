
% Function to estimate head motion using FSL's MCFLIRT

function motionParams = estimateHeadMotion(inputPath,outputDir,moCorrTargetVolNum)

if ~exist('moCorrTargetVolNum','var')
    [~,vols] = system(['fslval ' inputPath ' dim4']);
    vols = str2num(strtrim(vols));
    moCorrTargetVolNum = ceil(vols/2);
end
[~,inputName,inputExt] = filepartsGZ(inputPath);
mcBase = [outputDir '/' strrep(inputName,'_bold','') '_motion'];

% Run FSL's MCFLIRT
system(['mcflirt -in ' inputPath ' -out ' mcBase ...
    ' -mats -plots -refvol ' int2str(moCorrTargetVolNum)]);

% Rename outputs
system(['mv ' mcBase '.par ' mcBase '.tsv']);   % Needs column labels?
matFiles = dir([mcBase '.mat/MAT_*']);
for f=1:length(matFiles)
    system(['mv ' mcBase '.mat/' matFiles(f).name ' ' outputDir '/' ...
        bidsChangeEntity(inputName,{'from','to','mode'},{['NativeFunc' ...
        matFiles(f).name(5:end)],'NativeFunc','image'},'xfm','.mat')]);
end
system(['rm -rf ' mcBase '.nii.gz ' mcBase '.mat']);

motionParams = load([mcBase '.tsv']);

end