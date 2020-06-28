
% registerMRI(subject,targetFunc,highresHead,varargin)
%
% Computes registration files between various spaces, including the target
% functional space used for data analysis, the high-resolution anatomical 
% image, the same anatomical image in the orientation/resolution used by 
% Freesurfer, and MNI152 standard neurological space. Should be run after 
% preprocMRI and recon-all. Steps include:
% - Register target functional to anatomical using bbregister
% - Register anatomical to MNI space using FNIRT, and to Freesurfer
%   anatomical space using FLIRT
% - Concatenate registrations between multiple spaces where appropriate
% - Generates gray matter, white matter, and cerebrospinal fluid masks from
%   Freesurfer parcellations. Converts Freesurfer-based images (masks and
%   parcellations) to anatomical and target spaces. Masks anatomical using 
%   Freesurfer brain mask.
%
% Example usage: registerMRI('SUB1','/studyDirectory/SUB1/func/targetFunctional.nii.gz',...
%                    '/studyDirectory/SUB1/func/highres_anatomical.nii.gz');
%
% Output naming conventions (consistent with FSL/Freesurfer when possible):
% - target_func: functional image that was the registration target for
%   functional data in preproc, located in func directory
% - highres: anatomical image, located in anat directory
% - orig: anatomical image in Freesurfer space
% - standard: MNI space
% - imageA2imageB.mat: registration file to transform imageA to imageB
% - imageA2imageB.nii.gz: imageA registered to imageB space
%
% Required arguments:
% - subject (string): name of subject, used to specify Freesurfer subject
%       directory
% - targetFunc (string): path to target functional image
% - highresHead (string): path to whole-head T1-weighted anatomical image
%
% Variable arguments:
% - overwrite (boolean; default=0): whether to overwrite files that have
%       already been written by this function.
% - anatRegDir (string): output directory for anatomical registration files
% - funcRegDir (string): output direectory for functional registration files
% - fsDir (string): freesurfer directory
% - funcRoiDir (string): output directory for func-space Freesurfer ROIs
% - anatRoiDir (string): output directory for anat-space Freesurfer ROIs

function registerMRI(subject,targetFunc,highresHead,varargin)

addpath([strrep(mfilename('fullpath'),mfilename,'') '/utils']);

% Load/check config variables.
[configError, fslPrefix] = checkConfig;
if ~isempty(configError)
    fprintf('%s\n',configError);
    return;
end

% Check if input images exist
if exist(targetFunc,'file')==0
    fprintf('%s\n\n',['Target functional ' targetFunc ' does not exist.']);
    return;
elseif exist(highresHead,'file')==0
    fprintf('%s\n\n',['Anatomical image ' thighresHead ' does not exist.']);
    return;
end

overwrite = 0;
[anatDir,~,highresExt] = fileparts(highresHead);
[~,targetName,~] = fileparts(targetFunc);
targetName = strrep(targetName,'.nii','');
fsDir = [anatDir '/../../freesurfer'];
anatRegDir = [anatDir '/../reg'];
funcRegDir = [anatRegDir '/' targetName];
funcRoiDir = [anatDir '/../roi/' targetName];
anatRoiDir = anatDir;

% Edit variable arguments.  Note: optInputs checks for proper input.
varArgList = {'overwrite','anatRegDir','funcRegDir','fsDir','funcRoiDir','anatRoiDir'};
for i=1:length(varArgList)
    argVal = optInputs(varargin,varArgList{i});
    if ~isempty(argVal)
        eval([varArgList{i} ' = argVal;']);
    end
end

fsPrefix = ['export SUBJECTS_DIR=' fsDir '; '];
parcNames = {'aparc+aseg','aparc.a2009s+aseg'}; % Freesurfer parcellation names.  Note: script assumes #2 is aparc.a2009s+aseg.

% Freesurfer labels for CSF and WM (aparc.a2009s+aseg parcellation)
csfLabelNums = {4:5,14:15,43:44,72,75:76};
wmLabelNums = {2,7,41,46,77:79,85,100,110,155:158,250:255};

fslDir = getenv('FSL_DIR');
standard = [fslDir '/data/standard/MNI152_T1_2mm_brain.nii.gz'];
standardHead = [fslDir '/data/standard/MNI152_T1_2mm.nii.gz'];
standardMask = [fslDir '/data/standard/MNI152_T1_2mm_brain_mask_dil.nii.gz'];

if ~exist(anatRegDir,'dir'), mkdir(anatRegDir); end
if ~exist(funcRegDir,'dir'), mkdir(funcRegDir); end
if ~exist(funcRoiDir,'dir'), mkdir(funcRoiDir); end

% Convert highresHead to .nii.gz, if it's not already
if ~contains(lower(highresHead),'.nii.gz')
    system(['mri_convert ' highresHead ' ' strrep(highresHead,highresExt,'.nii.gz')]);
    highresHead = strrep(highresHead,highresExt,'.nii.gz');
end

highres = strrep(highresHead,'.nii.gz','_brain.nii.gz');
highres2mm = strrep(highresHead,'.nii.gz','_brain_2mm.nii.gz');
target = targetFunc;

% FS mri files
fsSubjDir = [fsDir '/' subject];
orig = [fsSubjDir '/mri/orig.nii.gz'];
origMGZ = strrep(orig,'.nii.gz','.mgz');
norm = [fsSubjDir '/mri/norm.nii.gz'];
normMGZ = strrep(norm,'.nii.gz','.mgz');
fsMask = [fsSubjDir '/mri/brainmask.nii.gz'];
fsMaskMGZ = strrep(fsMask,'.nii.gz','.mgz');
fsMaskHR = [anatRoiDir '/fsmask.nii.gz'];
fsMaskTarg = [funcRoiDir '/fsmask.nii.gz'];
% Parcellations
for p = 1:length(parcNames)
    parcs{p} = [fsSubjDir '/mri/' parcNames{p} '.nii.gz'];
    parcsMGZ{p} = strrep(parcs{p},'.nii.gz','.mgz');
    parcsHR{p} = [anatRoiDir '/' parcNames{p} '.nii.gz'];
    parcsTarg{p} = [funcRoiDir '/' parcNames{p} '.nii.gz'];
end
% GM, WM, CSF masks
maskNames = {'gmmask','wmmask','csfmask'};
for m=1:3
    masks{m} = [fsSubjDir '/mri/' maskNames{m} '.nii.gz'];
    masksTarg{m} = [funcRoiDir '/' maskNames{m} '.nii.gz'];
    masksHR{m} = [anatRoiDir '/' maskNames{m} '.nii.gz'];
    masksTargErode{m} = [funcRoiDir '/' maskNames{m} '_erode1.nii.gz'];
end

% freesurfer <> highres registrations
freesurfer2Highres = [anatRegDir '/freesurfer2highres.mat'];
highres2Freesurfer = [anatRegDir '/highres2freesurfer.mat'];
highres2FreesurferFS = [anatRegDir '/highres2freesurfer.dat'];

% Anatomical <> standard registrations
highres2Std = [anatRegDir '/highres2standard.mat'];
std2Highres = [anatRegDir '/standard2highres.mat'];
highres2StdImg = [anatRegDir '/highres2standard.nii.gz'];
highres2StdWarp = [anatRegDir '/highres2standard_warp.nii.gz'];
highres2StdJac = [anatRegDir '/highres2standard_jac.nii.gz'];
std2HighresWarp = [anatRegDir '/standard2highres_warp.nii.gz'];

% Target registrations
target2Highres = [funcRegDir '/target_func2highres.mat'];
target2Highres_flirt = [funcRegDir '/target_func2highres.init.mat'];   % FLIRT-based initializing transform.
target2FreesurferInit = [funcRegDir '/target_func2freesurfer.init.mat'];
target2FreesurferInitFS = [funcRegDir '/target_func2freesurfer.init.dat'];
highres2Target = [funcRegDir '/highres2target_func.mat'];
target2Freesurfer = [funcRegDir '/target_func2freesurfer.mat'];
freesurfer2Target = [funcRegDir '/freesurfer2target_func.mat'];
target2FreesurferFS = [funcRegDir '/target_func2freesurfer.dat'];
target2HighresImg = [funcRegDir '/target_func2highres.nii.gz'];
target2Highres_flirtImg = [funcRegDir '/target_func2highres.init.nii.gz'];
target2StdImg = [funcRegDir '/target_func2standard.nii.gz'];

% Check if Freesurfer recon has run to completion.
if ~exist(parcsMGZ{1},'file')
    fprintf('%s\n\n',['Freesurfer recon has not completed for subject ' subject '.']);
    return
end

fprintf('%s\n',['Running registration for ' subject '.']);

% Convert MGZ to NII
if ~exist(norm,'file'), system(['mri_convert ' normMGZ ' ' norm]); end
if ~exist(orig,'file'), system(['mri_convert ' origMGZ ' ' orig]); end
if ~exist(fsMask,'file'), system(['mri_convert ' fsMaskMGZ ' ' fsMask]); end

% Generate GM/WM/CSF masks in Freesurfer space.
if ~exist(masks{end},'file') || overwrite
    
    for p = 1:length(parcNames)
        system(['mri_convert ' parcsMGZ{p} ' ' parcs{p}]);
    end
    
    % GM mask
    system([fslPrefix 'fslmaths ' parcs{2} ' -thr 10000 -bin ' masks{1}]);
    
    % WM mask
    tmpPath = [fsSubjDir '/mri/tmp_gen_masks'];
    cmd = [fslPrefix 'fslmaths'];
    for n = 1:length(wmLabelNums)
        system([fslPrefix 'fslmaths ' parcs{2} ' -thr ' int2str(min(wmLabelNums{n})) ...
            ' -uthr ' int2str(max(wmLabelNums{n})) ' -bin ' tmpPath '-' int2str(n) '.nii.gz']);
        if n==1
            cmd = [cmd ' ' tmpPath '-' int2str(n) '.nii.gz'];
        else
            cmd = [cmd ' -add ' tmpPath '-' int2str(n) '.nii.gz'];
        end
    end
    cmd = [cmd ' -bin ' masks{2} '; rm -rf ' tmpPath '*'];
    system(cmd);
    
    % CSF mask
    cmd = [fslPrefix 'fslmaths'];
    for n = 1:length(csfLabelNums)
        system([fslPrefix 'fslmaths ' parcs{2} ' -thr ' int2str(min(csfLabelNums{n})) ...
            ' -uthr ' int2str(max(csfLabelNums{n})) ' -bin ' tmpPath '-' int2str(n) '.nii.gz']);
        if n==1
            cmd = [cmd ' ' tmpPath '-' int2str(n) '.nii.gz'];
        else
            cmd = [cmd ' -add ' tmpPath '-' int2str(n) '.nii.gz'];
        end
    end
    cmd = [cmd ' -bin ' masks{3} '; rm -rf ' tmpPath '*'];
    system(cmd);
    
end

% Highres/orig registrations (FLIRT)
if ~(exist(highres2Freesurfer,'file') && exist(fsMaskHR,'file')) || overwrite
    % Register orig to highres.
    fprintf('%s\n','Registering orig to highres...');
    system([fslPrefix 'flirt -in ' orig ' -ref ' highresHead ' -omat ' freesurfer2Highres ' -cost normcorr ' ...
        '-dof 6 -searchrx -180 180 -searchry -180 180 -searchrz -180 180']);
    system([fslPrefix 'convert_xfm -omat ' highres2Freesurfer ' -inverse ' freesurfer2Highres]);
    system(['tkregister2 --mov ' highresHead ' --targ ' orig ' --s ' subject ' --fsl ' highres2Freesurfer ' --reg ' highres2FreesurferFS ' --noedit']);
    
    % Mask anatomical using Freesurfer-based mask.
    fprintf('%s\n','Masking anatomical image using Freesurfer mask...');
    system([fslPrefix 'flirt -in ' fsMask ' -ref ' highresHead ' -applyxfm -init ' ...
        freesurfer2Highres ' -out ' fsMaskHR ' -interp nearestneighbour']);
    system([fslPrefix 'fslmaths ' fsMaskHR ' -bin ' fsMaskHR]);
    system([fslPrefix 'fslmaths ' highresHead ' -mul ' fsMaskHR ' ' highres]);
end

% Highres/standard registrations (FNIRT)
% NOTE: FNIRT command mysteriously failing
if ~(exist(highres2StdWarp,'file') && exist(std2HighresWarp,'file')) || overwrite
    
    % Anatomical to MNI registration using FLIRT/FNIRT.
    fprintf('%s\n','Registering anatomical image to MNI space...');
    system([fslPrefix 'flirt -ref ' standard ' -in ' highres ' -omat ' highres2Std ...
        ' -cost corratio -dof 12 -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -interp trilinear']);
%     system([fslPrefix 'fnirt --in=' highresHead ' --aff=' highres2Std ' --cout=' highres2StdWarp ' --iout=' ...
%         highres2StdImg ' --jout=' highres2StdJac ' --config=T1_2_MNI152_2mm --ref=' standardHead ...
%         ' --refmask=' standardMask ' --warpres=10,10,10']);
    system([fslPrefix 'convert_xfm -omat ' std2Highres ' -inverse ' highres2Std]);
    
%     fprintf('%s\n','Inverting warp field...');
    system([fslPrefix 'fslmaths ' highres ' -subsamp2 ' highres2mm]);
%     system([fslPrefix 'invwarp --ref=' highres2mm ' --warp=' highres2StdWarp ' --out=' std2HighresWarp]);
    
end

% Target/highres registrations (bbregister)
if ~exist(target2Highres,'file') || overwrite
    % Find initializing FLIRT registration, and convert to Freesurfer .dat format.
    fprintf('%s\n','Running FLIRT for target2Highres reg...');
    system([fslPrefix 'flirt -in ' target ' -ref ' highres ' -out ' target2Highres_flirtImg ...
        ' -omat ' target2Highres_flirt ' -dof 6 -searchrx -10 10 -searchry -50 50 -searchrz -10 10']);
    system([fslPrefix 'convert_xfm -omat ' target2FreesurferInit ' -concat ' highres2Freesurfer ' ' target2Highres_flirt]);
    system(['tkregister2 --mov ' target ' --targ ' orig ' --s ' subject ' --fsl ' target2FreesurferInit ' --reg ' target2FreesurferInitFS ' --noedit']);
    
    % Perform bbregister
    fprintf('%s\n','Running bbregister for target2Highres reg...');
    system([fsPrefix 'bbregister --s ' subject ' --mov ' target ' --reg ' target2FreesurferFS ' --bold --init-reg ' target2FreesurferInitFS]);
    system([fsPrefix 'tkregister2 --mov ' target ' --targ ' orig ' --reg ' target2FreesurferFS ' --fslregout ' target2Freesurfer ' --noedit']);
    system([fslPrefix 'convert_xfm -omat ' freesurfer2Target ' -inverse ' target2Freesurfer]);
    system([fslPrefix 'convert_xfm -omat ' target2Highres ' -concat ' freesurfer2Highres ' ' target2Freesurfer]);
    system([fslPrefix 'convert_xfm -omat ' highres2Target ' -inverse ' target2Highres]);
    
    % Generate target2Highres image
    system([fslPrefix 'flirt -in ' target ' -ref ' highres ' -applyxfm -init ' ...
        target2Highres ' -out ' target2HighresImg]);
    
    % Generate target2standard
%     system([fslPrefix 'applywarp --ref=' standard ' --in=' target ' --out=' ...
%         target2StdImg ' --warp=' highres2StdWarp ' --premat=' target2Highres]);
    
    system([fslPrefix 'flirt -in ' fsMask ' -ref ' target ' -applyxfm -init ' ...
        freesurfer2Target ' -out ' fsMaskTarg ' -interp nearestneighbour']);
    system([fslPrefix 'fslmaths ' fsMaskTarg ' -bin ' fsMaskTarg]);
end

% Register parcellations and masks to target space.
if ~exist(masksTarg{end},'file') || overwrite
    
    for p = 1:length(parcNames)
        fprintf('%s\n',['Transforming ' parcNames{p} ' to target space.']);
        system([fslPrefix 'flirt -in ' parcs{p} ' -ref ' target ' -applyxfm -init ' ...
            freesurfer2Target ' -out ' parcsTarg{p} ' -interp nearestneighbour']);
    end
    
    for m = 1:length(masks)
        fprintf('%s\n',['Transforming ' maskNames{m} ' to target space.']);
        system([fslPrefix 'flirt -in ' masks{m} ' -ref ' target ' -applyxfm -init ' ...
            freesurfer2Target ' -out ' masksTarg{m} ' -interp nearestneighbour']);
    end
    
    % Erode WM/CSF masks by one voxel
    for m=2:3
        system([fslPrefix 'fslmaths ' masksTarg{m} ' -ero ' masksTargErode{m}]);
    end
    
end

% Register parcellations and masks to highres space.
if ~exist(masksHR{end},'file') || overwrite

    for p = 1:length(parcNames)
        fprintf('%s\n',['Transforming ' parcNames{p} ' to highres space.']);
        system([fslPrefix 'flirt -in ' parcs{p} ' -ref ' highres ' -applyxfm -init ' ...
            freesurfer2Highres ' -out ' parcsHR{p} ' -interp nearestneighbour']);
    end

    for m = 1:length(masks)
        fprintf('%s\n',['Transforming ' maskNames{m} ' to highres space.']);
        system([fslPrefix 'flirt -in ' masks{m} ' -ref ' highres ' -applyxfm -init ' ...
            freesurfer2Highres ' -out ' masksHR{m} ' -interp nearestneighbour']);
    end

end

fprintf('%s\n\n',['Finished registration for ' subject ', ' targetFunc]);
    
end