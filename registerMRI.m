
% registerMRI(study,subjects,varargin)
%
% Computes registration files between various spaces, including the target
% functional space used for data analysis, the high-resolution anatomical 
% image, the same anatomical image in the orientation/resolution used by 
% Freesurfer, and MNI152 standard neurological space. Loops through
% subjects. Should be run after preprocMRI and recon-all. Steps include:
% - Register target functional to anatomical using bbregister
% - Register anatomical to MNI space using FNIRT, and to Freesurfer
%   anatomical space using FLIRT
% - Concatenate registrations between multiple spaces where appropriate
% - Generates gray matter, white matter, and cerebrospinal fluid masks from
%   Freesurfer parcellations. Converts Freesurfer-based images (masks and
%   parcellations) to anatomical and target spaces. Masks anatomical using 
%   Freesurfer brain mask.
%
% Example usage: registerMRI('studyName','SUB*','targetName','target_func_expt1')
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
% - study (string): name of the study to analyze, used to define the study 
%       directory as [getenv('FMRI_BASE_DIR') '/' study].
% - subjects (string or cell array of strings): list of subjects to
%       analyze. Can use asterisk-based regular expression. Examples:
%       {'SUB01','SUB02'} or 'SUB*'.
% 
% Variable arguments:
% - overwrite (boolean; default=0): whether to overwrite files that have
%       already been written by this function.
% - targetName (string; default='target_func'): name of image in func 
%       directory that was used for registration by preprocMRI, without 
%       file extension. If a custom image was specified for preprocMRI,
%       the same image should be used here. registerMRI can be run multiple
%       times for different target_func images.
% - highresName (string): name of anatomical image in
%       anat directory to register to, without file extension.  Should be
%       the same image used for surface reconstruction (recon-all).
%       Defaults to the first anatomical image listed in scanlogs.

function registerMRI(study,subjects,varargin)

addpath([strrep(mfilename('fullpath'),mfilename,'') '/utils']);

% Load/check config variables.
[configError, studyDir, fslPrefix] = checkConfig(study);
if ~isempty(configError)
    fprintf('%s\n',configError);
    return;
end;

overwrite = 0;
targetName = 'target_func';     % Functional-resolution image in func directory to register to
highresName = '';               % Anatomical-resolution image in anat directory to register to (e.g. "anat-2")
                                % If not specified, takes first image listed in scanlogs with "anat" label

fsDir = [studyDir '/freesurfer'];
fsPrefix = ['export SUBJECTS_DIR=' fsDir '; '];
parcNames = {'aparc+aseg','aparc.a2009s+aseg'}; % Freesurfer parcellation names.  Note: script assumes #2 is aparc.a2009s+aseg.

% Freesurfer labels for CSF and WM (aparc.a2009s+aseg parcellation)
csfLabelNums = {4:5,14:15,43:44,72,75:76};
wmLabelNums = {2,7,41,46,77:79,85,100,110,155:158,250:255};

fslDir = getenv('FSL_DIR');
standard = [fslDir '/data/standard/MNI152_T1_2mm_brain.nii.gz'];
standardHead = [fslDir '/data/standard/MNI152_T1_2mm.nii.gz'];
standardMask = [fslDir '/data/standard/MNI152_T1_2mm_brain_mask_dil.nii.gz'];

% Convert list of subjects from "dir" inputs to actual names
subjectsNew = {};
if ~iscell(subjects), subjects = {subjects}; end;
for s=1:length(subjects)
    if ~ischar(subjects{s})
        fprintf('%s\n','ERROR: subject argument must be a string or cell array of strings.');
        return;
    end;
    if strfind(subjects{s},'*')
        subjectsTmp = dir([studyDir '/' subjects{s}]);
        for j=1:length(subjectsTmp)
            subjectsNew{end+1} = subjectsTmp(j).name;
        end;
    else
        subjectsNew{end+1} = subjects{s};
    end;
end;
subjects = subjectsNew;

% Edit variable arguments.  Note: optInputs checks for proper input.
varArgList = {'overwrite','targetName','highresName'};
for i=1:length(varArgList)
    argVal = optInputs(varargin,varArgList{i});
    if ~isempty(argVal)
        eval([varArgList{i} ' = argVal;']);
    end;
end;

for s = 1:length(subjects)
    
    subject = subjects{s};
    subjDir = [studyDir '/' subject];
    analDir = [subjDir '/analysis'];
    anatDir = [subjDir '/anat'];
    funcDir = [subjDir '/func'];
    regDir = [subjDir '/reg'];
    regTargDir = [regDir '/' targetName];
    roiDir = [subjDir '/roi'];
    roiTargDir = [roiDir '/' targetName];
    
    scanlogFiles = regExpDir(subjDir,'.*scanlo.*([^~]$)');
    
    if exist(subjDir,'dir')==0
        fprintf('%s\n\n',['Subject directory does not exist for ' subject '. Skipping this subject.']);
        continue;
    elseif exist(analDir,'dir')==0
        fprintf('%s\n\n',['Analysis directory does not exist for ' subject '. Skipping this subject.']);
        continue;
    elseif exist(anatDir,'dir')==0
        fprintf('%s\n\n',['Anatomical directory does not exist for ' subject '. Skipping this subject.']);
        continue;
    elseif exist(funcDir,'dir')==0
        fprintf('%s\n\n',['Functional directory does not exist for ' subject '. Skipping this subject.']);
        continue;
    elseif isempty(scanlogFiles)
        fprintf('%s\n\n',['Missing scanlog file for subject ' subject '. Skipping this subject.']);
        continue;
    end;
    
    if ~exist(regDir,'dir'), mkdir(regDir); end;
    if ~exist(roiTargDir,'dir'), mkdir(roiTargDir); end;
    
    % Find anatomical image to be used for registration
    if ~isempty(highresName)
        if ~exist([anatDir '/' highresName '.nii.gz'])
            fprintf('%s\n\n',['Requested anatomical image does not exist for ' subject '. Skipping this subject.']);
            continue;
        end;
    else
        br = 0;
        for scan = 1:length(scanlogFiles)
            [acqs, expts, runs] = textread([subjDir '/' scanlogFiles(scan).name],'%d%s%d\n');
            for i = 1:length(acqs)
                runSuffix = ['-' int2str(runs(i))];
                if sum(regexpi(expts{i},'anat'))>0
                    if isempty(highresName), highresName = [expts{i} runSuffix]; end;
                    br = 1;
                end;
                if br, break; end;
            end;
            if br, break; end;
        end;
        if isempty(highresName)
            fprintf('%s\n\n',['No anatomical file found in scanlogs for ' subject '. Skipping this subject.']);
            continue;
        end;
    end;
    
    highres = [anatDir '/' highresName '_brain.nii.gz'];
    highresHead = [anatDir '/' highresName '.nii.gz'];
    highres2mm = [anatDir '/' highresName '_brain_2mm.nii.gz'];
    
    target = [funcDir '/' targetName '.nii.gz'];
    
    if ~exist(highresHead,'file')
        fprintf('%s\n\n',['Missing anatomical image for subject ' subject '. Skipping this subject.']);
        continue;
    elseif ~exist(target,'file')
        fprintf('%s\n\n',['Missing target functional image for subject ' subject '. Skipping this subject.']);
        continue;
    end;
    
    % FS mri files
    fsSubjDir = [fsDir '/' subject];
    orig = [fsSubjDir '/mri/orig.nii.gz'];
    origMGZ = strrep(orig,'.nii.gz','.mgz');
    norm = [fsSubjDir '/mri/norm.nii.gz'];
    normMGZ = strrep(norm,'.nii.gz','.mgz');
    fsMask = [fsSubjDir '/mri/brainmask.nii.gz'];
    fsMaskMGZ = strrep(fsMask,'.nii.gz','.mgz');
    fsMaskHR = [anatDir '/fsmask.nii.gz'];
    fsMaskTarg = [roiTargDir '/fsmask.nii.gz'];
    % Parcellations
    for p = 1:length(parcNames)
        parcs{p} = [fsSubjDir '/mri/' parcNames{p} '.nii.gz'];
        parcsMGZ{p} = strrep(parcs{p},'.nii.gz','.mgz');
        parcsHR{p} = [anatDir '/' parcNames{p} '.nii.gz'];
        parcsTarg{p} = [roiTargDir '/' parcNames{p} '.nii.gz'];
    end;
    % GM, WM, CSF masks
    maskNames = {'gmmask','wmmask','csfmask'};
    for m=1:3
        masks{m} = [fsSubjDir '/mri/' maskNames{m} '.nii.gz'];
        masksTarg{m} = [roiTargDir '/' maskNames{m} '.nii.gz'];
        masksHR{m} = [anatDir '/' maskNames{m} '.nii.gz'];
        masksTargErode{m} = [roiTargDir '/' maskNames{m} '_erode1.nii.gz'];
    end;
    
    % Orig <> highres registrations
    orig2Highres = [regDir '/orig2highres.mat'];
    highres2Orig = [regDir '/highres2orig.mat'];
    highres2OrigFS = [regDir '/highres2orig.dat'];
    
    % Anatomical <> standard registrations
    highres2Std = [regDir '/highres2standard.mat'];
    std2Highres = [regDir '/standard2highres.mat'];
    highres2StdImg = [regDir '/highres2standard.nii.gz'];
    highres2StdWarp = [regDir '/highres2standard_warp.nii.gz'];
    highres2StdJac = [regDir '/highres2standard_jac.nii.gz'];
    std2HighresWarp = [regDir '/standard2highres_warp.nii.gz'];
    
    % Target registrations
    target2Highres = [regTargDir '/target_func2highres.mat'];
    target2Highres_flirt = [regTargDir '/target_func2highres.init.mat'];   % FLIRT-based initializing transform.
    target2OrigInit = [regTargDir '/target_func2orig.init.mat'];
    target2OrigInitFS = [regTargDir '/target_func2orig.init.dat'];
    highres2Target = [regTargDir '/highres2target_func.mat'];
    target2Orig = [regTargDir '/target_func2orig.mat'];
    orig2Target = [regTargDir '/orig2target_func.mat'];
    target2OrigFS = [regTargDir '/target_func2orig.dat'];
    target2HighresImg = [regTargDir '/target_func2highres.nii.gz'];
    target2Highres_flirtImg = [regTargDir '/target_func2highres.init.nii.gz'];
    target2StdImg = [regTargDir '/target_func2standard.nii.gz'];
    
    % Check if Freesurfer recon has run to completion.
    if ~exist(parcsMGZ{1},'file')
        fprintf('%s\n\n',['Freesurfer recon has not completed for subject ' subject '. Skipping this subject.']);
        continue;
    end;
    
    fprintf('%s\n',['Running registration for ' subject '.']);
    
    % Convert MGZ to NII
    if ~exist(norm,'file'), system(['mri_convert ' normMGZ ' ' norm]); end;
    if ~exist(orig,'file'), system(['mri_convert ' origMGZ ' ' orig]); end;
    if ~exist(fsMask,'file'), system(['mri_convert ' fsMaskMGZ ' ' fsMask]); end;
    
    % Generate GM/WM/CSF masks in Freesurfer space.
    if ~exist(masks{end},'file') || overwrite
        
        for p = 1:length(parcNames)
            system(['mri_convert ' parcsMGZ{p} ' ' parcs{p}]);
        end;
        
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
            end;
        end;
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
            end;
        end;
        cmd = [cmd ' -bin ' masks{3} '; rm -rf ' tmpPath '*'];
        system(cmd);
        
    end;
    
    % Highres/orig registrations (FLIRT)
    if ~(exist(highres2Orig,'file') && exist(fsMaskHR,'file')) || overwrite
        % Register orig to highres.
        fprintf('%s\n','Registering orig to highres...');
        system([fslPrefix 'flirt -in ' orig ' -ref ' highresHead ' -omat ' orig2Highres ' -cost normcorr ' ...
            '-dof 6 -searchrx -180 180 -searchry -180 180 -searchrz -180 180']);
        system([fslPrefix 'convert_xfm -omat ' highres2Orig ' -inverse ' orig2Highres]);
        system(['tkregister2 --mov ' highresHead ' --targ ' orig ' --s ' subject ' --fsl ' highres2Orig ' --reg ' highres2OrigFS ' --noedit']);
        
        % Mask anatomical using Freesurfer-based mask.
        fprintf('%s\n','Masking anatomical image using Freesurfer mask...');
        system([fslPrefix 'flirt -in ' fsMask ' -ref ' highresHead ' -applyxfm -init ' ...
            orig2Highres ' -out ' fsMaskHR ' -interp nearestneighbour']);
        system([fslPrefix 'fslmaths ' fsMaskHR ' -bin ' fsMaskHR]);
        system([fslPrefix 'fslmaths ' highresHead ' -mul ' fsMaskHR ' ' highres]);
    end;
    
    % Highres/standard registrations (FNIRT)
    if ~(exist(highres2StdWarp,'file') && exist(std2HighresWarp,'file')) || overwrite
        
        % Anatomical to MNI registration using FLIRT/FNIRT.
        fprintf('%s\n','Registering anatomical image to MNI space...');
        system([fslPrefix 'flirt -ref ' standard ' -in ' highres ' -omat ' highres2Std ...
            ' -cost corratio -dof 12 -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -interp trilinear']);
        system([fslPrefix 'fnirt --in=' highresHead ' --aff=' highres2Std ' --cout=' highres2StdWarp ' --iout=' ...
            highres2StdImg ' --jout=' highres2StdJac ' --config=T1_2_MNI152_2mm --ref=' standardHead ...
            ' --refmask=' standardMask ' --warpres=10,10,10']);
        system([fslPrefix 'convert_xfm -omat ' std2Highres ' -inverse ' highres2Std]);
        
        fprintf('%s\n','Inverting warp field...');
        system([fslPrefix 'fslmaths ' highres ' -subsamp2 ' highres2mm]);
        system([fslPrefix 'invwarp --ref=' highres2mm ' --warp=' highres2StdWarp ' --out=' std2HighresWarp]);
        
    end;
    
    % Target/highres registrations (bbregister)
    if ~exist(target2Highres,'file') || overwrite
        % Find initializing FLIRT registration, and convert to Freesurfer .dat format.
        fprintf('%s\n','Running FLIRT for target2Highres reg...');
        if ~exist(regTargDir,'dir'), mkdir(regTargDir); end;
        system([fslPrefix 'flirt -in ' target ' -ref ' highres ' -out ' target2Highres_flirtImg ...
            ' -omat ' target2Highres_flirt ' -dof 6 -searchrx -10 10 -searchry -50 50 -searchrz -10 10']);
        system([fslPrefix 'convert_xfm -omat ' target2OrigInit ' -concat ' highres2Orig ' ' target2Highres_flirt]);
        system(['tkregister2 --mov ' target ' --targ ' orig ' --s ' subject ' --fsl ' target2OrigInit ' --reg ' target2OrigInitFS ' --noedit']);
                            
        % Perform bbregister
        fprintf('%s\n','Running bbregister for target2Highres reg...');
        system([fsPrefix 'bbregister --s ' subject ' --mov ' target ' --reg ' target2OrigFS ' --bold --init-reg ' target2OrigInitFS]);
        system([fsPrefix 'tkregister2 --mov ' target ' --targ ' orig ' --reg ' target2OrigFS ' --fslregout ' target2Orig ' --noedit']);
        system([fslPrefix 'convert_xfm -omat ' orig2Target ' -inverse ' target2Orig]);
        system([fslPrefix 'convert_xfm -omat ' target2Highres ' -concat ' orig2Highres ' ' target2Orig]);
        system([fslPrefix 'convert_xfm -omat ' highres2Target ' -inverse ' target2Highres]);
        
        % Generate target2Highres image
        system([fslPrefix 'flirt -in ' target ' -ref ' highres ' -applyxfm -init ' ...
            target2Highres ' -out ' target2HighresImg]);
        
        % Generate target2standard
        system([fslPrefix 'applywarp --ref=' standard ' --in=' target ' --out=' ...
            target2StdImg ' --warp=' highres2StdWarp ' --premat=' target2Highres]);
        
        system([fslPrefix 'flirt -in ' fsMask ' -ref ' target ' -applyxfm -init ' ...
            orig2Target ' -out ' fsMaskTarg ' -interp nearestneighbour']);
        system([fslPrefix 'fslmaths ' fsMaskTarg ' -bin ' fsMaskTarg]);
    end;
    
    % Register parcellations and masks to target space.
    if ~exist(masksTarg{end},'file') || overwrite
        
        for p = 1:length(parcNames)
            fprintf('%s\n',['Transforming ' parcNames{p} ' to target space.']);
            system([fslPrefix 'flirt -in ' parcs{p} ' -ref ' target ' -applyxfm -init ' ...
            	orig2Target ' -out ' parcsTarg{p} ' -interp nearestneighbour']);
        end;
        
        for m = 1:length(masks)
            fprintf('%s\n',['Transforming ' maskNames{m} ' to target space.']);
            system([fslPrefix 'flirt -in ' masks{m} ' -ref ' target ' -applyxfm -init ' ...
            	orig2Target ' -out ' masksTarg{m} ' -interp nearestneighbour']);
        end;
        
        % Erode WM/CSF masks by one voxel
        for m=2:3
            system([fslPrefix 'fslmaths ' masksTarg{m} ' -ero ' masksTargErode{m}]);
        end;
        
    end;
    
    % Register parcellations and masks to highres space.
    if ~exist(masksHR{end},'file') || overwrite
        
        for p = 1:length(parcNames)
            fprintf('%s\n',['Transforming ' parcNames{p} ' to highres space.']);
            system([fslPrefix 'flirt -in ' parcs{p} ' -ref ' highres ' -applyxfm -init ' ...
            	orig2Highres ' -out ' parcsHR{p} ' -interp nearestneighbour']);
        end;
        
        for m = 1:length(masks)
            fprintf('%s\n',['Transforming ' maskNames{m} ' to highres space.']);
            system([fslPrefix 'flirt -in ' masks{m} ' -ref ' highres ' -applyxfm -init ' ...
            	orig2Highres ' -out ' masksHR{m} ' -interp nearestneighbour']);
        end;
        
    end;
    
    fprintf('%s\n\n',['Finished registration for ' subject '.']);
    
end;

end