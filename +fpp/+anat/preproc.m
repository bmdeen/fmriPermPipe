%
% fpp.anat.preproc(inputT1Paths,inputT2Paths,outputDir,varargin)
% 
% Preprocesses anatomical (T1- and T2-weighted) images for a single
% participant, prior to Freesurfer surface reconstruction. Register and
% average multiple images if they exist, register image to MNI-like
% ACPC alignment, register T2 to T1 image, and bias-correct.
%
% This script is a modified version of the Human Connectome Project's
% Prefreesurfer pipeline (https://www.humanconnectome.org/software/hcp-mr-pipelines;
% Glasser et al. 2013, "The minimal preprocessing pipelines for the Human 
% Connectome Project").
%
% For hi-res data with T1w and T2w images, follow this script with
% recon-all as follows (opts file should include "mris_inflate -n 30")
%
% recon-all -s ${subject} -i ${anatDir}/sub-${subject}_space-individual_desc-preproc_T1w.nii.gz
%     -all -hires -expert ${optsFile} -T2 ${anatDir}/sub-${subject}_space-individual_desc-preproc_T2w.nii.gz
%     -T2pial -deface
%
%
%
% TODO
% - N4 bias correction option when not using T2 images. (e.g., try to use
%   freesurfer's AntsN4BiasFieldCorrectionFs, to avoid requiring ANTS)
% - Alternative values for standardPath and standardName
% - Change naming of run # for xfms to the real run, if BIDS info is
%   available?
% - Why isn't T2 to T1 xfm concatenation method working?

function preproc(inputT1Paths,inputT2Paths,outputDir,varargin)

% Load/check config variables.
configError = fpp.util.checkConfig;
if ~isempty(configError)
    fprintf('%s\n',configError);
    return;
end

% Basic parameters
overwrite = 0;                  % Whether to overwrite output
[fppFuncDir,~,~]		= fileparts(mfilename('fullpath'));			% path to the directory containing this script
tmp = dir([fppFuncDir '/../../data']);
dataDir = tmp(1).folder;
% standardPath = [getenv('FSLDIR') '/data/standard/MNI152_T1_1mm.nii.gz']; % Alternative: 1mm MNI space
standardDir = [getenv('FSLDIR') '/data/standard'];
standardName = 'MNI152NLin6ASym';
standardPath = [dataDir '/MNI152_T1_0.7mm.nii.gz'];
standardPath2mm = [standardDir '/MNI152_T1_2mm.nii.gz'];
standardMask = [dataDir '/MNI152_T1_2mm_brain_mask_dil_edit.nii.gz'];    % Use dilated input to provide an inclusive initial brain mask

% Edit variable arguments.  Note: optInputs checks for proper input.
varArgList = {'overwrite'};
for i=1:length(varArgList)
    argVal = fpp.util.optInputs(varargin,varArgList{i});
    if ~isempty(argVal)
        eval([varArgList{i} ' = argVal;']);
    end
end

% Check if T2 images are included
if isempty(inputT2Paths)
    usingT2 = 0;
else
    usingT2 = 1;
end

% Define wrapper function for fpp.bids.removeBidsDir (for cellfun functionality)
removeBidsDir = @(x) fpp.bids.removeBidsDir(x);

% Define input/output directories
[inputDir,inputName,~] = fpp.util.fileParts(inputT1Paths{1});
inputNameGeneric = strrep(fpp.bids.changeName(inputName,'run',[]),'_T1w','');
if isempty(inputDir), inputDir = pwd; end
if strcmp(outputDir(end),'/'), outputDir = outputDir(1:end-1); end
anatPreprocDir = [outputDir '/anat'];

% Check if output exists.
finalOutputPath = fpp.bids.changeName([anatPreprocDir '/' inputName '.nii.gz'],{'desc','run','space'},{'preproc',[],'individual'});
if exist(finalOutputPath,'file') && ~overwrite
    return;
end

% Generate output directories
if ~exist(anatPreprocDir,'dir'), mkdir(anatPreprocDir); end

% Copy raw data/metadata, convert to .nii.gz if necessary
for i=1:length(inputT1Paths)
    [~,inputT1Names{i},~] = fpp.util.fileParts(inputT1Paths{i});
    outputT1Paths{i} = [anatPreprocDir '/' inputT1Names{i} '.nii.gz'];
    fpp.util.copyImageAndJson(inputT1Paths{i},outputT1Paths{i},'mri');
    fpp.bids.jsonChangeValue(fpp.bids.jsonPath(outputT1Paths{i}),{'Description','RawSources'},...
        {'Raw data copied to derivative directory.',fpp.bids.removeBidsDir(inputT1Paths{i})});
end
outputT2Paths = {};
for i=1:length(inputT2Paths)
    [~,inputT2Names{i},~] = fpp.util.fileParts(inputT2Paths{i});
    outputT2Paths{i} = [anatPreprocDir '/' inputT2Names{i} '.nii.gz'];
    fpp.util.copyImageAndJson(inputT2Paths{i},outputT2Paths{i},'mri');
    fpp.bids.jsonChangeValue(fpp.bids.jsonPath(outputT2Paths{i}),{'Description','RawSources'},...
        {'Raw data copied to derivative directory.',fpp.bids.removeBidsDir(inputT2Paths{i})});
end
inputT1PathsRaw = outputT1Paths;
inputT2PathsRaw = outputT2Paths;

pathsToDelete = [inputT1PathsRaw inputT2PathsRaw];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 1: Register all anatomicals to first run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 1, Register and average anatomicals       - ' inputNameGeneric]);
inputT1Paths = outputT1Paths;
inputT2Paths = outputT2Paths;
outputT1Path = fpp.bids.changeName(inputT1Paths{1},{'run','space'},{[],'nativeT1w'});
addCmdT1 = '';
if usingT2
    outputT2Path = fpp.bids.changeName(inputT2Paths{1},{'run','space'},{[],'nativeT2w'});
    addCmdT2 = '';
end
for i=2:length(inputT1Paths)
    outputT1Paths{i} = fpp.bids.changeName(inputT1Paths{i},'space','nativeT1w');
    xfmThisRun2FirstRun = fpp.bids.changeName(inputT1Paths{i},{'run','from','to','mode'},...
    {'',['T1wIm' fpp.util.numPad(i,2)],'nativeT1w','image'},'xfm','.mat');
    fpp.fsl.flirt(inputT1Paths{i},inputT1Paths{1},xfmThisRun2FirstRun,outputT1Paths{i},...
        'dof',6,'interp','sinc');
    fpp.bids.jsonChangeValue(fpp.bids.jsonPath(outputT1Paths{i}),{'Description','Sources'},...
        {'Raw data, registered to first anatomical with sinc interpolation.',...
        fpp.bids.removeBidsDir(inputT1Paths{i})});
    addCmdT1 = [addCmdT1 '-add ' outputT1Paths{i} ' '];
    pathsToDelete = [pathsToDelete outputT1Paths{i} xfmThisRun2FirstRun];
end
for i=2:length(inputT2Paths)
    outputT2Paths{i} = fpp.bids.changeName(inputT2Paths{i},'space','nativeT2w');
    xfmThisRun2FirstRun = fpp.bids.changeName(inputT2Paths{i},{'run','from','to','mode'},...
    {'',['T2wIm' fpp.util.numPad(i,2)],'nativeT2w','image'},'xfm','.mat');
    fpp.fsl.flirt(inputT2Paths{i},inputT2Paths{1},xfmThisRun2FirstRun,outputT2Paths{i},...
        'dof',6,'interp','sinc');
    fpp.bids.jsonChangeValue(fpp.bids.jsonPath(outputT2Paths{i}),{'Description','Sources'},...
        {'Raw data, registered to first anatomical with sinc interpolation.',...
        fpp.bids.removeBidsDir(inputT2Paths{i})});
    addCmdT2 = [addCmdT2 '-add ' outputT2Paths{i} ' '];
    pathsToDelete = [pathsToDelete outputT2Paths{i} xfmThisRun2FirstRun];
end
if length(inputT1Paths)>1
    addCmdT1 = [addCmdT1 '-div ' int2str(length(inputT1Paths))];
    fpp.fsl.fslMaths(inputT1Paths{1},addCmdT1,outputT1Path);
    fpp.bids.jsonChangeValue(fpp.bids.jsonPath(outputT1Path),{'Description','Sources','RawSources'},...
        {'Raw data, averaged across runs.',cellfun(removeBidsDir,outputT1Paths,'UniformOutput',false),...
        cellfun(removeBidsDir,inputT1PathsRaw,'UniformOutput',false)});
end
if usingT2 && length(inputT2Paths)>1
    addCmdT2 = [addCmdT2 '-div ' int2str(length(inputT2Paths))];
    fpp.fsl.fslMaths(inputT2Paths{1},addCmdT2,outputT2Path);
    fpp.bids.jsonChangeValue(fpp.bids.jsonPath(outputT2Path),{'Description','Sources','RawSources'},...
        {'Raw data, averaged across runs.',cellfun(removeBidsDir,outputT2Paths,'UniformOutput',false),...
        cellfun(removeBidsDir,inputT2PathsRaw,'UniformOutput',false)});
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 2: Realignment to "std"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 2, Rotate to LAS/RAS                      - ' inputNameGeneric]);
inputT1Path = outputT1Path;
outputT1Path = fpp.bids.changeName(inputT1Path,'space','realignedT1w');
fpp.util.reorientToStd(inputT1Path,outputT1Path);
fpp.bids.jsonChangeValue(fpp.bids.jsonPath(outputT1Path),{'Description','SpatialRef'},...
    {'Raw data, averaged and rotated to LAS/RAS alignment.',fpp.bids.removeBidsDir(outputT1Path)});
pathsToDelete = [pathsToDelete outputT1Path];
if usingT2
    inputT2Path = outputT2Path;
    outputT2Path = fpp.bids.changeName(inputT2Path,'space','realignedT2w');
    fpp.util.reorientToStd(inputT2Path,outputT2Path);
    fpp.bids.jsonChangeValue(fpp.bids.jsonPath(outputT2Path),{'Description','SpatialRef'},...
        {'Raw data, averaged and rotated to LAS/RAS alignment.',fpp.bids.removeBidsDir(outputT2Path)});
    pathsToDelete = [pathsToDelete outputT2Path];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 3: ACPC alignment via rigid-body MNI registration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 3, ACPC alignment                         - ' inputNameGeneric]);
inputT1Path = outputT1Path;
outputT1Path = fpp.bids.changeName(inputT1Path,'space','individual');
inputT1PathRobust = fpp.bids.changeName(inputT1Path,'space','realignedT1wRobustFOV');
inputT12RobustXfm = fpp.bids.changeName(inputT1Path,{'from','to','mode','space'},...
	{'realignedT1w','realignedT1wRobustFOV','image',''},'xfm','.mat');
robust2InputT1Xfm = fpp.bids.changeName(inputT1Path,{'from','to','mode','space'},...
	{'realignedT1wRobustFOV','realignedT1w','image',''},'xfm','.mat');
robust2StandardXfm = fpp.bids.changeName(inputT1Path,{'from','to','mode','space'},...
	{'realignedT1wRobustFOV',standardName,'image',''},'xfm','.mat');
inputT12StandardXfm = fpp.bids.changeName(inputT1Path,{'from','to','mode','space'},...
	{'realignedT1w',standardName,'image',''},'xfm','.mat');
inputT12StandardXfm6Dof = fpp.bids.changeName(inputT1Path,{'from','to','mode','space','desc'},...
	{'realignedT1w',standardName,'image','','6dof'},'xfm','.mat');
% Compute ACPC alignment, based on HCP ACPCAlignment.sh routine
fpp.fsl.robustFOV(inputT1Path,inputT1PathRobust,robust2InputT1Xfm,inputT12RobustXfm);
fpp.fsl.flirt(inputT1PathRobust,standardPath,robust2StandardXfm,[],'dof',12,...
    'searchrx',[-30 30],'searchry',[-30 30],'searchrz',[-30 30]);
fpp.fsl.concatXfm(robust2StandardXfm,inputT12RobustXfm,inputT12StandardXfm);
fpp.fsl.aff2Rigid(inputT12StandardXfm,inputT12StandardXfm6Dof);
fpp.fsl.moveImage(inputT1Path,standardPath,outputT1Path,inputT12StandardXfm6Dof,'interp','sinc');
fpp.bids.jsonChangeValue(outputT1Path,{'Description','SpatialRef'},...
    {'Raw data, averaged and ACPC-aligned.',fpp.bids.removeBidsDir(outputT1Path)});
pathsToDelete = [pathsToDelete inputT1PathRobust inputT12RobustXfm robust2InputT1Xfm...
    robust2StandardXfm inputT12StandardXfm inputT12StandardXfm6Dof];
if usingT2
    inputT2Path = outputT2Path;
    outputT2Path = fpp.bids.changeName(inputT2Path,'space','acpcT2w');
    inputT2PathRobust = fpp.bids.changeName(inputT2Path,'space','realignedT2wRobustFOV');
    inputT2_2RobustXfm = fpp.bids.changeName(inputT2Path,{'from','to','mode','space'},...
        {'realignedT2w','realignedT2wRobustFOV','image',''},'xfm','.mat');
    robust2InputT2Xfm = fpp.bids.changeName(inputT2Path,{'from','to','mode','space'},...
        {'realignedT2wRobustFOV','realignedT2w','image',''},'xfm','.mat');
    robust2StandardXfm = fpp.bids.changeName(inputT2Path,{'from','to','mode','space'},...
        {'realignedT2wRobustFOV',standardName,'image',''},'xfm','.mat');
    inputT2_2StandardXfm = fpp.bids.changeName(inputT2Path,{'from','to','mode','space'},...
        {'realignedT2w',standardName,'image',''},'xfm','.mat');
    inputT2_2StandardXfm6Dof = fpp.bids.changeName(inputT2Path,{'from','to','mode','space','desc'},...
        {'realignedT2w',standardName,'image','','6dof'},'xfm','.mat');
    % Compute ACPC alignment, based on HCP ACPCAlignment.sh routine
    fpp.fsl.robustFOV(inputT2Path,inputT2PathRobust,robust2InputT2Xfm,inputT2_2RobustXfm);
    fpp.fsl.flirt(inputT2PathRobust,standardPath,robust2StandardXfm,[],'dof',12,...
        'searchrx',[-30 30],'searchry',[-30 30],'searchrz',[-30 30]);
    fpp.fsl.concatXfm(robust2StandardXfm,inputT2_2RobustXfm,inputT2_2StandardXfm);
    fpp.fsl.aff2Rigid(inputT2_2StandardXfm,inputT2_2StandardXfm6Dof);
    fpp.fsl.moveImage(inputT2Path,standardPath,outputT2Path,inputT2_2StandardXfm6Dof,'interp','sinc');
    fpp.bids.jsonChangeValue(outputT2Path,{'Description','SpatialRef'},...
        {'Raw data, averaged and ACPC-aligned.',fpp.bids.removeBidsDir(outputT2Path)});
    pathsToDelete = [pathsToDelete inputT2PathRobust inputT2_2RobustXfm robust2InputT2Xfm...
        robust2StandardXfm inputT2_2StandardXfm inputT2_2StandardXfm6Dof outputT2Path];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 3.5: T2 to T1 registration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if usingT2
    fprintf('%s\n',['Step 3.5, T2 to T1 registration                - ' inputNameGeneric]);
    inputT1Path = outputT1Path;
    inputT2Path = outputT2Path;
    outputT2Path = fpp.bids.changeName(inputT2Path,'space','individual');
    inputT2_2InputT1Xfm = fpp.bids.changeName(inputT2Path,{'from','to','mode','space'},...
        {'acpcT2w','individual','image',''},'xfm','.mat');
%     stdT2_2InputT1Xfm = fpp.bids.changeName(inputT2Path,{'from','to','mode','space'},...
%         {'realignedT2w','individual','image',''},'xfm','.mat');
%     stdT2_2InputT2Xfm6Dof = inputT2_2StandardXfm6Dof;
    fpp.fsl.flirt(inputT2Path,inputT1Path,inputT2_2InputT1Xfm,outputT2Path,'dof',6,...
        'interp','sinc','searchrx',[-30 30],'searchry',[-30 30],'searchrz',[-30 30]);
%     fpp.fsl.concatXfm(inputT2_2InputT1Xfm,stdT2_2InputT2Xfm6Dof,stdT2_2InputT1Xfm);
%     fpp.fsl.moveImage(inputT2Path,inputT1Path,outputT2Path,stdT2_2InputT1Xfm);
% NOTE: concat xfm method not working!! Not sure why. (BD, 11/11/20)
    pathsToDelete = [pathsToDelete inputT2_2InputT1Xfm];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 4: Brain masking (rough, FNIRT-based)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 4: Rough brain masking                    - ' inputNameGeneric]);

% Define output paths
inputT1Path = outputT1Path;
outputT1PathStandard = fpp.bids.changeName(inputT1Path,{'space','desc'},{standardName,'nonlinearInit'});
outputMaskPath = fpp.bids.changeName(inputT1Path,{'desc'},{'nonlinearInitBrain'},'mask');
individual2StandardXfmLinear = fpp.bids.changeName(inputT1Path,{'from','to','mode','space','desc'},...
    {'individual',standardName,'image','','linearInit'},'xfm','.mat');
individual2StandardXfm = fpp.bids.changeName(inputT1Path,{'from','to','mode','space','desc'},...
    {'individual',standardName,'image','','nonlinearInit'},'xfm','.nii.gz');
standard2IndividualXfmLinear = fpp.bids.changeName(individual2StandardXfmLinear,{'from','to'},{standardName,'individual'});
standard2IndividualXfm = fpp.bids.changeName(individual2StandardXfm,{'from','to'},{standardName,'individual'});
% Compute initial FLIRT/FNIRT 2 standard transforms
fpp.fsl.flirt(inputT1Path,standardPath2mm,individual2StandardXfmLinear,[],'dof',12);
fpp.fsl.fnirt(inputT1Path,standardPath2mm,'aff',individual2StandardXfmLinear,'refmask',...
    standardMask,'cout',individual2StandardXfm,'iout',outputT1PathStandard,...
    'config',[getenv('FSLDIR') '/etc/flirtsch/T1_2_MNI152_2mm.cnf']);
fpp.fsl.invertXfm(individual2StandardXfmLinear,standard2IndividualXfmLinear);
fpp.fsl.invertWarp(individual2StandardXfm,standard2IndividualXfm,standardPath2mm);  % Using 2mm reference for speed, since images are in same space
% Define rough brain mask, by transforming dilated MNI space mask
fpp.fsl.moveImage(standardMask,inputT1Path,outputMaskPath,[],'warp',...
    standard2IndividualXfm,'rel',1,'interp','nn');
fpp.bids.jsonReconstruct(inputT1Path,outputMaskPath,'mri');
fpp.bids.jsonChangeValue(outputMaskPath,{'Description','Sources','Type'},...
    {'Initial brain mask for anatomical, generated by whole-head FNIRT registration to MNI.',...
    fpp.bids.removeBidsDir(inputT1Path),'Brain'});



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 5: Bias correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 5: Bias correction                        - ' inputNameGeneric]);
if usingT2
    inputT1Path = outputT1Path;
    inputT2Path = outputT2Path;
	outputT1Path = fpp.bids.changeName(inputT1Path,'desc','preproc');
	outputT2Path = fpp.bids.changeName(inputT2Path,'desc','preproc');
    outputBiasPath = fpp.bids.changeName(inputT1Path,[],[],'bias');
    
    % Compute bias field with sqrt(T1w*T2w) method
    fpp.wb.computeBiasField(inputT1Path,inputT2Path,outputMaskPath,outputBiasPath);
    
    % Bias-field correct T1 and T2 images
    fpp.util.copyImageAndJson(inputT1Path,outputT1Path,'mri');
    fpp.bids.jsonChangeValue(fpp.bids.jsonPath(outputT1Path),'Description',...
        'Anatomical data preprocessed by fmriPermPipe.');
    fpp.fsl.fslMaths(inputT1Path,['-div ' outputBiasPath],outputT1Path);
    fpp.util.copyImageAndJson(inputT2Path,outputT2Path,'mri');
    fpp.bids.jsonChangeValue(fpp.bids.jsonPath(outputT2Path),'Description',...
        'Anatomical data preprocessed by fmriPermPipe.');
    fpp.fsl.fslMaths(inputT2Path,['-div ' outputBiasPath],outputT2Path);
    
else
    % TO ADD HERE: N4 BIAS CORRECTION OPTION
    
    inputT1Path = outputT1Path;
	outputT1Path = fpp.bids.changeName(inputT1Path,'desc','preproc');
    fpp.util.copyImageAndJson(inputT1Path,outputT1Path,'mri');
    fpp.bids.jsonChangeValue(fpp.bids.jsonPath(outputT1Path),'Description',...
        'Anatomical data preprocessed by fmriPermPipe.');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CLEANUP: Delete unneeded intermediate images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(pathsToDelete)
    fpp.util.system(['rm -rf ' pathsToDelete{i}]);
    if exist(fpp.bids.jsonPath(pathsToDelete{i}),'file')
        fpp.util.system(['rm -rf ' fpp.bids.jsonPath(pathsToDelete{i})]);
    end
end

end