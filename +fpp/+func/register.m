
% fpp.func.register(subjID,funcTemplatePath,inputT1Path,fsSubDir)
% 
% Script to register functional to anatomical images using bbregister, and
% move parcellations/masks to functional template space. Should be run
% after fpp.anat.postproc and fpp.func.defineTemplate.
% 
% Arguments:
% - subjID (string): subject ID
% - funcTemplatePath (string): path to functional template image
% - inputT1Path (string): path to preprocesseed T1w image
% - fsSubDir (string): path to Freesurfer subject directory
% 
% Variable arguments:
% - overwrite (boolean; default=0): whether to overwrite files that have
%       already been written by this function.
% - funcTemplateSpace (string): space of functional template, needed if
%       using non-standard funcTemplatePath that lacks a BIDS space entity

function register(subjID,funcTemplatePath,inputT1Path,fsSubDir,varargin)

% Check system configuration
fpp.util.checkConfig;

% Variable arguments
overwrite = 0;
funcTemplateSpace = '';          % Space of functional template image

% Edit variable arguments.  Note: optInputs checks for proper input.
varArgList = {'overwrite'};
for i=1:length(varArgList)
    argVal = fpp.util.optInputs(varargin,varArgList{i});
    if ~isempty(argVal)
        eval([varArgList{i} ' = argVal;']);
    end
end

% Check that inputs exist
if ~exist(funcTemplatePath,'file')
    error('funcTemplatePath does not exist.')
end
if ~exist(inputT1Path,'file')
    error('inputT1Path does not exist.')
end
if ~exist(fsSubDir,'dir')
    error('fsSubDir does not exist.')
end

% Extract names and directories
[funcPreprocDir,funcTemplateName,~] = fpp.util.fileParts(funcTemplatePath);
if isempty(funcPreprocDir), funcPreprocDir = pwd; end
[anatPreprocDir,inputT1Name,~] = fpp.util.fileParts(inputT1Path);
if isempty(anatPreprocDir), anatPreprocDir = pwd; end

% Extract func template and input T1 key values
funcTemplateSpace = fpp.bids.checkNameValue(funcTemplateName,'space');
if isempty(funcTemplateSpace)
    error('funcTemplateSpace must be specified when funcTemplatePath lacks a BIDS space entity.')
end
funcRes = fpp.bids.checkNameValue(funcTemplateName,'res');
anatRes = fpp.bids.checkNameValue(inputT1Path,'res');

% Check if output exists.
finalOutputPath = fpp.bids.changeName(funcTemplatePath,{'desc','echo'},{'brain',[]},'mask','.nii.gz');
if exist(finalOutputPath,'file') && ~overwrite
    return;
end

% Prefix for Freesurfer commands
fsPrefix = ['export SUBJECTS_DIR=' fsSubDir '/..; '];

% Define paths
inputT1BrainPath = fpp.bids.changeName(inputT1Path,'desc','preprocBrain');
origPath = [fsSubDir '/mri/orig.mgz'];
individual2FsnativeXfm = fpp.bids.changeName(inputT1Path,{'desc','space','res','from','to','mode'},...
    {[],[],[],'individual','fsnative','image'},'xfm','.mat');
fsnative2IndividualXfm = fpp.bids.changeName(individual2FsnativeXfm,{'from','to'},{'fsnative','individual'});
func2IndividualXfm = fpp.bids.changeName(funcTemplatePath,{'desc','space','res','echo','from','to','mode'},...
    {[],[],[],[],funcTemplateSpace,'individual','image'},'xfm','.mat');
individual2FuncXfm = fpp.bids.changeName(func2IndividualXfm,{'from','to'},{'individual',funcTemplateSpace});
func2FsnativeXfm = fpp.bids.changeName(func2IndividualXfm,'to','fsnative');
fsnative2FuncXfm = fpp.bids.changeName(individual2FuncXfm,'from','fsnative');
func2FsnativeXfmDat = strrep(func2FsnativeXfm,'.mat','.dat');
fsnative2FuncXfmDat = strrep(fsnative2FuncXfm,'.mat','.dat');
func2IndividualXfmInit = fpp.bids.changeName(func2IndividualXfm,'desc','FlirtInit');
func2FsnativeXfmInit = fpp.bids.changeName(func2FsnativeXfm,'desc','FlirtInit');
func2FsnativeXfmInitDat = strrep(func2FsnativeXfmInit,'.mat','.dat');
pathsToDelete = {func2IndividualXfmInit,func2FsnativeXfmInit,func2FsnativeXfmInitDat};

% Check if fpp.anat.postproc has been run.
if ~exist(individual2FsnativeXfm)
    error('Must run fpp.anat.postproc before fpp.func.register.');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 1: Register func template to individual
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist(func2IndividualXfm,'file') || overwrite
    fprintf('%s\n',['Step 1, Register func template to individual   - ' subjID]);
    fpp.fsl.flirt(funcTemplatePath,inputT1BrainPath,func2IndividualXfmInit,[],'dof',6);   % Initial FLIRT-based xfm
    fpp.fsl.concatXfm(individual2FsnativeXfm,func2IndividualXfmInit,func2FsnativeXfmInit);
    fpp.util.system([fsPrefix 'tkregister2 --mov ' funcTemplatePath ' --targ ' origPath ' --s sub-' subjID ...
        ' --fsl ' func2FsnativeXfmInit ' --reg ' func2FsnativeXfmInitDat ' --noedit']);     % Convert .mat to .dat
    cmd = [fsPrefix 'bbregister --s sub-' subjID ' --mov ' funcTemplatePath ' --reg ' ...
        func2FsnativeXfmDat ' --bold --init-reg ' func2FsnativeXfmInitDat];                 % Refine xfm with bbregister
    fpp.util.system(cmd);
    fpp.util.system([fsPrefix 'tkregister2 --mov ' funcTemplatePath ' --targ ' origPath ' --reg ' ...
        func2FsnativeXfmDat ' --fslregout ' func2FsnativeXfm ' --noedit']);
    fpp.bids.jsonReconstruct(func2FsnativeXfmInit,func2FsnativeXfm,'xfm');                  % Generate json file for bbreg xfm
    [~,fsVersion] = fpp.util.system('recon-all -version');
    fpp.bids.jsonChangeValue(func2FsnativeXfm,{'Description','Software','SoftwareVersion','CommandLine'},...
        {'Affine registration file generated by bbregister.','bbregister',strtrim(fsVersion),cmd});
    fpp.fsl.invertXfm(func2FsnativeXfm,fsnative2FuncXfm);
    fpp.util.system(['tkregister2 --mov ' origPath ' --targ ' funcTemplatePath ' --s sub-' subjID ...
        ' --fsl ' fsnative2FuncXfm ' --reg ' fsnative2FuncXfmDat ' --noedit']);             % Convert .mat to .dat
    fpp.fsl.concatXfm(fsnative2IndividualXfm,func2FsnativeXfm,func2IndividualXfm);
    fpp.fsl.invertXfm(func2IndividualXfm,individual2FuncXfm);
    fpp.fsl.moveImage(funcTemplatePath,inputT1Path,[anatPreprocDir '/' ...
        fpp.bids.changeName(funcTemplateName,{'space','res'},{'individual',anatRes})],func2IndividualXfm);
    fpp.fsl.moveImage(inputT1BrainPath,funcTemplatePath,[funcPreprocDir '/' fpp.bids.changeName(...
        inputT1Name,{'space','res','desc'},{funcTemplateSpace,funcRes,'preprocBrain'})],individual2FuncXfm);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 2: Move masks/parcs to func space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 2, Move masks/parcs to func space         - ' subjID]);
% Move brain mask to func template space
maskPath = fpp.bids.changeName(inputT1Path,'desc','brainFSdil1','mask','.nii.gz');
maskPathFunc = fpp.bids.changeName(funcTemplatePath,{'desc','echo'},{'brain',[]},'mask','.nii.gz');
fpp.fsl.moveImage(maskPath,funcTemplatePath,maskPathFunc,individual2FuncXfm,'interp','nn');
% Move segment ROIs to func template space, and erode WM/CSF masks
roiNames = {'gm','csf','wm','gmcortical','gmsubcortical'};
for r=1:length(roiNames)
    roiPath = fpp.bids.changeName(inputT1Path,'desc',roiNames{r},'mask','.nii.gz');
    roiPathFunc = fpp.bids.changeName(funcTemplatePath,{'desc','echo'},{roiNames{r},[]},'mask','.nii.gz');
    if exist(roiPath,'file') && (~exist(roiPathFunc,'file') || overwrite)
        fpp.fsl.moveImage(roiPath,funcTemplatePath,roiPathFunc,individual2FuncXfm,'interp','nn');
        if sum(strcmp(roiNames{r},{'csf','wm'}))>0
            fpp.fsl.maths(roiPathFunc,'-ero',fpp.bids.changeName(roiPathFunc,'desc',[roiNames{r} 'ero1']));
        end
    end
end
% Move parcellations to func template space
parcs = {'wmparc','aparc+aseg','aparc.a2009s+aseg','MMP','RSN','Gordon'};
for p=1:length(parcs)
    parcPath = fpp.bids.changeName(inputT1Path,'desc',parcs{p},'dseg','.nii.gz');
    parcPathFunc = fpp.bids.changeName(funcTemplatePath,{'desc','echo'},{parcs{p},[]},'dseg','.nii.gz');
    if exist(parcPath,'file') && (~exist(parcPathFunc,'file') || overwrite)
        fpp.fsl.moveImage(parcPath,funcTemplatePath,parcPathFunc,individual2FuncXfm,'interp','nn');
    end
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
fpp.util.system(['rm -rf ' func2FsnativeXfmDat '.*']);  % Additional files generated by bbregister

end