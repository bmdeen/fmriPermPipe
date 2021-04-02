
% fpp.func.model2ndPerm(inputPaths,varargin)
% 
% Step 2 of a two-step process (modelPerm, model2ndPerm) to perform a
% General Linear Model based analysis of fMRI data, computing statistics
% using a permutation test.  Step 2 combines results across runs within a
% subject, and computes voxelwise permutation-based statistics. Should be
% run after fpp.func.modelPerm.
%
% NOTE: CIFTI functionality has not been added. Currently only works with
% volumetric (NIFTI) inputs.
%
% Example usage: fpp.func.model2ndPerm({'/pathToData/sub-01_task-faceloc_run-01_space-individual_desc-preproc_bold.nii.gz',...
%   '/pathToData/sub-01_task-faceloc_run-02_space-individual_desc-preproc_bold.nii.gz',...
%   '/pathToData/sub-01_task-faceloc_run-03_space-individual_desc-preproc_bold.nii.gz'})
%
% Arguments:
% - inputPaths (cell array of strings): paths to preprocessed data from
%       each run to be included in cross-run stats.
% 
% Variable arguments:
% - overwrite (boolean; default=0): whether to overwrite files that have
%       already been written by this function.
% - inputSuffix (string): suffix for 1st-level analysis directories
% - outputSuffix (string): suffix for output directory
% - analysisDir (string): analysis output dir will be written in this dir.
%   Must match analysisDir used for fpp.func.modelPerm
% - weightRuns (boolean, default 0): Whether to use inverse variance
%   weighting across runs. This is statistically optimal, and important
%       when runs are different lengths or have different amounts of
%       collinearity between task and nuisance regressors, or outlier
%       volumes removed. However, it roughly doubles compute time, an
%       isn't critical when runs are similar in design, length, and noise
%       properties.
%
% Critical outputs:
% - contrast.nii.gz: cross-run contrast image
% - zstat.nii.gz: cross-run permutation-based voxelwise statistics,
%   transformed to z-scores
% - perms: cross-run contrast maps for each permutation
%

function model2ndPerm(inputPaths,varargin)

% TODO
% - Add CIFTI functionality. Will need to switch all fslmaths commands to
%   wb_command for CIFTI files. Maybe generate function for weighted
%   averaging of NIFTI/CIFTI files?

% Check system configuration
fpp.util.checkConfig;

overwrite = 0;              % Whether to overwrite output
inputSuffix = '';           % Suffix of directories 
outputSuffix = '';          % New suffix for output dir
analysisDir = '';           % Directory for analysis outputs
weightRuns = 0;             % Whether to weight runs by inverse variance factor

charLimit = 100000;         % Character limit for bash commands
outputExt = '.nii.gz';      % Output file extension

% Edit variable arguments.  Note: optInputs checks for proper input.
varArgList = {'overwrite','inputSuffix','outputSuffix','analysisDir','weightRuns'};
for i=1:length(varArgList)
    argVal = fpp.util.optInputs(varargin,varArgList{i});
    if ~isempty(argVal)
        eval([varArgList{i} ' = argVal;']);
    end
end

% Remove non-alphanumeric characters from suffices
outputSuffix = regexprep(outputSuffix,'[^a-zA-Z0-9]','');
inputSuffix = regexprep(inputSuffix,'[^a-zA-Z0-9]','');

% Define analysis directory
[funcPreprocDir,~,~] = fpp.util.fileParts(inputPaths{1});
if isempty(analysisDir)
    analysisDir = [funcPreprocDir '/../analysis'];
end

% Check if input exists, define input modelPerm directories
nRuns = length(inputPaths);
for r=1:nRuns
    if ~exist(inputPaths{r},'file')
        error(['Input path ' inputPaths{r} ' does not exist.']);
    end
    [~,inputNames{r},~] = fpp.util.fileParts(inputPaths{r});
    inputDirs{r} = [analysisDir '/' fpp.bids.changeName(inputNames{r},{'sub','desc'},...
        {[],inputSuffix},'modelperm','')];
    if ~exist(inputDirs{r},'dir')
        error(['Input analysis directory ' inputDirs{r} ' does not exist.']);
    end
    inputMat = [inputDirs{r} '/' fpp.bids.changeName(inputNames{r},'desc',inputSuffix,'RegressionData','.mat')];
    regrData{r} = load(inputMat);
    permItersVec(r) = length(regrData{r}.conVarBasePerm);
end
contrastNames = regrData{1}.contrastNames;
nContrasts = length(contrastNames);
permIters = min(permItersVec);          % # of permutation iterations

% Define and create output directory
outputNameGeneric = fpp.bids.changeName(inputNames{1},{'run','desc'},{'',[inputSuffix outputSuffix]},'model2perm','');
outputName = fpp.bids.changeName(outputNameGeneric,'sub','','model2perm','');
outputDirBase = [analysisDir '/' outputName];
if exist(outputDirBase,'dir')
    if overwrite
        fpp.util.system(['rm -rf ' outputDirBase]);
    else
        return;
    end
end
mkdir(outputDirBase);


% Main loop: compute permutation-based stats
for c=1:nContrasts
    
    outputDir = [outputDirBase '/' contrastNames{c}];
    permsDir = [outputDir '/perms'];
    if ~exist(permsDir,'dir')
        mkdir(permsDir);
    end
    
    % Output files
    outputContrastPath = [outputDir '/' fpp.bids.changeName(outputNameGeneric,'desc',...
        [inputSuffix outputSuffix contrastNames{c}],'contrast',outputExt)];
    outputZStatPath = [outputDir '/' fpp.bids.changeName(outputNameGeneric,'desc',...
        [inputSuffix outputSuffix contrastNames{c}],'zstat',outputExt)];
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% SUBSTEP 1: Average contrast images across runs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Unpermuted data
    if weightRuns
        % Inverse variance weighting across runs
        weightEquation = '(';
        flagText = '';
        weightSum = 0;
        for r=1:nRuns
            inputContrastPaths{r} = [inputDirs{r} '/' fpp.bids.changeName(inputNames{r},'desc',...
                [inputSuffix contrastNames{c}],'contrast',outputExt)];
            cStr = ['c' int2str(r)];
            flagText = [flagText ' -var ' cStr ' ' inputContrastPaths{r}];
            if r==1
                weightEquation = [weightEquation cStr '*' num2str(1/regrData{r}.conVarBase(c))];
            else
                weightEquation = [weightEquation '+' cStr '*' num2str(1/regrData{r}.conVarBase(c))];
            end
            weightSum = weightSum + 1/regrData{r}.conVarBase(c);
        end
        weightEquation = [weightEquation ')/' num2str(weightSum)];
        fpp.wb.command('volume-math',[],weightEquation,outputContrastPath,flagText);
    else
        % Unweighted average
        cmd = '';
        for r=1:nRuns
            inputContrastPaths{r} = [inputDirs{r} '/' fpp.bids.changeName(inputNames{r},'desc',...
                [inputSuffix contrastNames{c}],'contrast',outputExt)];
            if r>1
                cmd = [cmd '-add ' inputContrastPaths{r} ' '];
            end
        end
        cmd = [cmd ' -div ' int2str(nRuns)];
        fpp.fsl.maths(inputContrastPaths{1},cmd,outputContrastPath);
    end
    % Permuted data
    for iter=1:permIters
        iterSuffix = ['iter' int2str(iter)];
        outputContrastPathPerm = [permsDir '/' fpp.bids.changeName(outputNameGeneric,'desc',...
            [iterSuffix inputSuffix outputSuffix contrastNames{c}],'contrast','.nii.gz')];
        if weightRuns
            % Inverse variance weighting across runs
            weightEquation = '(';
            flagText = '';
            weightSum = 0;
            for r=1:nRuns
                inputContrastPaths{r} = [inputDirs{r} '/perms/' iterSuffix '/' fpp.bids.changeName(inputNames{r},'desc',...
                    [iterSuffix inputSuffix contrastNames{c}],'contrast',outputExt)];
                cStr = ['c' int2str(r)];
                flagText = [flagText ' -var ' cStr ' ' inputContrastPaths{r}];
                if r==1
                    weightEquation = [weightEquation cStr '*' num2str(1/regrData{r}.conVarBase(c))];
                else
                    weightEquation = [weightEquation '+' cStr '*' num2str(1/regrData{r}.conVarBase(c))];
                end
                weightSum = weightSum + 1/regrData{r}.conVarBase(c);
            end
            weightEquation = [weightEquation ')/' num2str(weightSum)];
            fpp.wb.command('volume-math',[],weightEquation,outputContrastPathPerm,flagText);
        else
            % Unweighted average
            cmd = '';
            for r=1:nRuns
                inputContrastPaths{r} = [inputDirs{r} '/perms/' iterSuffix '/' fpp.bids.changeName(inputNames{r},'desc',...
                    [iterSuffix inputSuffix contrastNames{c}],'contrast',outputExt)];
                if r>1
                    cmd = [cmd '-add ' inputContrastPaths{r} ' '];
                end
            end
            cmd = [cmd ' -div ' int2str(nRuns)];
            fpp.fsl.maths(inputContrastPaths{1},cmd,outputContrastPathPerm);
        end
        
        if mod(iter,10)==0
            fprintf('%s\n',['Iter ' int2str(iter) ' - ' contrastNames{c} ' ' outputNameGeneric]);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% SUBSTEP 2: Merge permuted contrasts across iterations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    concatContrastPath = [permsDir '/' fpp.bids.changeName(outputNameGeneric,'desc',...
    	[inputSuffix outputSuffix contrastNames{c} 'Permutations'],'contrast','.nii.gz')];
    outputContrastMeanPath = [permsDir '/' fpp.bids.changeName(outputNameGeneric,'desc',...
        [inputSuffix outputSuffix contrastNames{c} 'PermutationMean'],'contrast',outputExt)];
    outputContrastStdDevPath = [permsDir '/' fpp.bids.changeName(outputNameGeneric,'desc',...
        [inputSuffix outputSuffix contrastNames{c} 'PermutationStdDev'],'contraststddev',outputExt)];
    mergeCmd = ['fslmerge -t ' concatContrastPath ' '];
    
    % Split up merge command to avoid Linux command character limits
    cLength = length([permsDir '/' fpp.bids.changeName(outputNameGeneric,'desc',...
            ['iter' permIters inputSuffix outputSuffix contrastNames{c}],'contrast','.nii.gz')]);
    mLength = length(mergeCmd);
    contrastsPerMerge = floor((charLimit-mLength)/cLength);
    mergeIters = ceil(permIters/contrastsPerMerge);
    for i=1:mergeIters
        concatContrastPaths{i} = [permsDir '/' fpp.bids.changeName(outputNameGeneric,'desc',...
            [inputSuffix outputSuffix contrastNames{c} 'Permutations' int2str(i)],'contrast','.nii.gz')];
        mergeCmd = [mergeCmd ' ' concatContrastPaths{i}];
        mergeCmd2 = ['fslmerge -t ' concatContrastPaths{i} ' '];
        for j=(1+(i-1)*contrastsPerMerge):min(permIters,i*contrastsPerMerge)
            outputContrastPathPerm = [permsDir '/' fpp.bids.changeName(outputNameGeneric,'desc',...
                ['iter' int2str(j) inputSuffix outputSuffix contrastNames{c}],'contrast','.nii.gz')];
            mergeCmd2 = [mergeCmd2 ' ' outputContrastPathPerm];
        end
        fpp.util.system(mergeCmd2);
    end
    fpp.util.system(mergeCmd);
    for i=1:mergeIters
        fpp.util.system(['rm -rf ' concatContrastPaths{i}]);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% SUBSTEP 3: Compute permutation-based stats
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute zstat by fitting a voxelwise Gaussian to permuted contrasts
    fpp.util.system(['fslmaths ' concatContrastPath ' -Tmean ' outputContrastMeanPath]);
    fpp.util.system(['fslmaths ' concatContrastPath ' -Tstd ' outputContrastStdDevPath]);
    fpp.util.system(['fslmaths ' outputContrastPath ' -sub ' outputContrastMeanPath ...
        ' -div ' outputContrastStdDevPath ' ' outputZStatPath]);
    
    fprintf('%s\n',[outputNameGeneric ' ' contrastNames{c} ' - Computed permutation stats']);

end

end
