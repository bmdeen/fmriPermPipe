
% fpp.func.model2ndFeat(inputDirs,varargin)
% 
% Step 2 of a two-step process (modelFeat, model2ndFeat) to perform a
% General Linear Model based analysis of fMRI data, FSL's FEAT/FILM.  Step 
% 2 combines results across runs within a subject (fixed effects analysis).
% Should be run after fpp.func.modelFeat.
%
% Example usage: fpp.func.model2ndFeat({'/pathToAnalyses/sub-01_task-faceloc_run-01_space-individual_modelfeat',...
%   '/pathToAnalyses/sub-01_task-faceloc_run-02_space-individual_modelfeat',...
%   '/pathToAnalyses/sub-01_task-faceloc_run-03_space-individual_modelfeat'})
%
% Arguments:
% - inputDirs (cell array of strings): paths to input modelfeat directories
%       for each run to be included in cross-run stats.
% 
% Variable arguments:
% - overwrite (boolean; default=0): whether to overwrite files that have
%       already been written by this function.
% - outputSuffix (string): suffix for output directory
% - analysisDir (string): analysis output dir will be written in this dir.
%   Must match analysisDir used for fpp.func.modelFeat
%
% Critical outputs:
% - contrast.nii.gz: cross-run contrast image
% - contrastvariance.nii.gz: cross-run contrast variance image
% - tstat.nii.gz: cross-run t-statistic
% - zstat.nii.gz: cross-run z-statistic


function model2ndFeat(inputDirs,varargin)

% Check system configuration
fpp.util.checkConfig;

overwrite = 0;              % Whether to overwrite output
outputSuffix = '';          % New suffix for output dir
analysisDir = '';           % Directory for analysis outputs

% Edit variable arguments.  Note: optInputs checks for proper input.
varArgList = {'overwrite','outputSuffix','analysisDir'};
for i=1:length(varArgList)
    argVal = fpp.util.optInputs(varargin,varArgList{i});
    if ~isempty(argVal)
        eval([varArgList{i} ' = argVal;']);
    end
end

% Remove non-alphanumeric characters from outputSuffix
outputSuffix = regexprep(outputSuffix,'[^a-zA-Z0-9]','');

% Define analysis directory
if isempty(analysisDir)
    [analysisDir,~,~] = fpp.util.fileParts(inputDirs{1});
end

% Define inputSuffix
inputSuffix = fpp.bids.checkNameValue(inputDirs{1},'desc');

% Check if inputs exist, load RegressionData.mat
nRuns = length(inputDirs);
for r=1:nRuns
    if ~exist(inputDirs{r},'file')
        error(['Input analysis directory ' inputDirs{r} ' does not exist.']);
    end
    [~,inputNames{r},~] = fpp.util.fileParts(inputDirs{r});
    inputMat = [inputDirs{r} '/' fpp.bids.changeName(inputNames{r},'desc',inputSuffix,'RegressionData','.mat')];
    regrData{r} = load(inputMat);
end
contrastNames = regrData{1}.contrastNames;
nContrasts = length(contrastNames);
nConds = length(regrData{1}.condNames);

% Define and create output directory
outputName = fpp.bids.changeName(inputNames{1},{'run','desc'},{'',[inputSuffix outputSuffix]},'model2feat','');
outputDir = [analysisDir '/' outputName];
if exist(outputDir,'dir')
    if overwrite
        fpp.util.system(['rm -rf ' outputDir]);
    else
        return;
    end
end
mkdir(outputDir);

% Check whether inputs are CIFTI or NIFTI
paths = dir([inputDirs{1} '/' fpp.bids.changeName(inputNames{1},'desc',...
    [inputSuffix contrastNames{1}],'contrast','') '.*nii*']);
if isempty(paths), error('Could not find contrast maps in first input directory.'); end
[~,~,outputExt] = fpp.util.fileParts(paths(1).name);
if ~ismember(lower(outputExt),{'.nii.gz','.nii','.dscalar.nii'})
    error('Inputs must be a NIFTI or CIFTI dscalar files.');
end
isCifti = 0;
if strcmpi(outputExt,'.dscalar.nii')
    isCifti = 1;
end
if isCifti
    imageType = 'cifti';
else
    imageType = 'volume';
end

% Compute contrast variance for weighting runs
for r=1:nRuns
    V = inv(regrData{r}.X'*regrData{r}.X);
    % Estimated contrast variance, up to error variance term
    regrData{r}.conVarBase = diag(regrData{r}.contrastMat*V(1:nConds,1:nConds)*regrData{r}.contrastMat');
end


% Main loop: combine stats across runs.
for c=1:nContrasts
    
    % Output paths
    outputContrastPath = [outputDir '/' fpp.bids.changeName(outputName,'desc',...
        [inputSuffix outputSuffix contrastNames{c}],'contrast',outputExt)];
    outputContrastVarPath = [outputDir '/' fpp.bids.changeName(outputName,'desc',...
        [inputSuffix outputSuffix contrastNames{c}],'contrastvariance',outputExt)];
    outputTStatPath = [outputDir '/' fpp.bids.changeName(outputName,'desc',...
        [inputSuffix outputSuffix contrastNames{c}],'tstat',outputExt)];
    outputZStatPath = [outputDir '/' fpp.bids.changeName(outputName,'desc',...
        [inputSuffix outputSuffix contrastNames{c}],'zstat',outputExt)];
    outputDOFPath = [outputDir '/' fpp.bids.changeName(outputName,'desc',...
        [inputSuffix outputSuffix],'dof','')];
    
    % Combine contrast maps across runs, with inverse variance weighting
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
    fpp.wb.command([imageType '-math'],[],weightEquation,outputContrastPath,flagText);
    
    % Combine contrast variance maps across runs
    weightEquation = '(';
    flagText = '';
    weightSum = 0;
    for r=1:nRuns
        inputContrastVarPaths{r} = [inputDirs{r} '/' fpp.bids.changeName(inputNames{r},'desc',...
            [inputSuffix contrastNames{c}],'contrastvariance',outputExt)];
        cStr = ['c' int2str(r)];
        flagText = [flagText ' -var ' cStr ' ' inputContrastVarPaths{r}];
        if r==1
            weightEquation = [weightEquation cStr '*' num2str(1/regrData{r}.conVarBase(c)^2)];
        else
            weightEquation = [weightEquation '+' cStr '*' num2str(1/regrData{r}.conVarBase(c)^2)];
        end
        weightSum = weightSum + 1/regrData{r}.conVarBase(c);
    end
    weightSum = weightSum^2;
    weightEquation = [weightEquation ')/' num2str(weightSum)];
    fpp.wb.command([imageType '-math'],[],weightEquation,outputContrastVarPath,flagText);
    
    % Compute t-statistic
    equation = 'con/sqrt(convar)';
    flagText = ['-var con ' outputContrastPath ' -var convar ' outputContrastVarPath];
    fpp.wb.command([imageType '-math'],[],equation,outputTStatPath,flagText);
    
    % Compute total dof
    for r=1:nRuns
        inputDOFPaths{r} = [inputDirs{r} '/' fpp.bids.changeName(inputNames{r},'desc',...
            inputSuffix,'dof','')];
        dof(r) = load(inputDOFPaths{r});
    end
    fid = fopen(outputDOFPath,'w');
    fprintf(fid,'%d',sum(dof));
    fclose(fid);
    
    % Convert t- to z-statistic.
    fpp.util.convertTtoZ(outputTStatPath,outputZStatPath,sum(dof));
end

disp(['Finished FEAT 2nd - ' outputName]);

end