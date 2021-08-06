
% fpp.func.model2ndLevel(inputDirs,varargin)
% 
% Step 2 of a two-step process to perform a General Linear Model based
% analysis of fMRI data for one subject. Step 1 (modelPerm, modelArma, or
% modelFilm) computes results within individual runs. Step 2 combines
% results across runs within an experiment and subject (a "fixed effects"
% analysis). Should be run after fpp.func.modelPerm, Arma, or Film.
% 
% Example usage: fpp.func.model2ndLevel({'/pathToAnalyses/sub-01_task-faceloc_run-01_space-individual_modelarma',...
%   '/pathToAnalyses/sub-01_task-faceloc_run-02_space-individual_modelarma',...
%   '/pathToAnalyses/sub-01_task-faceloc_run-03_space-individual_modelarma'})
% 
% Arguments:
% - inputDirs (cell array of strings): paths to input model directories for
%     	each run to be included in cross-run stats.
% 
% Variable arguments:
% - overwrite (boolean; default=0): whether to overwrite files that have
%       already been written by this function.
% - outputSuffix (string): suffix for output directory/file desc field
% - outputTask (string): task for output directory/files. Should only be
%       specified if results from different experiments are being combined.
% - contrastNames (cell array of strings): array of contrast names to use
% - analysisDir (string): analysis output dir will be written in this dir.
%       Must match analysisDir used for first-level analyses.
% - fdrThresh (vector of values in (0,1)): FDR q-threshold values
% - fdrTails (vector of 1s and 2s): whether FDR thresholding should be
%       one- or two-tailed, for each contrast (default = 2)
% 
% Critical outputs:
% - contrast.nii.gz: cross-run contrast image
% - contrastvariance.nii.gz: cross-run contrast variance image
% - tstat.nii.gz: cross-run t-statistic
% - zstat.nii.gz: cross-run z-statistic


function model2ndLevel(inputDirs,varargin)

% Check system configuration
fpp.util.checkConfig;

overwrite = 0;              % Whether to overwrite output
outputSuffix = '';          % New suffix for output dir
outputTask = '';            % New task name for output dir
contrastNames = {};         % Contrast names
analysisDir = '';           % Directory for analysis outputs
fdrThresh = [.05 .01];      % FDR thresholds
fdrTails = [];              % Whether FDR thresholding for each contrast should be 1- or 2- tailed

% Edit variable arguments.  Note: optInputs checks for proper input.
varArgList = {'overwrite','outputSuffix','analysisDir','fdrThresh','fdrTails','contrastNames','outputTask'};
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

% Determine type of input analysis: perm, arma, or film
modelOptions = {'modelperm','modelarma','modelfilm'};
modelType = '';
for m=1:length(modelOptions)
    if strcmp(inputDirs{1}(end-8:end),modelOptions{m})
        modelType = modelOptions{m}(end-3:end);
        break;
    end
end
if isempty(modelType)
    error('Input directories must have modelperm, modelarma, or modelfilm suffix.');
end

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
customContrast = 0;
if isempty(contrastNames)
    contrastNames = regrData{1}.contrastNames;
    customContrasts = 1;
end
nContrasts = length(contrastNames);
condPath = [inputDirs{1} '/' fpp.bids.changeName(inputNames{1},'desc',inputSuffix,'conditions','.tsv')];
condTSV = bids.util.tsvread(condPath);
nConds = length(condTSV.cond_names);

% Define whether FDR thresholding is 1- or 2-tailed
if isempty(fdrTails)
    fdrTails = 2*ones(1,nContrasts);
elseif length(fdrTails)~=nContrasts
    error('fdrTails must have same length = nContrasts');
end

% Define and create output directory
if isempty(outputTask)
    outputTask = fpp.bids.checkNameValue(inputNames{1},'task');
end
outputName = fpp.bids.changeName(inputNames{1},{'run','desc','task'},{'',...
    [inputSuffix outputSuffix],outputTask},['model2' modelType],'');
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
    % Permute conVarBase indices based on contrastNames varargin
    if customContrast
        for con=1:nContrasts
            if ~ismember(contrastNames{con},regrData{r}.contrastNames)
                error(['Contrast ' contrastNames{con} ' was not found in run '...
                    int2str(r) ' - ' outputName]);
            end
            permInd(con) = find(strcmp(contrastNames{con},regrData{r}.contrastNames));
        end
        regrData{r}.conVarBase = regrData{r}.conVarBase(permInd);
    end
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
    
    if strcmp(modelType,'perm')
        % Compute z-statistic
        equation = 'con/sqrt(convar)';
        flagText = ['-var con ' outputContrastPath ' -var convar ' outputContrastVarPath];
        fpp.wb.command([imageType '-math'],[],equation,outputZStatPath,flagText);
    else
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
        
        if isCifti
            fpp.wb.command('cifti-palette',outputTStatPath,'MODE_USER_SCALE',outputTStatPath,...
                '-pos-user 2.5 6 -neg-user -2.5 -6 -palette-name FSL -disp-pos true');
        end
    end
    if isCifti
        fpp.wb.command('cifti-palette',outputZStatPath,'MODE_USER_SCALE',outputZStatPath,...
            '-pos-user 2.5 6 -neg-user -2.5 -6 -palette-name FSL -disp-pos true');
    end
    
    % FDR threshold
    for f=1:length(fdrThresh)
        outputZStatThreshPath = [outputDir '/' fpp.bids.changeName(outputName,'desc',...
            [inputSuffix outputSuffix contrastNames{c} 'FDR' num2str(fdrThresh(f))],'zstat',outputExt)];
        outputThreshTextPath = [outputDir '/' fpp.bids.changeName(outputName,'desc',...
            [inputSuffix outputSuffix contrastNames{c} 'FDR' num2str(fdrThresh(f))],'fdrthresh','')];
        [critZ,~] = fpp.func.analysis.fdrCorrect(outputZStatPath,outputZStatThreshPath,[],fdrThresh(f),fdrTails(c));
        fid = fopen(outputThreshTextPath,'w+');
        fprintf(fid,'%f',critZ);
        fclose(fid);
    end
end

disp(['Finished model 2nd - ' outputName]);

end