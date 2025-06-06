
% fpp.func.modelArma(inputPath,eventsPath,contrastMatrixPath,varargin)
% 
% Step 1 of a two-step process (modelArma, model2ndLevel) to perform a
% General Linear Model based analysis of fMRI data using AFNI's 3dREMLfit,
% modeling autocorrelation with an ARMA(1,1) model. Should be run after
% fpp.func.preproc.
% 
% Example usage: modelArma('/pathToData/sub-01_task-faceloc_run-01_space-individual_desc-preproc_bold.nii.gz',...
%   '/pathToData/sub-01_task-faceloc_run-01_events.tsv','/pathToData/task-faceloc_contrastmatrix.tsv')
% 
% Arguments:
%   - inputPath (string): path to input preprocessed data (NIFTI/CIFTI)
%   - eventsPath (string): path to events.tsv file
%   - contrastMatrixPath (string): path to contrastmatrix.tsv file.
% 
% Example contrastmatrix.tsv file:
%   Name    Coef1   Coef2
%   FacesVsObjects  1   -1
%   ObjectsVsFaces  -1  1
% 
% Variable arguments:
%   - overwrite (boolean; default=0): whether to overwrite files that have
%       already been written by this function.
%   - maskPath (string): path to brain mask image to use for analysis
%   - outputSuffix (string): suffix for output directory/file desc field.
%       Can be used to run analyses multiple times with different options.
%   - analysisDir (string): analysis output dir will be written in this dir
%   - condNames (cell array of strings): array of condition names. Needed
%       to determine cond order if events.tsv lacks trial_type_id field
%   - contrastNames (cell array of strings): array of contrast names to use
%   - confoundPath (string): path to confound.tsv file with nuisance
%       regressors
%   - confoundNames (cell array of strings): fields in confound file to use
%   - confoundFilt (boolean vector): which confound regressors to
%       temporally filter
%   - outlierPath (string): path to outliers.tsv specifying time points to
%       exclude
%   - useTedana (boolean): whether to include TEDANA removed components as
%       confound regressors (unless data has "NoTedana" in desc)
%   - tedanaPath (string): sample tedana directory/file name, to determine
%       space and res (if any) to use for tedana input dir. Default: for
%       NIFTI, taken from inputPath; for CIFTI, space-session is assumed.
%   - hrfType (1 or 2; default=1): which hemodynamic response function to
%       use for regressor convolution. 1 = double-gamma (FSL/SPM default 
%       parameters); 2 = gamma (FS-FAST default parameters).
%   - upsampledTR (scalar in (0,Inf); default=.01): regressor sampling rate 
%       (s) used for convolution. Code requires mod(tr/2,upsampledTR)==0.
%   - tempFilt (boolean; default=0): whether to highpass filter regressors.
%       Should be used if data were filtered by preprocMRI.
%   - filtCutoff (scalar in (0,Inf); default=.01): highpass filter cutoff
%       (Hz)
%   - filtOrder (integer in (0,Inf); default=60): Hanning-window FIR filter
%       order.
%   - deleteAfni (boolean; default true): whether to delete directory
%       containing raw 3dREMLfit outputs, after copying relevant files to 
%       modelarma dir
%   - plotResults (boolean; default=0): whether to display result plots
%   - writeResiduals (boolean; default=0): whether to output GLM residuals
% 
% 
% Critical outputs:
% - contrast: contrast ("contrast of parameter estimate") images
% - beta: parameter estimate (beta) images
% - zstat: z-statistic images
% - desc-Regressors_image.png: image showing regressor time series
% - desc-RegressorCorrelations_image.png: image showing correlations
%   between regressors
% - desc-RegressorVarianceRemoved_image.png: image showing the proportion
%   of variance in task regressors that is explained by nuisance regressors

function modelArma(inputPath,eventsPath,contrastMatrixPath,varargin)

% Check system configuration
fpp.util.checkConfig;

% Basic parameters
overwrite = 0;              % Whether to overwrite output
maskPath = '';              % Path to brain mask image
outputSuffix = '';          % Suffix for output dir
condNames = {};             % Condition names (supercedes events.trial_type_id)
contrastNames = {};         % Contrast names
analysisDir = '';           % Base dir for analysis output directories

% Confound file parameters
confoundPath = fpp.bids.changeName(inputPath,'desc',[],'confounds','.tsv'); % Path to confounds.tsv file
confoundNames = {};         % Confound file fields to use
confoundFilt = [];          % Which confound variables to temporally filter (indexed by confoundNames)
outlierPath = fpp.bids.changeName(inputPath,{'space','desc'},{[],[]},'outliers','.tsv');    % Path to outliers.tsv file
useTedana = 1;              % Whether to include TEDANA removed components as confound regressors (unless data has "NoTedana" in desc)
tedanaPath = '';            % Sample TEDANA path

% Modeling parameters
hrfType = 1;                % 1 = double-gamma (FSL/SPM default parameters)
                            % 2 = gamma (FS-FAST default parameters)
                            % 3 = gamma (MION imaging)
upsampledTR = .01;          % Regressor sampling rate (s) used for convolution. Code requires mod(tr/2,upsampledTR)==0
subtractHalfTR = 1;         % Whether to subtract .5*TR to regressor onsets, to account for slice timing correction
                            % This should only be set to zero if STC method corrected to beginning (not middle) of TR

% Regressor filtering parameters
tempFilt = 0;               % Whether to highpass filter regressors (use if data were filtered)
filtCutoff = 1/100;         % Highpass filter Cutoff (Hz)
filtOrder = 60;             % Hanning-window FIR filter order
filtType = '';              % Type of filter (high, low, bandpass)

% Data generation parameters
deleteAfni = 1;             % Whether to delete output afni directory
plotResults = 0;            % Whether to display result plots
writeResiduals = 0;         % Whether to write 4-D residual image



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 1: Update and validate parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Edit variable arguments.  Note: optInputs checks for proper input.
varArgList = {'overwrite','outputSuffix','permuteRest','tempFilt','filtCutoff',...
    'filtOrder','hrfType','upsampledTR','writeResiduals','permIters','plotResults',...
    'randSeed','condNames','confoundPath','confoundNames','confoundFilt','outlierPath',...
    'contrastNames','filtType','analysisDir','useTedana','maskPath','deleteAfni',...
    'tedanaPath'};
for i=1:length(varArgList)
    argVal = fpp.util.optInputs(varargin,varArgList{i});
    if ~isempty(argVal)
        eval([varArgList{i} ' = argVal;']);
    end
end

% Check confoundNames input
if size(confoundNames,1)>size(confoundNames,2), confoundNames = confoundNames'; end
% Check confoundFilt input
if tempFilt && ~isempty(confoundNames) && (isempty(confoundFilt) ||...
        length(confoundFilt)~=length(confoundNames))
    error(['If temporal filtering and confound regression are both being used, '...
        'a valid confoundFilt argument must be provided.']);
end

% Remove non-alphanumeric characters from outputSuffix
outputSuffix = regexprep(outputSuffix,'[^a-zA-Z0-9]','');

% Check input
if ~exist(inputPath,'file')
    error(['Input path ' inputPath ' does not exist.']);
end
[inputDir,inputName,inputExt] = fpp.util.fileParts(inputPath);
if isempty(inputDir), inputDir = pwd; end
outputName = fpp.bids.changeName(inputName,'desc',outputSuffix,'modelarma','');
outputPrefix = strrep(outputName,'_modelarma','');
numVols = fpp.util.checkMRIProperty('vols',inputPath);
tr = fpp.util.checkMRIProperty('tr',inputPath);
if mod(tr/2,upsampledTR)~=0
    error(['TR/2 must be a multiple of upsampledTR - ' outputName]);
end
if ~ismember(lower(inputExt),{'.nii.gz','.nii','.dtseries.nii'})
    error(['inputPath must be a NIFTI or CIFTI dtseries file - ' outputName]);
end
isCifti = 0;
if strcmpi(inputExt,'.dtseries.nii')
    isCifti = 1;
end

% Define event timing and condition names
events = bids.util.tsvread(eventsPath);
if exist(condNames)
    for b=1:length(events.trial_type)
        events.trial_type_id{b} = find(strcmpi(events.trial_type{b},condNames));
        if isempty(events.trial_type_id{b})
            error(['condNames input does not match events TSV file - ' outputName]);
        end
    end
    if size(condNames,1)>size(condNames,2), condNames = condNames'; end
elseif isfield(events,'trial_type_id')
    for c=1:max(events.trial_type_id)
        ind = find(events.trial_type_id==c);
        if isempty(ind)
            error(['All conditions must be included in each run - ' outputName]);
        end
        condNames{c} = events.trial_type{ind(1)};
    end
else
    error(['If events.tsv does not contain a trial_type_id field, condNames must be '...
        'defined to specify condition order - ' outputName]);
end
nConds = length(condNames);
for c=1:nConds  % 3-column event format
    blockInd = find(events.trial_type_id==c);
    regressors3Col{c} = [events.onset(blockInd) events.duration(blockInd) ones(length(blockInd),1)];
end

% Define contrast matrix
contrastMat = [];
contrastData = bids.util.tsvread(contrastMatrixPath);
if ~isempty(contrastNames)
    contrastInd = find(ismember(contrastData.Name,contrastNames));
else
    contrastNames = contrastData.Name';
    contrastInd = 1:length(contrastNames);
end
for c=1:nConds
    contrastMat(:,end+1) = contrastData.(['Coef' int2str(c)])(contrastInd);
end
nContrasts = size(contrastMat,1);

% Define outlier indices
outlierInd = [];
if exist(outlierPath,'file') && ~isempty(bids.util.tsvread(outlierPath))
    outlierInd = fpp.util.readOutlierTSV(outlierPath);
end

% Define nuisance regressors
nuisRegrMat = [];
if exist(confoundPath,'file') && ~isempty(confoundNames)
    confound = bids.util.tsvread(confoundPath);
    for c=1:length(confoundNames)
        nuisRegrMat(:,end+1) = confound.(confoundNames{c});
    end
    
    nuisRegrMat = bsxfun(@minus,nuisRegrMat,mean(nuisRegrMat));
    
    % Filter nuisance regressors
    if tempFilt
        nuisRegrMat(:,confoundFilt==1) = fpp.util.firFilter(...
            nuisRegrMat(:,confoundFilt==1),1/tr,filtCutoff,filtType,filtOrder);
    end
    
    % Redefine confound TSV after temporal filtering.
    confoundNew = struct();
    for c=1:length(confoundNames)
        confoundNew.(confoundNames{c}) = nuisRegrMat(:,c);
    end
    confound = confoundNew;
else
    confound = struct();
end
regrNames = [condNames,confoundNames];

% Load and add TEDANA regressors to confound struct
% NOTE: Not designed to handle tedana dirs with desc field
inputDesc = fpp.bids.checkNameValue(inputName,'desc');
if useTedana && sum(regexpi(inputDesc,'NoTedana'))==0
    if ~isempty(tedanaPath)
        tedanaSpace = fpp.bids.checkNameValue(tedanaPath,'space');
        tedanaRes = fpp.bids.checkNameValue(tedanaPath,'res');
    else
        if isCifti  % Can't use CIFTI space info, b/c TEDANA is volumetric
            tedanaSpace = 'session';
            tedanaRes = [];
        else
            tedanaSpace = fpp.bids.checkNameValue(inputPath,'space');
            tedanaRes = fpp.bids.checkNameValue(inputPath,'res');
        end
    end
    tedanaTSVPath = [inputDir '/' fpp.bids.changeName(inputName,{'desc','space','res','den'},...
        {'',tedanaSpace,tedanaRes,[]},'tedana','') '/' fpp.bids.changeName(inputName,...
        {'desc','space','res','den'},{'tedanaICARejected', tedanaSpace,tedanaRes,[]},'mixing','.tsv')];
    if exist(tedanaTSVPath,'file')
        tedana = bids.util.tsvread(tedanaTSVPath);
        tedanaNames = fieldnames(tedana)';
        for i=1:length(tedanaNames)
            confound.(tedanaNames{i}) = tedana.(tedanaNames{i});
            nuisRegrMat(:,end+1) = tedana.(tedanaNames{i});
        end
        regrNames = [regrNames,tedanaNames];
    else
        warning(['Couldn''t find tedana path ' tedanaTSVPath ' - tedana regressors not included']);
    end
end

% Demean nuisance regressors, after adding tedana regressors
nuisRegrMat = bsxfun(@minus,nuisRegrMat,mean(nuisRegrMat));

% Define and create output directory
if isempty(analysisDir)
    analysisDir = strrep(inputDir,'/func','/analysis'); % Note: using inputDir/../analysis doesn't work for AFNI
end
outputDir = [analysisDir '/' outputName];
outputAfniDir = [outputDir '/afni'];
if ~exist(analysisDir,'dir'), mkdir(analysisDir); end
if exist(outputDir,'dir')
    if overwrite
        fpp.util.system(['rm -rf ' outputDir]);
    else
        return;
    end
end
mkdir(outputAfniDir);
curDir = pwd;
cd(outputAfniDir);
outputMat = [outputDir '/' fpp.bids.changeName(inputName,'desc',outputSuffix,'RegressionData','.mat')];

% Convert CIFTI to NIFTI for subsequent processing
if isCifti
    inputPathOrig = inputPath;
    inputPath = [outputDir '/' inputName '.nii.gz'];
    fpp.wb.command('cifti-convert',[],[],[],['-to-nifti ' inputPathOrig ' ' inputPath]);
    if ~isempty(maskPath)
        maskPathOrig = maskPath;
        [~,maskName,~] = fpp.util.fileParts(maskPath);
        maskPath = [outputDir '/' maskName '.nii.gz'];
        fpp.wb.command('cifti-convert',[],[],[],['-to-nifti ' maskPathOrig ' ' maskPath]);
    end
end

% Write full confound TSV file
confoundOutputPath = [outputDir '/' fpp.bids.changeName(outputName,'desc',...
    [outputSuffix 'Arma'],'confounds','.tsv')];
if ~isempty(fieldnames(confound))
    confoundOutputPath = [outputDir '/' fpp.bids.changeName(outputName,'desc',...
        [outputSuffix 'Arma'],'confounds','.tsv')];
    bids.util.tsvwrite(confoundOutputPath,confound);
end

% Define and create condition TSV file
condPath = [outputDir '/' fpp.bids.changeName(inputName,'desc',outputSuffix,'conditions','.tsv')];
for c=1:nConds
    condTSV.cond_names{c} = condNames{c};
end
bids.util.tsvwrite(condPath,condTSV);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 2: Define task regressors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define task regressors by convolving boxcar with double-gamma HRF.
taskRegrMat = [];
for r=1:nConds
    taskRegrMat(:,r)= fpp.func.analysis.constructRegressor(regressors3Col{r},...
        hrfType,numVols,tr,upsampledTR,subtractHalfTR);
end
taskRegrMat = bsxfun(@minus,taskRegrMat,mean(taskRegrMat)); % Demean

% Filter task regressors
if tempFilt
    taskRegrMat = fpp.util.firFilter(taskRegrMat,1/tr,filtCutoff,filtType,filtOrder);
end

% Check design matrix X
X = [taskRegrMat nuisRegrMat];
nRegrs = size(X,2);
dof = numVols - nRegrs - 1 - length(outlierInd);    % Regression degrees of freedom
if rank(X)<nRegrs
    error(['Design matrix is rank-deficient - ' outputName]);
end

% Write regressors to .1D files
for r=1:nRegrs
    regrPaths{r} = [outputDir '/' fpp.bids.changeName(outputName,'desc',...
        [outputSuffix regrNames{r}],'regressor','.1D')];
    fid = fopen(regrPaths{r},'w');
    fprintf(fid,'%f\n',X(:,r));
    fclose(fid);
end

% Plot regressor information.
fpp.func.analysis.plotRegressors(taskRegrMat,nuisRegrMat,regrNames,outputDir,...
    outputName,~plotResults);

% Save outputs
save(outputMat,'regressors3Col','X','condNames','tr','inputPath',...
    'confoundPath','outlierPath','hrfType','tempFilt','filtCutoff',...
    'filtOrder','contrastNames','contrastMat');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 3: Run 3dDeconvolve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3dDeconvolve command
designPrefix = [outputAfniDir '/' outputPrefix '_design'];
cmd = ['3dDeconvolve -nodata ' int2str(numVols) ' ' int2str(tr)...
    ' -x1D_stop -polort 0 -bucket ' designPrefix];

% Define regressor info string for 3dDeconvolve command
regrStr = [' -num_stimts ' int2str(nRegrs)];
for r=1:nRegrs
    regrStr = [regrStr ' -stim_file ' int2str(r) ' ' regrPaths{r} ...
        ' -stim_label ' int2str(r) ' ' regrNames{r}];
end
cmd = [cmd regrStr];

% Define contrast info string for 3dDeconvolve command
conStr = [' -num_glt ' int2str(nContrasts)];
for con=1:nContrasts
    conStr = [conStr ' -gltsym ''SYM:'];
    for c=1:nConds
        if contrastMat(con,c)~=0
            conStr = [conStr ' ' int2str(contrastMat(con,c)) '*' condNames{c}];
        end
    end
    conStr = [conStr ''' -glt_label ' int2str(con) ' ' contrastNames{con}];
end
cmd = [cmd conStr];

% Run 3dDeconvolve
disp(['Running 3dDeconvolve - ' outputName]);
fpp.util.system(cmd);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 4: Run 3dREMLfit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

designPath = [designPrefix '.xmat.1D'];
remlPrefix = [outputAfniDir '/' outputPrefix '_remlbuck'];
remlVarPrefix = [outputAfniDir '/' outputPrefix '_remlvar'];
cmd = ['3dREMLfit -matrix ' designPath ' -input ' inputPath...
    ' -tout -noFDR -nobout -Rbuck ' remlPrefix ' -Rvar ' remlVarPrefix];
if writeResiduals
    residPrefix = [outputAfniDir '/' outputPrefix '_residuals'];
    cmd = [cmd ' -Rwherr ' residPrefix];
end
if ~isempty(maskPath)
    cmd = [cmd ' -mask ' maskPath];
end

% Run 3dREMLfit
disp(['Running 3dREMLfit    - ' outputName]);
fpp.util.system(cmd);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 5: Copy/rename output files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
remlPath = [remlPrefix '.nii.gz'];
remlVarPath = [remlVarPrefix '.nii.gz'];

% Determine AFNI space
afniSpaceOptions = {'orig','acpc','tlrc'};
for i=1:length(afniSpaceOptions)
    if exist([remlPrefix '+' afniSpaceOptions{i} '.BRIK'],'file')
        afniSpace = afniSpaceOptions{i};
    end
end

% Convert REMLfit outputs from AFNI to .nii.gz
fpp.util.system(['3dAFNItoNIFTI ' remlPrefix '+' afniSpace '.BRIK']);
fpp.fs.mriConvert(remlPath(1:end-3),remlPath);
fpp.util.system(['rm -rf ' remlPath(1:end-3)]);
fpp.util.system(['3dAFNItoNIFTI ' remlVarPrefix '+' afniSpace '.BRIK']);
fpp.fs.mriConvert(remlVarPath(1:end-3),remlVarPath);
fpp.util.system(['rm -rf ' remlVarPath(1:end-3)]);
% Move remlvar file to outputDir
remlVarNewPath = [outputDir '/' outputPrefix '_remlvar.nii.gz'];
fpp.util.system(['mv ' remlVarPath ' ' remlVarNewPath]);
remlVarPath = remlVarNewPath;
outputMapPaths = {};

% Write mean functional
meanPath = [outputDir '/' fpp.bids.changeName(inputName,'desc',...
    [outputSuffix 'Mean'],'bold','.nii.gz')];
fpp.wb.command('volume-reduce',inputPath,'MEAN',meanPath);
if exist(fpp.bids.jsonPath(meanPath))   % Not yet dealing with json files for analysis outputs!
    fpp.util.system(['rm -rf ' fpp.bids.jsonPath(meanPath)]);
end

for c=1:nConds
    % Beta values
    outputBetaPath = [outputDir '/' fpp.bids.changeName(inputName,'desc',...
        [outputSuffix condNames{c}],'beta','.nii.gz')];
    outputPSCPath = fpp.bids.changeName(outputBetaPath,[],[],'psc');
    fpp.util.system(['fslroi ' remlPath ' ' outputBetaPath ' ' int2str(1+(c-1)*2) ' 1']);
    fpp.wb.command('volume-math',[],'100*beta/mean',outputPSCPath,...
        ['-var beta ' outputBetaPath ' -var mean ' meanPath]);
    outputMapPaths = [outputMapPaths outputBetaPath outputPSCPath];
end
for con=1:nContrasts
    % Contrast values
    outputConPath = [outputDir '/' fpp.bids.changeName(inputName,'desc',...
        [outputSuffix contrastNames{con}],'contrast','.nii.gz')];
    % Contrast std dev
    outputConStdDevPath = [outputDir '/' fpp.bids.changeName(inputName,'desc',...
        [outputSuffix contrastNames{con}],'contraststddev','.nii.gz')];
    % Contrast variance
    outputConVarPath = [outputDir '/' fpp.bids.changeName(inputName,'desc',...
        [outputSuffix contrastNames{con}],'contrastvariance','.nii.gz')];
    % T-stats
    outputTStatPath = [outputDir '/' fpp.bids.changeName(inputName,'desc',...
        [outputSuffix contrastNames{con}],'tstat','.nii.gz')];
    % Z-stats
    outputZStatPath = [outputDir '/' fpp.bids.changeName(inputName,'desc',...
        [outputSuffix contrastNames{con}],'zstat','.nii.gz')];
    
    fpp.util.system(['fslroi ' remlPath ' ' outputConPath ' ' int2str(2*nRegrs+1+(con-1)*2) ' 1']);
    fpp.util.system(['fslroi ' remlPath ' ' outputTStatPath ' ' int2str(2*nRegrs+con*2) ' 1']);
    fpp.fsl.maths(outputConPath,['-div ' outputTStatPath],outputConStdDevPath);
    fpp.fsl.maths(outputConStdDevPath,'-sqr',outputConVarPath);
    fpp.util.system(['rm -rf ' outputConStdDevPath]);
    fpp.util.convertTtoZ(outputTStatPath,outputZStatPath,dof);
    outputMapPaths = [outputMapPaths outputConPath outputConVarPath outputTStatPath outputZStatPath];
end

% DoF
outputDOFPath = [outputDir '/' fpp.bids.changeName(inputName,'desc',outputSuffix,'dof','')];
fid = fopen(outputDOFPath,'w');
fprintf(fid,'%d',dof);
fclose(fid);

% Error variance
outputErrorStdDevPath = [outputDir '/' fpp.bids.changeName(inputName,'desc',outputSuffix,'errorstddev','.nii.gz')];
outputErrorVarPath = [outputDir '/' fpp.bids.changeName(inputName,'desc',outputSuffix,'errorvariance','.nii.gz')];
fpp.util.system(['fslroi ' remlVarPath ' ' outputErrorStdDevPath ' 3 1']);
fpp.fsl.maths(outputErrorStdDevPath,'-sqr',outputErrorVarPath);
fpp.util.system(['rm -rf ' outputErrorStdDevPath]);
outputMapPaths = [outputMapPaths outputErrorVarPath];

% Autocorrelation
outputACPath = [outputDir '/' fpp.bids.changeName(inputName,'desc',outputSuffix,'autocorrelation','.nii.gz')];
fpp.util.system(['fslroi ' remlVarPath ' ' outputACPath ' 2 1']);
outputMapPaths = [outputMapPaths outputACPath];

% Residuals
if writeResiduals
    residPath = [residPrefix '.nii.gz'];
    outputResidPath = [outputDir '/' fpp.bids.changeName(inputName,...
        'desc',[outputSuffix 'Residuals'],'bold','.nii.gz')];
    fpp.util.system(['3dAFNItoNIFTI ' residPrefix '+' afniSpace '.BRIK ' residPath(1:end-3)]);
    fpp.fs.mriConvert(residPath(1:end-3),residPath);
    fpp.util.system(['rm -rf ' residPath(1:end-3)]);
    fpp.util.system(['mv ' residPath ' ' outputResidPath]);
    
    % Convert NIFTI to CIFTI
    outputCiftiPath = fpp.bids.changeName(outputResidPath,[],[],[],'.dscalar.nii');
    fpp.wb.command('cifti-convert',[],[],[],['-from-nifti ' outputResidPath ' '...
        inputPathOrig ' ' outputCiftiPath]);
end

% Design matrix
fpp.util.system(['mv ' designPath ' ' outputDir '/' outputPrefix '_design.xmat.1D']);

% Convert NIFTI to CIFTI
if isCifti
    % Generate mean func CIFTI dscalar file - template for output stat CIFTIs
    inputMeanPath = fpp.bids.changeName(inputPath,'desc','tmpAfniMean',[],'.dscalar.nii');
    fpp.wb.command('cifti-reduce',inputPathOrig,'MEAN',inputMeanPath);
    
    % Generate multi-volume CIFTI dscalar file - template for remlvar
    inputMeanRepeatPath = fpp.bids.changeName(inputMeanPath,'desc','tmpAfniMeanRepeat');
    remlVarVols = fpp.util.checkMRIProperty('Vols',remlVarPath);
    repeatCmd = '';
    for i=1:remlVarVols
        repeatCmd = [repeatCmd '-cifti ' inputMeanPath ' '];
    end
    fpp.wb.command('cifti-merge',[],[],inputMeanRepeatPath,repeatCmd);
    
    % Convert NIFTI to CIFTI
    for i=1:length(outputMapPaths)
        outputCiftiPath = fpp.bids.changeName(outputMapPaths{i},[],[],[],'.dscalar.nii');
        fpp.wb.command('cifti-convert',[],[],[],['-from-nifti ' outputMapPaths{i}...
            ' ' inputMeanPath ' ' outputCiftiPath]);
        if sum(strcmp(outputCiftiPath,{'zstat.dscalar.nii','tstat.dscalar.nii'}))>0
            fpp.wb.command('cifti-palette',outputCiftiPath,'MODE_USER_SCALE',outputCiftiPath,...
                '-pos-user 2.5 6 -neg-user -2.5 -6 -palette-name FSL -disp-pos true');
        end
    end
    outputCiftiPath = fpp.bids.changeName(remlVarPath,[],[],[],'.dscalar.nii');
    fpp.wb.command('cifti-convert',[],[],[],['-from-nifti ' remlVarPath ' '...
        inputMeanRepeatPath ' ' outputCiftiPath]);
    
    % Delete NIFTI and template CIFTI files
    fpp.util.deleteImageAndJson({inputMeanPath,inputMeanRepeatPath});
    fpp.util.system(['rm -rf ' outputDir '/*.nii.gz']);
end

% Delete AFNI directory
cd(curDir);
if deleteAfni, fpp.util.system(['rm -rf ' outputAfniDir]); end

end