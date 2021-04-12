
% fpp.func.modelPerm(inputPath,eventsPath,contrastMatrixPath,varargin)
%
% Step 1 of a two-step process (modelPerm, model2ndPerm) to perform a
% General Linear Model based analysis of fMRI data, computing statistics
% using a permutation test.  Step 1 performs voxelwise time series 
% regressions within runs, for a large number of permuted block orderings.
% Note that this step on its own does not output permutation-based 
% statistics, which is done by fpp.func.model2ndPerm, combining results
% across runs. Should be run after fpp.func.preproc.
% 
% Example usage: modelPerm('/pathToData/sub-01_task-faceloc_run-01_space-individual_desc-preproc_bold.nii.gz',...
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
%   - outputSuffix (string): suffix for output directory. Can be used to
%       run the script multiple times with different options.
%   - analysisDir (string): analysis output dir will be written in this dir
%   - condNames (cell array of strings): array of condition names. Needed
%       to determine cond order if events.tsv lacks trial_type_id field
%   - contrastNames (cell array of strings): array of contrast names to use
%   - maskPath (string): path to volumetric mask image (NIFTI)
%   - confoundPath (string): path to confound.tsv file with nuisance
%       regressors
%   - confoundNames (cell array of strings): fields in confound file to use
%   - outlierPath (string): path to outliers.tsv specifying time points to
%       exclude
%   - outlierInd (vector): vector of time points to exclude (supercedes
%       outlierPath)
%   - randSeed (scalar): seed for random number generator
%   - permIters (integer in (0,Inf), default=5000): iterations of
%       permutation test. It is recommended to use at least 5000. Can set
%       to zero to just compute unpermuted OLS betas/contrasts.
%   - permuteRest (boolean; default=0): whether to permute baseline blocks
%       as well as task blocks.  See README (Notes on Theory and Practice
%       section) for further discussion.
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
%   - plotResults (boolean; default=0): whether to display result plots
%   - writeResiduals (boolean; default=0): whether to output GLM residuals
%
%
% Critical outputs:
% - contrast: contrast ("contrast of parameter estimate") images
% - beta: parameter estimate (beta) images
% - desc-OLS_zstat: z-statistics computed based on ordinary least squares.
%   These are not valid statistics, due to the presence of temporal
%   autocorrelation, but can give a rough sense of where effects are
%   located.
% - desc-Regressors_image.png: image showing regressor time series
% - desc-RegressorCorrelations_image.png: image showing correlations
%   between regressors
% - desc-RegressorVarianceRemoved_image.png: image showing the proportion
%   of variance in task regressors that is explained by nuisance regressors
% - perms: directory containing contrast images for each permutation

function modelPerm(inputPath,eventsPath,contrastMatrixPath,varargin)

% Check system configuration
fpp.util.checkConfig;

% Basic parameters
overwrite = 0;              % Whether to overwrite output
outputSuffix = '';          % Suffix for output dir
permIters = 5000;           % Iterations of permutation test
condNames = {};             % Condition names (supercedes events.trial_type_id)
contrastNames = {};         % Contrast names
randSeed = sum(100*clock);  % Random number generator seed
maskPath = '';              % Path to brain mask (NIFTI/CIFTI)
analysisDir = '';           % Base dir for analysis output directories

% Confound file parameters
confoundPath = fpp.bids.changeName(inputPath,'desc',[],'confounds','.tsv'); % Path to confounds.tsv file
confoundNames = {};         % Confound file fields to use
outlierPath = fpp.bids.changeName(inputPath,{'space','desc'},{[],[]},'outliers','.tsv');    % Path to outliers.tsv file
outlierInd = [];            % Outlier time point indices (supercedes outlierPath)

% Modeling parameters
permuteRest = 0;            % Whether to permute baseline blocks as well as task blocks
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
plotResults = 0;            % Whether to display result plots
writeResiduals = 0;         % Whether to write 4-D residual image



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 1: Update and validate parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Edit variable arguments.  Note: optInputs checks for proper input.
varArgList = {'overwrite','outputSuffix','permuteRest','tempFilt','filtCutoff',...
    'filtOrder','hrfType','upsampledTR','writeResiduals','permIters','plotResults',...
    'randSeed','condNames','maskPath','confoundPath','confoundNames','outlierPath',...
    'outlierInd','contrastNames','filtType','analysisDir'};
for i=1:length(varArgList)
    argVal = fpp.util.optInputs(varargin,varArgList{i});
    if ~isempty(argVal)
        eval([varArgList{i} ' = argVal;']);
    end
end

% Update random number generator
if exist('rng','file')
    rng(randSeed,'twister');
else
    rand('twister',randSeed);
end

% Remove non-alphanumeric characters from outputSuffix
outputSuffix = regexprep(outputSuffix,'[^a-zA-Z0-9]','');

% Check input
if ~exist(inputPath,'file')
    error(['Input path ' inputPath ' does not exist.']);
end
[inputDir,inputName,inputExt] = fpp.util.fileParts(inputPath);
if isempty(inputDir), inputDir = pwd; end
outputNameGeneric = fpp.bids.changeName(inputName,'desc',outputSuffix,'bold','');
numVols = fpp.util.checkMRIProperty('vols',inputPath);
tr = fpp.util.checkMRIProperty('tr',inputPath);
exptDuration = tr*numVols;
if mod(tr/2,upsampledTR)~=0
    error(['TR/2 must be a multiple of upsampledTR - ' outputNameGeneric]);
end
if ~ismember(lower(inputExt),{'.nii.gz','.nii','.dtseries.nii'})
    error(['inputPath must be a NIFTI or CIFTI dtseries file - ' outputNameGeneric]);
end
isCifti = 0;
if strcmpi(inputExt,'.dtseries.nii')
    isCifti = 1;
end
if isCifti
    outputExt = '.dscalar.nii';
    outputExtSeries = '.dtseries.nii';
    imageType = 'cifti';
else
    outputExt = '.nii.gz';
    outputExtSeries = '.nii.gz';
    imageType = 'volume';
end

% Define event timing and condition names
events = bids.util.tsvread(eventsPath);
if exist(condNames)
    for b=1:length(events.trial_type)
        events.trial_type_id{b} = find(strcmpi(events.trial_type{b},condNames));
        if isempty(events.trial_type_id{b})
            error(['condNames input does not match events TSV file - ' outputNameGeneric]);
        end
    end
    if size(condNames,1)>size(condNames,2), condNames = condNames'; end
elseif isfield(events,'trial_type_id')
    for c=1:max(events.trial_type_id)
        ind = find(events.trial_type_id==c);
        if isempty(ind)
            error(['All conditions must be included in each run - ' outputNameGeneric]);
        end
        condNames{c} = events.trial_type{ind(1)};
    end
else
    error(['If events.tsv does not contain a trial_type_id field, condNames must be '...
        'defined to specify condition order - ' outputNameGeneric]);
end
nConds = length(condNames);
for c=1:nConds  % 3-column event format
    blockInd = find(events.trial_type_id==c);
    regressors3Col{c} = [events.onset(blockInd) events.duration(blockInd) ones(length(blockInd),1)];
end 
regressors3ColOrig = regressors3Col;

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
    contrastMat(:,end+1) = eval(['contrastData.Coef' int2str(c) '(contrastInd)']);
end
nContrasts = size(contrastMat,1);

% Define outlier indices
if isempty(outlierInd) && exist(outlierPath,'file')
    outlier = bids.util.tsvread(outlierPath);
    if isempty(outlier)
        outlierInd = [];
    else
        outlierNames = fieldnames(outlier);
        for i=1:length(outlierNames)
            outlierInd(end+1) = eval(['find(outlier.' outlierNames{i} ')']);
        end
    end
end
goodVolInd = setdiff(1:numVols,outlierInd);

% Define nuisance regressors
nuisRegrMat = [];
if exist(confoundPath,'file')
    confound = bids.util.tsvread(confoundPath);
    for c=1:length(confoundNames)
        nuisRegrMat(:,end+1) = eval(['confound.' confoundNames{c}]);
    end
    
    if ~isempty(nuisRegrMat)
        nuisRegrMat = nuisRegrMat(goodVolInd,:);
        nuisRegrMat = bsxfun(@minus,nuisRegrMat,mean(nuisRegrMat));
        
        % Filter nuisance regressors
        if tempFilt
            nuisRegrMat = fpp.util.firFilter(nuisRegrMat,1/tr,filtCutoff,filtType,filtOrder);
        end
    end
end
regrNames = [condNames,confoundNames];

% Define and create output directory
outputName = fpp.bids.changeName(inputName,'desc',outputSuffix,'modelperm','');
if isempty(analysisDir)
    analysisDir = [inputDir '/../analysis'];
end
outputDir = [analysisDir '/' outputName];
if ~exist(analysisDir,'dir'), mkdir(analysisDir); end
if exist(outputDir,'dir')
    if overwrite
        fpp.util.system(['rm -rf ' outputDir]);
    else
        return;
    end
end
mkdir(outputDir);
outputDirOrig = outputDir;
outputMat = [outputDirOrig '/' fpp.bids.changeName(inputName,'desc',outputSuffix,'RegressionData','.mat')];

% Define and create condition TSV file
condPath = [outputDir '/' fpp.bids.changeName(inputName,'desc',outputSuffix,'conditions','.tsv')];
for c=1:nConds
    condTSV.cond_names{c} = condNames{c};
end
bids.util.tsvwrite(condPath,condTSV);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 2: Load data, run regression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load mask data
maskVol = [];
if ~isempty(maskPath)
    mask = fpp.util.mriRead(maskPath);
    maskVol = mask.vol;
end

[funcMat,hdr] = fpp.util.readDataMatrix(inputPath,maskVol);     % Read NIFTI/CIFTI input to time by coordinate matrix
funcMat = funcMat';
funcMat = funcMat(goodVolInd,:);                                % Remove artifact time points
funcMat = bsxfun(@minus,funcMat,mean(funcMat));                 % Subtract time series mean

for iter=0:permIters
    
    % Redefine regressors3Col for new iteration
    regressors3Col = regressors3ColOrig;
    
    iterSuffix = '';
    if iter>0
        outputDir = [outputDirOrig '/perms/iter' int2str(iter)];
        mkdir(outputDir);
        iterSuffix = ['iter' int2str(iter)];
    end
    
    % Concatenate all regressors, for block permutation
    allRegr3Col = [];
    for r=1:nConds
        allRegr3Col = [allRegr3Col; regressors3Col{r}];
    end
    
    % Add rest blocks to regressors3Col if we want to permute them
    if permuteRest && iter>0
        regressors3Col{nConds+1} = fpp.func.analysis.defineRestBlocks(allRegr3Col,exptDuration,1);
        allRegr3Col = [allRegr3Col; regressors3Col{end}];
    end
    
    % Permute block order
    if iter>0
        blockPerm{iter} = randperm(size(allRegr3Col,1));
        allRegr3Col = allRegr3Col(blockPerm{iter},:);
        blockCount = 0;
        for r=1:length(regressors3Col)
            regressors3Col{r} = allRegr3Col(blockCount+1:...
                blockCount+size(regressors3Col{r},1),:);
            blockCount = blockCount+size(regressors3Col{r},1);
        end
    end
    if permuteRest && iter>0
        regressors3Col = regressors3Col(1:end-1);
    end
    
    % Define task regressors by convolving boxcar with double-gamma HRF.
    taskRegrMat = [];
    for r=1:nConds
        taskRegrMat(:,r)= fpp.func.analysis.constructRegressor(regressors3Col{r},...
            hrfType,numVols,tr,upsampledTR,subtractHalfTR);
    end
    
    % Remove artifact time points and demean
    taskRegrMat = taskRegrMat(goodVolInd,:);
    taskRegrMat = bsxfun(@minus,taskRegrMat,mean(taskRegrMat));
    
    % Filter task regressors
    if tempFilt
        taskRegrMat = fpp.util.firFilter(taskRegrMat,1/tr,filtCutoff,filtType,filtOrder);
    end
    
    % On first iteration, plot regressor information.
    if iter==0
        fpp.func.analysis.plotRegressors(taskRegrMat,nuisRegrMat,regrNames,outputDir,outputNameGeneric,~plotResults)
    end
    
    % Design matrix X
    X = [taskRegrMat nuisRegrMat];
    if rank(X)<size(X,2)
        error(['Design matrix is rank-deficient - ' outputNameGeneric]);
    end
    
    % Compute beta values and residuals
    V = inv(X'*X);
    betas = V*X'*funcMat;
    if iter==0
        resids = funcMat-X*betas;
        errorDOF = size(X,1)-size(X,2)-1;
        errorVar = sum(resids.^2)/errorDOF;
        rSquareds = ones(1,size(funcMat,2))-sum(resids.^2)./sum(funcMat.^2);
    end
    
    % Compute contrasts and statistics
    contrasts = contrastMat*betas(1:nConds,:);
    conVarBase = diag(contrastMat*V(1:nConds,1:nConds)*contrastMat');  % Estimated contrast variance, up to error variance term
    if iter==0
        conVars = repmat(errorVar,[size(contrasts,1) 1]).*repmat(conVarBase,[1 size(contrasts,2)]);
        tStats = contrasts./sqrt(conVars);
        zStats = zeros(size(tStats));
        zStats(tStats>0) = -icdf('norm',tcdf(-tStats(tStats>0),errorDOF),0,1);
        zStats(tStats<=0) = icdf('norm',tcdf(tStats(tStats<=0),errorDOF),0,1); % Prevents p-values near 1 from rounding to 1
    end
    
    % Save regression information
    if iter==0
        save(outputMat,'regressors3Col','X','V','betas','contrastMat',...
            'contrasts','conVarBase','resids','errorDOF','errorVar',...
            'conVars','tStats','zStats','rSquareds','regrNames','tr',...
            'inputPath','maskPath','confoundPath','outlierPath','outlierInd',...
            'randSeed','permuteRest','hrfType','tempFilt','filtCutoff',...
            'filtOrder','contrastNames');
    else
        conVarBasePerm{iter} = conVarBase;
    end
    
    % Write contrast maps
    for c=1:nContrasts
        outputPath = [outputDir '/' fpp.bids.changeName(inputName,'desc',...
            [iterSuffix outputSuffix contrastNames{c}],'contrast',outputExt)];
        fpp.util.writeDataMatrix(contrasts(c,:)',hdr,outputPath,maskVol);
    end
    
    % For unpermuted analysis, write beta/variance/zstat/etc images
    if iter==0
        % Mean functional
        meanPath = [outputDir '/' fpp.bids.changeName(inputName,'desc',...
            [outputSuffix 'Mean'],'bold',outputExt)];
        fpp.wb.command([imageType '-reduce'],inputPath,'MEAN',meanPath);
        
        % Beta values and percent signal change
        for r=1:size(betas,1)
            betaPath = [outputDir '/' fpp.bids.changeName(inputName,'desc',...
                [outputSuffix condNames{r}],'beta',outputExt)];
            pscPath = fpp.bids.changeName(betaPath,[],[],'psc',outputExt);
            fpp.util.writeDataMatrix(betas(r,:)',hdr,betaPath,maskVol);
            fpp.wb.command([imageType '-math'],[],'100*beta/(mean*mask)',pscPath,...
                ['-var beta ' betaPath ' -var mean ' meanPath ' -var mask ' maskPath]);
        end
        
        % Contrast values
        for c=1:nContrasts
            outputPath = [outputDir '/' fpp.bids.changeName(inputName,'desc',...
                [outputSuffix 'OLS' contrastNames{c}],'zstat',outputExt)];
            fpp.util.writeDataMatrix(zStats(c,:)',hdr,outputPath,maskVol);
            outputPath = [outputDir '/' fpp.bids.changeName(inputName,'desc',...
                [outputSuffix 'OLS' contrastNames{c}],'contrastvariance',outputExt)];
            fpp.util.writeDataMatrix(conVars(c,:)',hdr,outputPath,maskVol);
        end
        
        % OLS error variance and R^2
        outputPath = [outputDir '/' fpp.bids.changeName(inputName,'desc',...
            [outputSuffix 'OLS'],'errorvariance',outputExt)];
        fpp.util.writeDataMatrix(errorVar',hdr,outputPath,maskVol);
        outputPath = [outputDir '/' fpp.bids.changeName(inputName,'desc',...
            [outputSuffix 'OLS'],'rsquared',outputExt)];
        fpp.util.writeDataMatrix(rSquareds',hdr,outputPath,maskVol);
        
        % Residual data
        if writeResiduals
            outputPath = [outputDir '/' fpp.bids.changeName(inputName,'desc',...
                [outputSuffix 'Residuals'],'bold',outputExtSeries)];
            fpp.util.writeDataMatrix(resids',hdr,outputPath,maskVol);
        end
    end
    
    if mod(iter,10)==0
        fprintf('%s\n',['Iter ' int2str(iter) ' - ' outputNameGeneric]);
    end
end

if permIters>0, save(outputMat,'-append','blockPerm','conVarBasePerm'); end

end