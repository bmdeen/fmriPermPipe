
% fpp.func.modelFeat(inputPath,eventsPath,contrastMatrixPath,varargin)
%
% Step 1 of a two-step process (modelFeat, model2ndFeat) to perform a
% General Linear Model based analysis of fMRI data using FSL's FEAT. Should
% be run after fpp.func.preproc.
% 
% Example usage: modelFeat('/pathToData/sub-01_task-faceloc_run-01_space-individual_desc-preproc_bold.nii.gz',...
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
%   - confoundPath (string): path to confound.tsv file with nuisance
%       regressors
%   - confoundNames (cell array of strings): fields in confound file to use
%   - confoundFilt (boolean vector): which confound regressors to
%       temporally filter
%   - outlierPath (string): path to outliers.tsv specifying time points to
%       exclude
%   - useTedana (boolean): whether to include TEDANA removed components as
%       confound regressors (unless data has "NoTedana" in desc)
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
% - zstat: z-statistic images
% - desc-Regressors_image.png: image showing regressor time series
% - desc-RegressorCorrelations_image.png: image showing correlations
%   between regressors
% - desc-RegressorVarianceRemoved_image.png: image showing the proportion
%   of variance in task regressors that is explained by nuisance regressors

function modelFeat(inputPath,eventsPath,contrastMatrixPath,varargin)

% Check system configuration
fpp.util.checkConfig;

% Basic parameters
overwrite = 0;              % Whether to overwrite output
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
plotResults = 0;            % Whether to display result plots
writeResiduals = 0;         % Whether to write 4-D residual image



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 1: Update and validate parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Edit variable arguments.  Note: optInputs checks for proper input.
varArgList = {'overwrite','outputSuffix','permuteRest','tempFilt','filtCutoff',...
    'filtOrder','hrfType','upsampledTR','writeResiduals','permIters','plotResults',...
    'randSeed','condNames','confoundPath','confoundNames','confoundFilt','outlierPath',...
    'contrastNames','filtType','analysisDir','useTedana'};
for i=1:length(varArgList)
    argVal = fpp.util.optInputs(varargin,varArgList{i});
    if ~isempty(argVal)
        eval([varArgList{i} ' = argVal;']);
    end
end

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
outputNameGeneric = fpp.bids.changeName(inputName,'desc',outputSuffix,'bold','');
numVols = fpp.util.checkMRIProperty('vols',inputPath);
tr = fpp.util.checkMRIProperty('tr',inputPath);
totalVoxels = prod(fpp.util.checkMRIProperty('dims',inputPath));
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
if exist(outlierPath,'file')
    outlier = bids.util.tsvread(outlierPath);
    outlierNames = fieldnames(outlier);
else
    outlier = struct();
    outlierNames = {};
end
nOutliers = length(outlierNames);

% Define nuisance regressors
nuisRegrMat = [];
if exist(confoundPath,'file')
    confound = bids.util.tsvread(confoundPath);
    for c=1:length(confoundNames)
        nuisRegrMat(:,end+1) = eval(['confound.' confoundNames{c}]);
    end
    
    if ~isempty(nuisRegrMat)
        nuisRegrMat = bsxfun(@minus,nuisRegrMat,mean(nuisRegrMat));
        
        % Filter nuisance regressors
        if tempFilt
            nuisRegrMat(:,confoundFilt==1) = fpp.util.firFilter(...
                nuisRegrMatnuisRegrMat(:,confoundFilt==1),1/tr,filtCutoff,filtType,filtOrder);
        end
        
        % Redefine confound TSV after temporal filtering.
        confoundNew = struct();
        for c=1:length(confoundNames)
            eval(['confoundNew.' confoundNames{c} ' = nuisRegrMat(:,' int2str(c) ');']);
        end
        confound = confoundNew;
    end
else
    confound = struct();
end
regrNames = [condNames,confoundNames];

% Add outlier to confound struct
for i=1:length(outlierNames)
    eval(['confound.' outlierNames{i} ' = outlier.' outlierNames{i} ');']);
    nuisRegrMat(:,end+1) = eval(['outlier.' outlierNames{i}]);
end

% Load and add TEDANA regressors to confound struct
% NOTE: Not designed to handle tedana dirs with desc field
inputDesc = fpp.bids.checkNameValue(inputName,'desc');
if useTedana && sum(regexpi(inputDesc,'NoTedana'))==0
    tedanaTSVPath = [inputDir '/' fpp.bids.changeName(inputName,'desc','','tedana','')...
        '/' fpp.bids.changeName(inputName,'desc','tedanaICARejected','mixing','.tsv')];
    tedana = fpp.bids.tsvread(tedanaTSVPath);
    tedanaNames = fieldnames(tedana);
    for i=1:length(tedanaNames)
        eval(['confound.' outlierNames{i} ' = tedana.' tedanaNames{i} ');']);
        nuisRegrMat(:,end+1) = eval(['tedana.' tedanaNames{i}]);
    end
end

% Define and create output directory
outputName = fpp.bids.changeName(inputName,'desc',outputSuffix,'modelfeat','');
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
outputMat = [outputDirOrig '/' fpp.bids.changeName(inputName,'desc',outputSuffix,'RegressionData','.mat')];

% Define FEAT output directory
outputFeatDir = [analysisDir '/' fpp.bids.changeName(inputName,'desc',outputSuffix,'model','.feat')];

% Define and create condition TSV file
condPath = [outputDir '/' fpp.bids.changeName(inputName,'desc',outputSuffix,'conditions','.tsv')];
for c=1:nConds
    condTSV.cond_names{c} = condNames{c};
end
bids.util.tsvwrite(condPath,condTSV);

% Define MNI space path
stdPath = [getenv('FSLDIR') '/MNI152_T1_2mm_brain.nii.gz'];

% Define FPP data directory
[fppFuncDir,~,~]		= fileparts(mfilename('fullpath'));			% path to the directory containing this script
tmp = dir([fppFuncDir '/../../data']);
dataDir = tmp(1).folder;

% Define template and output design.fsf files
designTemplatePath = [dataDir '/desc-fpptemplate_design.fsf'];
designPath = [outputDir '/' fpp.bids.changeName(outputName,[],[],'design','.fsf')];

% Load FSF template
fid = fopen(template);
fsfTemplate = [];
while 1
    line = fgetl(fid);
    if ~isnumeric(line)
        fsfTemplate = [fsfTemplate '\n' line];
    else, break; end
end
fclose(fid);



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

% Write task regressors to .txt files
for c=1:nConds
    regrPaths{c} = [outputDir '/' fpp.bids.changeName(outputName,'desc',...
        [outputSuffix condNames{c}],'regressor','.txt')];
    fid = fopen(regrPaths{c});
    fprintf(fid,'%f\n',taskRegrMat(:,c));
    fclose(fid);
end

% Write full confound TSV file
confoundOutputPath = [outputDir '/' fpp.bids.changeName(outputName,'desc',...
    [outputSuffix 'feat'],'confounds','.tsv')];
if ~isempty(fieldnames(confound))
    fpp.bids.tsvwrite(confoundOutputPath,confound);
    confoundFile = ['set confoundev_files(1) "' confoundOutputPath '"'];
    useConfound = 1;
else
    confoundFile = '';
    useConfound = 0;
end

% Plot regressor information.
fpp.func.analysis.plotRegressors(taskRegrMat,nuisRegrMat,regrNames,outputDir,...
    outputNameGeneric,~plotResults);

% Check design matrix X
X = [taskRegrMat nuisRegrMat];
if rank(X)<size(X,2)
    error(['Design matrix is rank-deficient - ' outputNameGeneric]);
end

% Save outputs
save(outputMat,'regressors3Col','X','condNames','tr','inputPath',...
    'confoundPath','outlierPath','hrfType','tempFilt','filtCutoff',...
    'filtOrder','contrastNames','contrastMat');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 3: Write design.fsf file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define EV info for .fsf file
evInfo = '';
for c = 1:nConds
    evInfo = [evInfo sprintf('%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n\n',['# EV ' int2str(c)],...
        ['set fmri(evtitle' int2str(c) ') "' condNames{c} '"'], ...
        ['set fmri(shape' int2str(c) ') 2'],...
        ['set fmri(convolve' int2str(c) ') 0'],...
        ['set fmri(convolve_phase' int2str(c) ') 0'],...
        ['set fmri(tempfilt_yn' int2str(c) ') 0'],...
        ['set fmri(deriv_yn' int2str(c) ') 0'],...
        ['set fmri(custom' int2str(c) ') "' regrPaths{c} '"'])];
    % Set orthogonalization parameters to 0.
    for c2 = 0:nConds
        evInfo = [evInfo sprintf('%s\n',['set fmri(ortho' int2str(c) '.' int2str(c2) ') 0'])];
    end
    evInfo = [evInfo sprintf('%s\n\n','')];
end

% Define con info for .fsf file
conInfo = '';
for con = 1:nContrasts

    conInfo = [conInfo sprintf('%s\n%s\n\n',...
        ['set fmri(conpic_real.' int2str(con) ') 1'],...
        ['set fmri(conname_real.' int2str(con) ') "' contrastNames{con} '"'],...
        ['set fmri(conpic_orig.' int2str(con) ') 1'],...
        ['set fmri(conname_orig.' int2str(con) ') "' contrastNames{con} '"'])];

    for c = 1:nConds
        conInfo = [conInfo sprintf('%s\n\n',...
            ['set fmri(con_real' int2str(con) '.' int2str(c) ') ' int2str(contrastMat(con,c))])];
    end
    for c = 1:nConds
        conInfo = [conInfo sprintf('%s\n\n',...
            ['set fmri(con_orig' int2str(con) '.' int2str(c) ') ' int2str(contrastMat(con,c))])];
    end
    conInfo = [evinfo sprintf('%s\n\n','')];
    
end

conInfo = [conInfo sprintf('%s\n%s\n\n','set fmri(conmask_zerothresh_yn) 0','set fmri(conmask1_1) 0')];

% Define FEAT .fsf variables to change.
names = {'OutputDir',            'TR', 'Vols', 'nConds', 'nContrasts', 'TotalVoxels', 'UseConfound', 'ConfoundFile', 'EvInfo', 'ConInfo', 'StdPath'};
vals =   {outputFeatDir(1:end-5), tr,   vols,   nConds,   nContrasts,   totalVoxels,   useConfoud,    confoundFile,   evInfo,   conInfo,   stdPath};

% Write new values to fsf file
fsf = fsfTemplate;
for j = 1:length(vals)
    curval = vals{j};
    if isnumeric(curval), curval = num2str(curval); end
    fsf = strrep(fsf,['{' names{j} '}'],curval);
end
fid = fopen(designPath,'w');
fprintf(fid,fsf);
fclose(fid);

% Run FEAT analysis via command line
% fpp.util.system(['feat ' designPath]);


% Move files to new location






end