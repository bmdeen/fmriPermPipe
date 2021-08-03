
% [psc,condNames,runNames] = fpp.func.roiExtract(extractResponseDir,defineROIDir,contrastName,searchPath,varargin)
%
% Function to extract region-of-interest responses from a given task and
% participant, using ROIs defined by maximally responsive coordinates from
% a specific contrast in the same or another task, within a search space. 
% If the same tasks are used to define ROIs and extract responses, a
% leave-one-run-out analysis is performed.
%
% Arguments:
% - extractResponseDir (string): model directory for any run of the task
%       used to extract responses (modelarma, modelfilm, or modelperm)
% - defineROIDir (string or cell array of strings): model directory for
%       any one run of the task used to define ROIs. The function assumes
%       that extractResponseDir and defineROIDir share a parent directory.
%       Multiple tasks can be used in combination for defining ROIs, with
%       statistics averaged across tasks. In this case, statCoefs specifies
%       how they should be averaged, and the tasks must have the same
%       number of runs. 
% - contrastName (string or cell array of strings): name of the contrast(s)
%       to use for ROI definition.
% - searchPath (string): path to the search space
%
% Variable arguments:
% - overwrite (boolean): whether to overwrite existing ROIs
% - roiSize (scalar): size of ROI, in % or # of coords in a search space
% - sizeType ('pct' or 'num'): whether size is measured in # or % of coords
% - invertStats (boolean, default = 0): whether to invert statistical map
% - maskPath (string): path to brain mask, to intersect with search space
% - statCoefs (numeric vector): coefficients for statistical map averaging
%       across tasks
% - roiDesc (string): ROI description. If specified, this replaces
%       [defineROITask defineROIDesc contrastName] in the output ROI desc
%
% Outputs:
% - psc: run x condition matrix of PSC (% signal change)
% - condNames: cell array of condition names
% - runNames: cell array of run names

function [psc,condNames,runNames] = roiExtract(extractResponseDir,defineROIDir,contrastName,searchPath,varargin)

psc = [];

% Convert defineROIDir and contrastName to cell arrays, if not already
if ~iscell(defineROIDir), defineROIDir = {defineROIDir}; end
nTasks = length(defineROIDir);
if ~iscell(contrastName), contrastName = {contrastName}; end
if length(contrastName)==1 && nTasks>1
    contrastName = repmat(contrastName,[1 nTasks]);
end

% Variable argument defaults
overwrite = 0;
roiSize = 5;
sizeType = 'pct';
invertStats = 0;
maskPath = '';
statCoefs = [];
roiDesc = '';

% Edit variable arguments.  Note: optInputs checks for proper input.
varArgList = {'roiSize','sizeType','invertStats','maskPath','overwrite','statCoefs','roiDesc'};
for i=1:length(varArgList)
    argVal = fpp.util.optInputs(varargin,varArgList{i});
    if ~isempty(argVal)
        eval([varArgList{i} ' = argVal;']);
    end
end

% Process variable arguments
if ~(strcmpi(sizeType,'pct') || strcmpi(sizeType,'num'))
    error('sizeType must be set to pct or num');
end
if strcmpi(sizeType,'pct') && roiSize>100
    roiSize = 100;
end
if strcmpi(sizeType,'num')
    roiSize = round(roiSize);
    numSuffix = ['Top' num2str(roiSize)];
else
    numSuffix = ['Top' num2str(roiSize) 'Pct'];
end
if invertStats
    invertSuffix = 'Inverted';
else
    invertSuffix = '';
end
if isempty(statCoefs)
    statCoefs = ones(length(defineROIDir),1);
end

% Get search space info
[~,searchName,~] = fpp.util.fileParts(searchPath);
searchDesc = fpp.bids.checkNameValue(searchPath,'desc');

% Get model directory desc/task info
[analysisDir,extractResponseName,~] = fpp.util.fileParts(extractResponseDir);
extractResponseDesc = fpp.bids.checkNameValue(extractResponseName,'desc');
extractResponseTask = fpp.bids.checkNameValue(extractResponseName,'task');
for t=1:nTasks
    [~,defineROIName{t},~] = fpp.util.fileParts(defineROIDir{t});
    defineROIDesc{t} = fpp.bids.checkNameValue(defineROIName{t},'desc');
    defineROITask{t} = fpp.bids.checkNameValue(defineROIName{t},'task');
end

% If the same task is used to define ROIs and extract data, use a
% leave-one-run-out analysis to avoid non-independence errors.
if ismember(extractResponseTask,defineROITask)
    loro = 1;
else
    loro = 0;
end

% Find all extractResponse directories, across runs
dirs = fpp.util.regExpDir(analysisDir,regexprep(extractResponseName,'run-[a-zA-Z0-9]+_','run-[a-zA-Z0-9]+_'));
if isempty(dirs)
    error(['Could not find analysis directories with name ' extractResponseName]);
end
for d=1:length(dirs)
    extractResponseDirs{d} = [analysisDir '/' dirs(d).name];
    extractResponseRuns{d} = fpp.bids.checkNameValue(dirs(d).name,'run');
end
[~,sortInd] = sort(cellfun(@str2num,extractResponseRuns));  % Sort runs by run #
extractResponseDirs = extractResponseDirs(sortInd);
extractResponseRuns = extractResponseRuns(sortInd);
nRuns = length(extractResponseRuns);
runNames = extractResponseRuns;

% Find all defineROI directories, across runs, sort by run number
nRunsEachTask = [];
for t=1:nTasks
    dirs = fpp.util.regExpDir(analysisDir,regexprep(defineROIName{t},'run-[a-zA-Z0-9]+_','run-[a-zA-Z0-9]+_'));
    for d=1:length(dirs)
        defineROIDirs{t}{d} = [analysisDir '/' dirs(d).name];
        defineROIRuns{t}{d} = fpp.bids.checkNameValue(dirs(d).name,'run');
    end
    [~,sortInd] = sort(cellfun(@str2num,defineROIRuns{t}));  % Sort runs by run #
    defineROIDirs{t} = defineROIDirs{t}(sortInd);
    defineROIRuns{t} = defineROIRuns{t}(sortInd);
    nRunsEachTask(end+1) = length(defineROIRuns{t});
end
if length(unique(nRunsEachTask))~=1
    error('Different # of runs were found for defineROI model directories.');
end

% For leave-one-run-out analysis, make sure the same runs were identified
% for defineROI and extractResponse directories
if loro && (~isempty(setdiff(extractResponseRuns,defineROIRuns{1})) || ...
        ~isempty(setdiff(defineROIRuns{1},extractResponseRuns)))
    error('Different # of runs were found for defineROI and extractResponse model directories.');
end

% Determine input z-stat map extension (can be NIFTI or CIFTI)
paths = dir([defineROIDir{1} '/' fpp.bids.changeName(defineROIName{1},'desc',...
    [defineROIDesc{1} contrastName{1}],'zstat','') '.*nii*']);
if isempty(paths), error('Could not find z-stat map in defineROIDir.'); end
[~,~,inputExt] = fpp.util.fileParts(paths(1).name);

% Find all z-stat maps to use for defining ROI
for t=1:nTasks
    for r=1:length(defineROIDirs{t})
        zStatPaths{t}{r} = [defineROIDirs{t}{r} '/' fpp.bids.changeName(defineROIName{t},{'run','desc'},...
            {defineROIRuns{t}{r},[defineROIDesc{t} contrastName{t}]},'zstat',inputExt)];
        if ~exist(zStatPaths{t}{r},'file')
            error(['Expected input z-statistic path ' zStatPaths{t}{r} ' does not exist.']);
        end
    end
end

% Define ROI directory
roiDir = [analysisDir '/../roi'];
if ~exist(roiDir,'dir'), mkdir(roiDir); end

% Define task/contrast component of ROI description
if isempty(roiDesc)
    for t=1:nTasks
        roiDesc = [roiDesc defineROITask{t} defineROIDesc{t} contrastName{t}];
    end
end

% Extract responses
if loro
    for r=1:nRuns
        % Define ROI, using average statistical map from all but one run
        loroSuffix = ['LORO' defineROIRuns{1}{r}];
        roiPath = [roiDir '/' fpp.bids.changeName(searchName,'desc',...
            [searchDesc roiDesc invertSuffix numSuffix loroSuffix]) inputExt];
        if ~exist(roiPath,'file') || overwrite
            zStatPathsInput = {};
            statCoefsInput = [];
            for t=1:nTasks
                zStatPathsInput = [zStatPathsInput zStatPaths{t}(setdiff(1:nRuns,r))];
                statCoefsInput = [statCoefsInput repmat(statCoefs(t),[1 nRuns-1])/(nRuns-1)];
            end
            
            fpp.util.defineROI(zStatPathsInput,searchPath,roiPath,'roiSize',roiSize,...
                'sizeType',sizeType,'invertStats',invertStats,'maskPath',maskPath,...
                'statCoefs',statCoefsInput);
        end
        
        [psc(r,:),condNames] = fpp.func.roiExtractWorker(extractResponseDirs{r},roiPath);
    end
else
    % Define ROI, using average statistical map
    roiPath = [roiDir '/' fpp.bids.changeName(searchName,'desc',...
        [searchDesc roiDesc invertSuffix numSuffix]) inputExt];
    if ~exist(roiPath,'file') || overwrite
        
        % Define inputs for each run
        statCoefsInput = [];
        zStatPathsInput = {};
        for t=1:nTasks
            for r=1:nRunsEachTask(t)
                zStatPathsInput{end+1} = zStatPaths{t}{r};
                statCoefsInput(end+1) = statCoefs(t)/length(defineROIRuns{t});
            end
        end
        fpp.util.defineROI(zStatPathsInput,searchPath,roiPath,'roiSize',roiSize,...
            'sizeType',sizeType,'invertStats',invertStats,'maskPath',maskPath,...
            'statCoefs',statCoefsInput);
    end
    
    % Extract responses for each run
    for r=1:nRuns
        [psc(r,:),condNames] = fpp.func.roiExtractWorker(extractResponseDirs{r},roiPath);
    end
    
end

end