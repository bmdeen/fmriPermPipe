
% [psc,condNames,runNames] = fpp.func.roiExtract(extractResponseDir,defineROIDir,contrastName,searchPath,varargin)
%
% Function to extract region-of-interest responses from a given task and
% participant, using ROIs defined by maximally responsive coordinates from
% a specific contrast in another task, within a search space. If the two
% tasks are the same, a leave-one-run-out analysis is performed, to ensure
% that separate data are used to define ROIs and extract responses.
%
% Arguments:
% - extractResponseDir (string): modelperm directory for any run of the
%       task used to extract responses.
% - defineROIDir (string): modelperm directory for any run of the task used
%       to define ROIs. The function assumes that extractResponseDir and
%       defineROIDir share a parent directory.
% - contrastName (string): name of the contrast to use for ROI definition
% - searchPath (string): path to the search space
%
% Variable arguments:
% - roiSize (scalar): size of ROI, in % or # of coords in a search space
% - sizeType ('pct' or 'num'): whether size is measured in # or % of coords
% - invertStats (boolean, default = 0): whether to invert statistical map
% - maskPath (string): path to brain mask, to intersect with search space
%
% Outputs:
% - psc: run x condition matrix of PSC (% signal change)
% - condNames: cell array of condition names
% - runNames: cell array of run names

function [psc,condNames,runNames] = roiExtract(extractResponseDir,defineROIDir,contrastName,searchPath,varargin)

psc = [];

% Variable argument defaults
roiSize = 5;
sizeType = 'pct';
invertStats = 0;
maskPath = '';

% Edit variable arguments.  Note: optInputs checks for proper input.
varArgList = {'roiSize','sizeType','invertStats','maskPath'};
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

% Get search space info
[~,searchName,~] = fpp.util.fileParts(searchPath);
searchDesc = fpp.bids.checkNameValue(searchPath,'desc');

% Get modelperm directory desc/task info
[analysisDir,extractResponseName,~] = fpp.util.fileParts(extractResponseDir);
[~,defineROIName,~] = fpp.util.fileParts(defineROIDir);
extractResponseDesc = fpp.bids.checkNameValue(extractResponseName,'desc');
defineROIDesc = fpp.bids.checkNameValue(defineROIName,'desc');
extractResponseTask = fpp.bids.checkNameValue(extractResponseName,'task');
defineROITask = fpp.bids.checkNameValue(defineROIName,'task');

% If the same task is used to define ROIs and extract data, use a
% leave-one-run-out analysis to avoid non-independence errors.
if strcmp(extractResponseTask,defineROITask)
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
dirs = fpp.util.regExpDir(analysisDir,regexprep(defineROIName,'run-[a-zA-Z0-9]+_','run-[a-zA-Z0-9]+_'));
for d=1:length(dirs)
    defineROIDirs{d} = [analysisDir '/' dirs(d).name];
    defineROIRuns{d} = fpp.bids.checkNameValue(dirs(d).name,'run');
end
[~,sortInd] = sort(cellfun(@str2num,defineROIRuns));  % Sort runs by run #
defineROIDirs = defineROIDirs(sortInd);
defineROIRuns = defineROIRuns(sortInd);

% For leave-one-run-out analysis, make sure the same runs were identified
% for defineROI and extractResponse directories
if loro && (~isempty(setdiff(extractResponseRuns,defineROIRuns)) || ...
        ~isempty(setdiff(defineROIRuns,extractResponseRuns)))
    error('Different runs were found for defineROI and extractResponse model directories.');
end

% Determine input z-stat map extension (can be NIFTI or CIFTI)
paths = dir([defineROIDir '/' fpp.bids.changeName(defineROIName,'desc',...
    [defineROIDesc 'OLS' contrastName],'zstat','') '.*nii*']);
if isempty(paths), error('Could not find OLS z-stat map in defineROIDir.'); end
[~,~,inputExt] = fpp.util.fileParts(paths(1).name);

% Find all OLS z-stat maps to use for defining ROI
for r=1:length(defineROIDirs)
    zStatPaths{r} = [defineROIDirs{r} '/' fpp.bids.changeName(defineROIName,{'run','desc'},...
        {defineROIRuns{r},[defineROIDesc 'OLS' contrastName]},'zstat',inputExt)];
    if ~exist(zStatPaths{r},'file')
        error(['Expected input z-statistic path ' zStatPaths{r} ' does not exist.']);
    end
end

% Define ROI directory
roiDir = [analysisDir '/../roi'];
if ~exist(roiDir,'dir'), mkdir(roiDir); end

% Extract responses
if loro
    for r=1:nRuns
        % Define ROI, using average statistical map from all but one run
        loroSuffix = ['LORO' defineROIRuns{r}];
        roiPath = [roiDir '/' fpp.bids.changeName(searchName,'desc',...
            [searchDesc defineROITask defineROIDesc contrastName invertSuffix numSuffix loroSuffix]) inputExt];
        if ~exist(roiPath,'file')
            fpp.util.defineROI(zStatPaths(setdiff(1:nRuns,r)),searchPath,roiPath,'roiSize',roiSize,...
                'sizeType',sizeType,'invertStats',invertStats,'maskPath',maskPath);
        end
        
        [psc(r,:),condNames] = fpp.func.roiExtractWorker(extractResponseDirs{r},roiPath);
    end
else
    % Define ROI, using average statistical map
    roiPath = [roiDir '/' fpp.bids.changeName(searchName,'desc',...
        [searchDesc defineROITask defineROIDesc contrastName invertSuffix numSuffix]) inputExt];
    if ~exist(roiPath,'file')
        fpp.util.defineROI(zStatPaths,searchPath,roiPath,'roiSize',roiSize,...
            'sizeType',sizeType,'invertStats',invertStats,'maskPath',maskPath);
    end
    
    % Extract responses for each run
    for r=1:nRuns
        [psc(r,:),condNames] = fpp.func.roiExtractWorker(extractResponseDirs{r},roiPath);
    end
    
end

end