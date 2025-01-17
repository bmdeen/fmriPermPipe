
% [psc,condNames,runNames] = fpp.func.roiExtractAnat(extractResponseDir,searchPath,varargin)
%
% Function to extract region-of-interest responses from a given task and
% participant, using anatomically defined ROIs. To use functional ROIs
% defined with an anatomical constraint, use fpp.func.roiExtract.
%
% Arguments:
% - extractResponseDir (string): model directory for any run of the task
%       used to extract responses (modelarma, modelfilm, or modelperm)
% - searchPath (string): path to the search space (i.e., anatomical ROI)
%
% Variable arguments:
% - overwrite (boolean): whether to overwrite existing ROIs
% - maskPath (string): path to brain mask, to intersect with search space
%
% Outputs:
% - psc: run x condition matrix of PSC (% signal change)
% - condNames: cell array of condition names
% - runNames: cell array of run names

function [psc,condNames,runNames] = roiExtractAnat(extractResponseDir,searchPath,varargin)

psc = [];

% Variable argument defaults
overwrite = 0;
maskPath = '';

% Edit variable arguments.  Note: optInputs checks for proper input.
varArgList = {'maskPath','overwrite'};
for i=1:length(varArgList)
    argVal = fpp.util.optInputs(varargin,varArgList{i});
    if ~isempty(argVal)
        eval([varArgList{i} ' = argVal;']);
    end
end

% Get model directory desc/task info
[analysisDir,extractResponseName,~] = fpp.util.fileParts(extractResponseDir);
% extractResponseDesc = fpp.bids.checkNameValue(extractResponseName,'desc');
% extractResponseTask = fpp.bids.checkNameValue(extractResponseName,'task');
% subjID = fpp.bids.checkNameValue(extractResponseName,'sub');

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

% Extract responses
roiPath = searchPath;
for r=1:nRuns
    [psc(r,:),condNames] = fpp.func.roiExtractWorker(extractResponseDirs{r},roiPath);
end

end