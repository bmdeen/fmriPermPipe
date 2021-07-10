
% Function to remove nuisance time series from fMRI dataset via linear
% regression.
%
% fpp.func.removeNuisance(inputPath,outputPath,varargin)
%
% Arguments:
% - inputPath (string): path to input data
% - outputPath (string): path to outputData
%
% Variable arguments:
% - maskPath (string): mask to mask image
% - confoundPath (string): path to confound TSV file
% - confoundNames (cell array of strings): names of signals in confound TSV
%       to remove
% - outputDescription: new value for JSON Description field
% - appendDescription (boolean): whether to append outputDescription
%       to existing description
% - outlierPath (string): path to outlier TSV file
% - disdaqPath (string): path to text file with disdaq volume indices
% - disdaqs (numeric vector): disdaq indices; overrides disdaqPath
% - removeBadVols (boolean, default=1): whether to disregard outlier time
%       points when computing nuisance regression beta weights
% 
% NOTE: this script computes beta weights for removing nuisance signals
% without taking into account outlier time points (by default). However, 
% values at artifact time points are influenced by these regressions, and 
% may artificially take on large/small values. If the resulting data are 
% used for any subsequent processing step that doesn't censor artifact time
% points, they should be re-interpolated!

function removeNuisance(inputPath,outputPath,varargin)

% Check system configuration
fpp.util.checkConfig;

% Variable arguments: general
overwrite = 0;           % Whether to overwrite output
maskPath = [];           % Path to brain mask
confoundPath = fpp.bids.changeName(inputPath,'desc',[],'confounds','.tsv'); % Path to confound TSV file
confoundNames = {};
outputDescription = '';
appendDescription = 0;

% Variable arguments: volumes to remove
outlierPath = fpp.bids.changeName(inputPath,{'space','desc'},{[],[]},'outliers','.tsv'); % Path to outlier TSV file
disdaqPath = '';            % Path to file with list of disdaqs
disdaqs = [];               % Number of disdaqs to remove from start of run, overrides disdaqPath if specified
removeBadVols = 1;          % Whether to remove artifact time points before computing nuisance regression betas
                            % Should be used if bad volumes are being censored at analysis step

% Edit variable arguments.  Note: optInputs checks for proper input.
varArgList = {'overwrite','confoundPath','confoundNames','outlierPath','disdaqPath',...
    'disdaqs','removeBadVols','maskPath','outputDescription','appendDescription'};
for i=1:length(varArgList)
    argVal = fpp.util.optInputs(varargin,varArgList{i});
    if ~isempty(argVal)
        eval([varArgList{i} ' = argVal;']);
    end
end

% Define default mask path
if isempty(maskPath) && exist(fpp.bids.changeName(inputPath,'desc','brainNonZero','mask'),'file')
    maskPath = fpp.bids.changeName(inputPath,'desc','brainNonZero','mask');
end

% Check that required inputs exist!
if ~exist(inputPath,'file')
    error(['Input data ' inputPath ' does not exist.']);
end
if ~exist(confoundPath,'file')
    error(['Brain mask ' confoundPath ' does not exist.']);
end
if exist(outputPath,'file') && ~overwrite
    return;
end

% Determine bad volumes, disdaqs, and total number of volumes
if removeBadVols
    badVols = fpp.util.readOutlierTSV(outlierPath);
else
    badVols = [];
end
if isempty(disdaqs)
    if exist(disdaqPath,'file')
        disdaqVols = load(disdaqPath);
    else
        disdaqVols = [];
    end
else
    disdaqVols = (1:disdaqs)';
end
numVols = fpp.util.checkMRIProperty('vols',inputPath);
goodVols = setdiff(1:numVols,union(disdaqVols,badVols));

% Load data/mask
if ~isempty(maskPath)
    maskData = fpp.util.mriRead(maskPath);
    [funcMat,hdr] = fpp.util.readDataMatrix(inputPath,maskData.vol);
else
    [funcMat,hdr] = fpp.util.readDataMatrix(inputPath);
    maskData.vol = ones(size(inputPath,1:3));
end

% Define nuisance regressors
confound = bids.util.tsvread(confoundPath);
if isempty(confoundNames)
    confoundNames = fieldnames(confound);
end
nuisRegrMat = [];
for c=1:length(confoundNames)
    nuisRegrMat(:,end+1) = eval(['confound.' confoundNames{c}]);
end

% Remove nuisance signals from dataset, using a regression strategy that
% omits artifact time points in determining regression weights
X = [ones(length(goodVols),1) nuisRegrMat(goodVols,:)];
funcMat = funcMat';
betas = inv(X'*X)*X'*funcMat(goodVols,:);
funcMat = funcMat - nuisRegrMat*betas(2:end,:);
funcMat = funcMat';

% Write output
fpp.util.writeDataMatrix(funcMat,hdr,outputPath,maskData.vol);

% Reconstruct json metadata
if ~isempty(fpp.bids.getMetadata(inputPath))
    if ~strcmp(inputPath,outputPath)
        fpp.bids.jsonReconstruct(inputPath,outputPath);
    end
    if ~isempty(outputDescription)
        fpp.bids.jsonChangeValue(outputPath,'Description',outputDescription,appendDescription);
    end
    if ~strcmp(inputPath,outputPath)
        fpp.bids.jsonChangeValue(outputPath,'Sources',fpp.bids.removeBidsDir(inputPath));
    end
end

end