
% Function to remove nuisance time series from fMRI dataset via linear
% regression.
% 
% By default: use white matter and CSF PCs, linear trend, and motion 
% regressors. Assume default paths, with option to customize.
% 
% Note: optional arguments for defineNuisanceRegressors can also be passed
% through this function.
% 
% Note: the current script computes weights for removing nuisance signals
% without taking into account artifact time points (by default). However, 
% values at artifact time points are influenced by these regressions, and 
% may artificially take on large/small values. If the resulting data are 
% used for any subsequent processing step that doesn't censor artifact time
% points, they should be re-interpolated!
%
%
% TODO: Modify method for reading artifact time point data.

function removeNuisance(inputPath,varargin)

% Load/check config variables.
configError = fpp.util.checkConfig;
if ~isempty(configError)
    fprintf('%s\n',configError);
    return;
end

[inputDir,inputName,inputExt] = fpp.util.fileParts(inputPath);

if ~ismember(lower(inputExt),{'.nii','.nii.gz'})
    fprintf('%s\n','ERROR: inputPath must be a NIFTI file.');
    return;
end

% If desc field exists, add NuisRemoved to the end.
[startInd endInd] = regexp(inputName,'_desc-.*_');
if ~isempty(startInd)
    newDesc = [inputName(startInd:end-1) 'NuisRemoved_'];
    outputName = regexprep(inputName,'_desc-.*_',newDesc);
else
    outputName = fpp.bids.changeName(inputName,'desc','NuisRemoved','bold');
end
outputPath = [inputDir '/' outputName '.nii.gz'];

% Variable arguments: general
overwrite = 0;          % Whether to overwrite output
outputMat = '';            
maskPath = [inputDir '/mask.nii.gz'];                                        % Path to brain mask                                         % Path to output MAT file (optional)

% Variable arguments: volumes to remove
badVolPath = [inputDir '/../art/badvols'];                                   % Path to file with list of bad volumes.
disdaqPath = [inputDir '/../art/disdaqs'];                                   % Path to file with list of disdaqs
disdaqs = [];               % Number of disdaqs to remove from start of run, overrides disdaqPath if specified
removeBadVols = 1;          % Whether to remove artifact time points before computing nuisance regression betas
                            % Should be used if bad volumes are being censored at analysis step

% Edit variable arguments.  Note: optInputs checks for proper input.
varArgList = {'overwrite','outputPath','badVolPath','disdaqPath','disdaqs','removeBadVols','maskPath'};
for i=1:length(varArgList)
    argVal = fpp.util.optInputs(varargin,varArgList{i});
    if ~isempty(argVal)
        eval([varArgList{i} ' = argVal;']);
    end
end

% Check that required inputs exist!
if ~exist(inputPath,'file')
    fprintf('%s\n',['ERROR: input data ' inputPath ' does not exist.']);
    return;
end
if ~exist(maskPath,'file')
    fprintf('%s\n',['ERROR: brain mask ' maskPath ' does not exist.']);
    return;
end
if exist(outputPath,'file') && ~overwrite
    fprintf('%s\n',['Output path ' outputPath ' already exists.']);
    return;
end

% Determine bad volumes, disdaqs, and total number of volumes
if removeBadVols
    badVols = load(badVolPath);
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
funcData = fpp.util.mriRead(inputPath);
dims = size(funcData.vol);
maskData = fpp.util.mriRead(maskPath);
funcMat = reshape(funcData.vol,[prod(dims(1:3)) dims(4)])';
funcMat = funcMat(:,maskData.vol==1);
newData = funcData;
newData.vol = zeros(size(newData.vol));

% Define nuisance regressors
nuisRegrMat = fpp.func.defineNuisanceRegressors(inputPath,varargin{:});

% Remove nuisance signals from dataset, using a regression strategy that
% omits artifact time points
X = [ones(length(goodVols),1) nuisRegrMat(goodVols,:)];
betas = inv(X'*X)*X'*funcMat(goodVols,:);
funcMat = funcMat - nuisRegrMat*betas(2:end,:);

% Write output data
funcMatAllVoxels = zeros(dims(4),prod(dims(1:3)));
funcMatAllVoxels(:,maskData.vol==1) = funcMat;
newData.vol = reshape(funcMatAllVoxels',dims);
fpp.util.mriWrite(newData,outputPath);

% Save output .mat file with regressor info, if specified
if ~isempty(outputMat)
    save(outputMat,'nuisRegrMat','badVols','disdaqVols');
end

end