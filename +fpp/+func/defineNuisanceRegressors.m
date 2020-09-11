
% Function to define nuisance time series for fMRI dataset.
%
%%% CURRENTLY UNDERGOING MAJOR OVERHAUL
%%% CONVERTING TO BIDS NAMING CONVENTION
%
%
%
%
%
%
% Incorporates artifact time point removal (these time points are ignored
% in computing noise PCs), but output time series include values at
% artifact time points.
% 
% By default: use white matter and CSF principal components, linear trend, 
% and motion regressors. Assume default paths, with option to customize.
%
%
% TODO:
% - Load regressors via BIDS confounds file by default

function nuisRegrMat = defineNuisanceRegressors(dataPath,varargin)

% Load/check config variables.
configError = fpp.util.checkConfig;
if ~isempty(configError)
    fprintf('%s\n',configError);
    return;
end

[dataDir,~,~] = fileparts(dataPath);

% Assumed paths based on fMRIPermPipe directory structure. Assuming we're using data from the .prep/reg_target directory
roiPathWM = [dataDir '/../../../roi/target_func/wmmask_erode1.nii.gz'];     % Path to white matter ROI
roiPathCSF = [dataDir '/../../../roi/target_func/wmmask_erode1.nii.gz'];    % Path to CSF ROI
maskPath = [dataDir '/mask.nii.gz'];                                        % Path to brain mask
motionParPaths = dir([dataDir '/../mc/prefiltered_func_data*_mcf.par']);
if ~isempty(motionParPaths)
    motionParPath = [dataDir '/../mc/' motionParPaths(1).name];             % Path to motion parameter file
else
    motionParPath = [];
end
badVolPath = [dataDir '/../art/badvols'];                                   % Path to file with list of bad volumes.
disdaqPath = [dataDir '/../art/disdaqs'];                                   % Path to file with list of disdaqs
unsmoothedDataPath = [dataDir '/prefiltered_func_data_bet.nii.gz'];         % Path to unsmoothed data. Must have same # of time points as smoothed data

% Nuisance regressor parameters
useNuisLinear = 1;          % Whether to remove linear trend regressor
useNuisWMPCA = 1;           % Whether to remove white matter PCs
useNuisCSFPCA = 1;          % Whether to remove CSF PCs
pcaOrder = 5;               % Number of PCA regressors to use, across WM and/or CSF
pcaUnsmoothed = 1;          % Use unsmoothed data for PCA
useNuisMotion = 1;          % Whether to remove motion parameters
useNuisGMR = 0;             % Whether to remove the global mean signal (time series averaged across brain mask)
customNuisRegr = [];        % Custom nuisance regressors
erodeROI = 0;               % Whether to erode WM/CSF ROIs before extracting signals (default: assume ROIs are already eroded)
disdaqs = [];               % Number of disdaqs to remove from start of run, overrides disdaqPath if specified
removeBadVols = 1;          % Whether to remove artifact time points before computing nuisance regression betas
                            % Should be used if bad volumes are being censored at analysis step

% Regressor filtering parameters
tempFilt = 0;               % Whether to highpass filter regressors (use if data were filtered, and badVolsRemoved==0)
filtCutoff = [.01 .1];      % Bandpass filter cutoffs (Hz)
filtOrder = 60;             % Hanning-window FIR filter order

% Edit variable arguments.  Note: optInputs checks for proper input.
varArgList = {'roiPathWM','roiPathCSF','maskPath','motionParPath',...
    'badVolPath','useNuisLinear','useNuisWMPCA','useNuisCSFPCA','pcaOrder','useNuisMotion',...
    'useNuisGMR','customNuisRegr','disdaqPath','disdaqs','tempFilt','filtCutoff','filtOrder',...
    'removeBadVols','unsmoothedDataPath','outputMat'};
for i=1:length(varArgList)
    argVal = fpp.util.optInputs(varargin,varArgList{i});
    if ~isempty(argVal)
        eval([varArgList{i} ' = argVal;']);
    end
end

% Check that required inputs exist!
if ~exist(dataPath,'file')
    fprintf('%s\n',['ERROR: input data ' dataPath ' does not exist.']);
    return;
end
if ~exist(maskPath,'file')
    fprintf('%s\n',['ERROR: brain mask ' maskPath ' does not exist.']);
    return;
end
if useNuisWMPCA && ~exist(roiPathWM,'file')
    fprintf('%s\n',['ERROR: white matter mask ' roiPathWM ' does not exist.']);
    return;
end
if useNuisCSFPCA && ~exist(roiPathCSF,'file')
    fprintf('%s\n',['ERROR: CSF mask ' roiPathCSF ' does not exist.']);
    return;
end
if useNuisMotion && (isempty(motionParPath) || ~exist(motionParPath,'file'))
    fprintf('%s\n',['ERROR: motion parameter file ' motionParPath ' does not exist.']);
    return;
end
if pcaUnsmoothed && (useNuisCSFPCA || useNuisWMPCA) && ~exist(unsmoothedDataPath,'file')
    fprintf('%s\n',['ERROR: unsmoothed data ' unsmoothedDataPath ' does not exist.']);
    return;
end
if removeBadVols && ~exist(badVolPath,'file')
    fprintf('%s\n',['ERROR: bad volume file ' badVolPath ' does not exist.']);
    return;
end
if isempty(disdaqs) && ~exist(disdaqPath,'file')
    disdaqs = 0;
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
vols = fpp.util.checkMRIProperty('vols',dataPath);
goodVols = setdiff(1:numVols,union(disdaqVols,badVols));

% Load data/mask
funcData = fpp.util.readMRI(dataPath);
dims = size(funcData.vol);
maskData = fpp.util.readMRI(maskPath);
funcMat = reshape(funcData.vol,[prod(dims(1:3)) dims(4)])';
funcMat = funcMat(:,maskData.vol==1);
newData = funcData;
newData.vol = zeros(size(newData.vol));

% Define nuisance regressors
nuisRegrMat = [];

% Linear trend
if useNuisLinear
    % Linear trend regressor
    linRegr = linspace(0,1,numVols);
    linRegr = linRegr - mean(linRegr);
    nuisRegrMat = [nuisRegrMat linRegr'];
end

% Six motion parameters
if useNuisMotion
    moRegr = load(motionParPath);
    nuisRegrMat = [nuisRegrMat moRegr];
end

% Global mean
if useNuisGMR
    gmrRegr = mean(funcMat,2);
    gmrRegr = gmrRegr - mean(gmrRegr);
    nuisRegrMat = [nuisRegrMat gmrRegr];
end

% Add custom nuisance regressors
if ~isempty(customNuisRegr)
    if size(customNuisRegr,1)~=numVols
        fprintf('%s\n',['NOTE: Custom nuisance regressor matrix must '...
            'have first dimension equal to # of volumes in raw '...
            'dataset.  Not using these regressors.']);
    end
    nuisRegrMat(:,end+1:end+size(customNuisRegr,2)) = customNuisRegr;
end

% Extract WM/CSF PCA regressors.
if useNuisWMPCA || useNuisCSFPCA
    
    % Load unsmoothed functional data
    if pcaUnsmoothed
        funcDataTmp = fpp.util.readMRI(unsmoothedDataPath);
        funcMatTmp = reshape(funcDataTmp.vol,[prod(dims(1:3)) dims(4)])';
    else
        funcMatTmp = funcMat;
    end
    
    noiseMat = [];  % Time point by voxel matrix, corresponding to noise volume
    
    if useNuisWMPCA
        wmData = fpp.util.readMRI(roiPathWM);
        if erodeROI==1, wmData.vol = erodeROIImg(wmData.vol); end
        wmInd = find(wmData.vol==1);
        wmMat = funcMatTmp(:,wmInd);
        noiseMat = [noiseMat wmMat];
    end
    
    if useNuisCSFPCA
        csfData = fpp.util.readMRI(roiPathCSF);
        if erodeROI==1, csfData.vol = erodeROIImg(csfData.vol); end
        csfInd = find(csfData.vol==1);
        csfMat = funcMatTmp(:,csfInd);
        noiseMat = [noiseMat csfMat];
    end
    
    noiseMat = bsxfun(@minus,noiseMat,mean(noiseMat));  % Remove mean at each voxel
    
    % Orthogonalize w.r.t. existing nuisance regressors (doesn't take into account disdaqs)
    % noiseMat = noiseMat - nuisRegrMat*inv(nuisRegrMat'*nuisRegrMat)*nuisRegrMat'*noiseMat;
    
    % Old PCA method, doesn't remove disdaq time points before computing PCs
    %[pcaRegr,~] = svd(noiseMat,'econ');
    
    [~,~,V] = svd(noiseMat(setdiff(1:numVols,disdaqs),:),'econ');   % Find spatial weights of PCs after removing disdaq time points
    pcaRegr = (noiseMat*V)./sqrt(sum((noiseMat*V).^2));             % Project full dataset (all time points) onto spatial components V
                                                                    % noiseMat = pcaRegr*diag(sqrt(sum((noiseMat*V).^2)))*V'
    
    nuisRegrMat(:,end+1:end+pcaOrder) = pcaRegr(:,1:pcaOrder);
    pcaRegrInds = size(nuisRegrMat,2)+1:size(nuisRegrMat,2)+pcaOrder;

end

% Demean nuisance regressors
nuisRegrMat = bsxfun(@minus,nuisRegrMat,mean(nuisRegrMat));

% Filter nuisance regressors
if tempFilt
    % If we're using the data itself to compute PCs, assume that it has already been temporally filtered
    regrIndsToFilt = 1:size(nuisRegrMat,2);
    if ~pcaUnsmoothed
        regrIndsToFilt = setdiff(regrIndsToFilt,pcaRegrInds);
    end
    
    tr = fpp.utils.checkMRIProperty('tr',dataPath);
    nuisRegrMat(:,regrIndsToFilt) = firFilter(nuisRegrMat(:,regrIndsToFilt),filtCutoff,1/tr);
end

end