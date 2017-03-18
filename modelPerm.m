
% modelPerm(study,subjects,expt,conMat,varargin)
% 
% Step 1 of a two-step process (modelPerm, model2ndPerm) to perform a
% General Linear Model based analysis of fMRI data, computing statistics
% using a permutation test.  Step 1 performs voxelwise time series 
% regressions within runs, for a large number of permuted block orderings.
% Linear and PCA-based (aCompCorr) nuisance regressors are included by
% default.  Note that this step on its own does not output
% permutation-based statistics, which is done by model2ndPerm. Should be 
% run after registerMRI.  Loops through subjects and runs of a given 
% experiment by default.
% 
% Example usage: modelPerm('studyName','SUB*','FaceLoc',[1 0 -1; -1 0 1])
% 
% Critical output files (in *.perm directory):
% - cope*.nii.gz: contrast ("contrast of parameter estimate") images
% - pe*.nii.gz: parameter estimate (beta) images
% - zstat*_ols.nii.gz: z-statistics computed based on ordinary least
%   squares.  These are not valid statistics, due to the presence of
%   temporal autocorrelation, but can give a rough sense of where effects
%   are located.
% - regr_plots.png: image showing regressor time series
% - regr_correlations.png: image showing correlations between regressors
% - regr_variance_removed.png: image showing the proportion of variance in
%   task regressors that is explained by nuisance regressors
% - perms: directory containing contrast images for each permutation
%
% Required arguments:
% - study (string): name of the study to analyze, used to define the study 
%       directory as [getenv('FMRI_BASE_DIR') '/' study].
% - subjects (string or cell array of strings): list of subjects to
%       analyze. Can use asterisk-based regular expression.  Examples:
%       {'SUB01','SUB02'} or 'SUB*'.
% - expt (string): name of experiment to analyze (e.g. 'FaceLocalizer')
% - conMat: matrix of contrast coefficients (# contrasts by # regressors).
% 
% Variable arguments:
% - overwrite (boolean; default=0): whether to overwrite files that have
%       already been written by this function.
% - runList (vector of integers in (0,Inf)): list of runs to process.
% - inputSuffix (string): suffix for input .prep directory.
% - outputSuffix (string): suffix for .perm directory. Can be used to run
%       the script multiple times with different options.
% - permIters (integer in (0,Inf), default=5000): iterations of permutation
%       test. It is recommended to use at least 5000.
% - useNuisLinear (boolean; default=1): whether to include linear trend 
%       nuisance regressor.
% - useNuisWMPCA (boolean; default=1): whether to include PCA-based
%       nuisance regressors from white matter.
% - useNuisCSFPCA (boolean; default=1): whether to include PCA-based
%       nuisance regressors from cerebrospinal fluid.
% - pcaOrder (integer in [0,Inf); default=5): number of PCA regressors to 
%       use, across WM and/or CSF.
% - useNuisMotion (boolean; default=0): whether to include motion 
%       parameters as nuisance regressors.
% - customNuisRegr (numeric matrix): matrix of custom nuisance regressor
%       time series.  Must be N time points (in raw data) by R regressors.
% - permuteRest (boolean; default=1): whether to permute baseline blocks as
%       well as task blocks.  See README (Notes on Theory and Practice
%       section) for further discussion.
% - hrfType (1 or 2; default=1): which hemodynamic response function to use
%       for regressor convolution. 1 = double-gamma (FSL/SPM default 
%       parameters); 2 = gamma (FS-FAST default parameters).
% - upsampledTR (scalar in (0,Inf); default=.01): regressor sampling rate 
%       (s) used for convolution. Code requires mod(tr/2,upsampledTR)==0.
% - forceTR (boolean; default=0): whether to force TR to default value, 
%       instead of reading from header.  Should be used only if TR
%       information in image header is corrupted.
% - defaultTR (scalar in (0,Inf); default=2): default TR (s), used if no TR
%       info is found in functional data header, or if forcetr==1.
% - tempFilt (boolean; default=0): whether to highpass filter regressors.
%       Should be used if data were filtered by preprocMRI.
% - filtCutoff (scalar in (0,Inf); default=.01: highpass filter cutoff (Hz)
% - filtOrder (integer in (0,Inf); default=60): Hanning-window FIR filter
%       order.
% - maskName (string): name of mask image in roi/$TARGET_FUNC directory,
%       without file extension.  By default, a whole-brain mask will be
%       used, but this option can be used to specify a more restrictive
%       mask, to boost power for cluster-level inference.
% - targetName (string; default='target_func'): name of image in func 
%       directory that was used for registration by preprocMRI, without 
%       file extension. If a custom image was specified for preprocMRI,
%       the same image should be used here.
% - plotResults (boolean; default=1): whether to display result plots
% - writeResiduals (boolean; default=0): whether to write 4-D residual image

function modelPerm(study,subjects,expt,conMat,varargin)

addpath([strrep(mfilename('fullpath'),mfilename,'') '/utils']);

% Load/check config variables.
[configError, studyDir, fslPrefix] = checkConfig(study);
if ~isempty(configError)
    fprintf('%s\n',configError);
    return;
end;

% Basic parameters
overwrite = 0;              % Whether to overwrite output
runList = [];               % Runs to analyze.
inputSuffix = '';           % Suffix for preproc dir
outputSuffix = '';          % New suffix for modeling dir
permIters = 5000;           % Iterations of permutation test

% Nuisance regressor parameters
useNuisLinear = 1;          % Whether to include linear trend regressor
useNuisWMPCA = 1;           % Whether to include WM for PCA-based denoising
useNuisCSFPCA = 1;          % Whether to include CSF for PCA-based denoising
pcaOrder = 5;               % Number of PCA regressors to use, across WM and/or CSF
useNuisMotion = 0;          % Whether to include motion parameters as nuisance regressors
customNuisRegr = [];        % Custom nuisance regressors

% Modeling parameters
permuteRest = 1;            % Whether to permute baseline blocks as well as task blocks
hrfType = 1;                % 1 = double-gamma (FSL/SPM default parameters); 2 = gamma (FS-FAST default parameters)
upsampledTR = .01;          % Regressor sampling rate (s) used for convolution. Code requires mod(tr/2,upsampledTR)==0
forceTR = 0;                % Whether to force TR to default value, instead of reading from header
defaultTR = 2;              % Default TR (s), assumed if no TR info in header, or if forcetr==1
subtractHalfTR = 1;         % Whether to subtract .5*TR to regressor onsets, to account for slice timing

% Regressor filtering parameters
tempFilt = 0;               % Whether to highpass filter regressors (use if data were filtered)
filtCutoff = 1/100;         % Highpass filter Cutoff (Hz)
filtOrder = 60;             % Hanning-window FIR filter order

% Mask/registration parameters
maskName = '';              % Name of mask image in roi/$TARGET_FUNC directory
targetName = 'target_func'; % Name of target_func image used for registration

% Data generation parameters
plotResults = 1;            % Whether to display result plots
writeResiduals = 0;         % Whether to write 4-D residual image


% Generate random seed based on clock time
if exist('rng','file')
    rng(sum(100*clock),'twister');
else
    rand('twister',sum(100*clock));
end;

% Convert list of subjects from "dir" inputs to actual names
subjectsNew = {};
if ~iscell(subjects), subjects = {subjects}; end;
for s=1:length(subjects)
    if ~ischar(subjects{s})
        fprintf('%s\n','ERROR: subject argument must be a string or cell array of strings.');
        return;
    end;
    if strfind(subjects{s},'*')
        subjectsTmp = dir([studyDir '/' subjects{s}]);
        for j=1:length(subjectsTmp)
            subjectsNew{end+1} = subjectsTmp(j).name;
        end;
    else
        subjectsNew{end+1} = subjects{s};
    end;
end;
subjects = subjectsNew;

% Edit variable arguments.  Note: optInputs checks for proper input.
varArgList = {'overwrite','runList','inputSuffix','outputSuffix','maskName',...
    'useNuisLinear','useNuisMotion','useNuisWMPCA','useNuisCSFPCA',...
    'pcaOrder','permuteRest','tempFilt','filtCutoff','filtOrder','hrfType',...
    'upsampledTR','defaultTR','forceTR','customNuisRegr','writeResiduals',...
    'permIters','plotResults'};
for i=1:length(varArgList)
    argVal = optInputs(varargin,varArgList{i});
    if ~isempty(argVal)
        eval([varArgList{i} ' = argVal;']);
    end;
end;

if ~(isnumeric(conMat) && ~isscalar(conMat) && numel(size(conMat))==2)
    fprintf('%s\n','ERROR: conMat must be a numeric matrix.');
    return;
end;
numCons = size(conMat,1);

if tempFilt && (useNuisLinear || useNuisWMPCA || useNuisCSFPCA)
    fprintf('%s\n',['WARNING: temporal filtering is not recommended when linear '...
        'or PCA-based nuisance regressors are in use.']);
end;

for s = 1:length(subjects)
    
    subject = subjects{s};
    subjDir = [studyDir '/' subject];
    funcDir = [subjDir '/func'];
    analDir = [subjDir '/analysis'];
    regrDir = [subjDir '/regressors'];
    roiTargDir = [subjDir '/roi/' targetName];
    wmMask = [roiTargDir '/wmmask_erode1.nii.gz'];
    csfMask = [roiTargDir '/csfmask_erode1.nii.gz'];
    
    scanlogFiles = regExpDir(subjDir,'.*scanlo.*([^~]$)');
    
    if exist(subjDir,'dir')==0
        fprintf('%s\n\n',['Subject directory does not exist for ' subject '. Skipping this subject.']);
        continue;
    elseif isempty(scanlogFiles)
        fprintf('%s\n\n',['Missing scanlog file for subject ' subject '. Skipping this subject.']);
        continue;
    elseif exist(funcDir,'dir')==0
        fprintf('%s\n\n',['Functional directory does not exist for ' subject '. Skipping this subject.']);
        continue;
    elseif exist(regrDir,'dir')==0
        fprintf('%s\n\n',['Regressor directory does not exist for ' subject '. Skipping this subject.']);
        continue;
    elseif ~exist(wmMask,'file') && useNuisWMPCA==1
        fprintf('%s\n\n',['White matter ROI does not exist for ' subject '. Skipping this subject.']);
        continue;
    elseif ~exist(csfMask,'file') && useNuisCSFPCA==1
        fprintf('%s\n\n',['CSF ROI does not exist for ' subject '. Skipping this subject.']);
        continue;
    end;
    
    if ~isempty(maskName)
        funcMaskInit = [roiTargDir '/' maskName '.nii.gz'];
        if ~exist(funcMaskInit,'file')
            fprintf('%s\n',['Requested mask file does not exist for ' subject ...
                '. Using default brain mask instead.']);
        end;
    end;
    
    for scan = 1:length(scanlogFiles)
        
        fileID = fopen([subjDir '/' scanlogFiles(scan).name]);
        fileData = textscan(fileID,'%d%s%d\n');
        fclose(fileID);
        dataLengths = cellfun(@length,fileData);
        if length(unique(dataLengths))>1 || sum(dataLengths==0)>1
            fprintf('%s\n',['Scanlog is invalid for ' subject ', scan ' int2str(scan) '. Skipping this scan.']);
            continue;
        end;
        [acqs, expts, runs] = deal(fileData{1},fileData{2},fileData{3});
        
        for i = 1:length(acqs)
            
            % Skip run if it's not the specified experiment/run number.
            if ~strcmpi(expts{i},expt) || (~isempty(runList) && ~ismember(runs(i),runList))
                continue;
            end;
            
            runSuffix = ['-' int2str(runs(i))];
            
            outputDir = [analDir '/' expts{i} runSuffix inputSuffix outputSuffix '.perm'];
            outputDirOrig = outputDir;
            if exist(outputDir,'dir')
                if overwrite, system(['rm -rf ' outputDir]);
                else continue; end;
            end;
            mkdir(outputDir);
            outputMat = [outputDirOrig '/regressiondata.mat'];
            
            % Load three-column regressor files from regrDir
            % regressors3Col{r}: r-th regressor matrix
            regrFiles = regExpDir(regrDir,[subject '\.' expt '\.' int2str(runs(i)) '\.[0-9]+']);
            if isempty(regrFiles)
                fprintf('%s\n',['No regressor files found for ' subject ', expt ' ...
                    expts{i} runSuffix '. Skipping this run.']);
                continue;
            end;
            for r=1:length(regrFiles)
                [~,~,fileExt] = fileparts(regrFiles(r).name);
                regrNums(r) = str2num(strrep(fileExt,'.',''));
            end;
            [~,indSort] = sort(regrNums);
            regrFiles = regrFiles(indSort);
            cont = 0;
            for r=1:length(regrFiles)
                fileID = fopen([regrDir '/' regrFiles(r).name]);
                fileData = textscan(fileID,'%f%f%f\n');
                fclose(fileID);
                dataLengths = cellfun(@length,fileData);
                if length(unique(dataLengths))>1 || sum(dataLengths==0)>1
                    cont=1;
                    break;
                end;
                for col=1:3, regressors3Col{r}(:,col) = fileData{col}; end;
            end;
            if cont==1
                fprintf('%s\n',['Invalid regressor files for ' subject ', expt ' ...
                    expts{i} runSuffix '. Skipping this run.']);
                continue;
            end;
            numTaskRegr = length(regressors3Col);   % Number of task regressors
            regressors3ColOrig = regressors3Col;
            regrNames = {};
            for r=1:numTaskRegr
                regrNames{r} = ['Task' int2str(r)];
            end;
            
            % Check that size of conMat corresponds to # of regressors
            if size(conMat,2)~=numTaskRegr
                fprintf('%s\n',['Invalid conMat for ' subject ', expt ' expts{i} ...
                    runSuffix '. Skipping this run.']);
                continue;
            end;
            
            prepDir = [analDir '/' expts{i} runSuffix inputSuffix '.prep'];
            prepRegDir = [prepDir '/reg'];
            prepRegTargDir = [prepDir '/reg_target'];
            
            if ~exist(prepDir,'dir')
                fprintf('%s\n',['Preproc dir does not exist for ' subject ', ' ...
                    expts{i} runSuffix inputSuffix '. Skipping this run.']);
                continue;
            end;
            
            if isempty(maskName) || ~exist(funcMaskInit,'file')
                funcMaskInit = [prepRegTargDir '/mask.nii.gz'];
            end;
            funcMask = [outputDir '/mask.nii.gz'];
            if isempty(maskName) || ~exist(funcMaskInit,'file')
                system(['cp ' funcMaskInit ' ' funcMask]);
            else
                system([fslPrefix 'fslmaths ' funcMaskInit ' -mul ' ...
                    prepRegTargDir '/mask.nii.gz ' funcMask]);
            end;
            
            rfd = [prepDir '/raw_func.nii.gz'];
            badVols = load([prepDir '/art/badvols']);
            [~,numVols] = system([fslPrefix 'fslval ' rfd ' dim4']);
            numVols = str2num(strtrim(numVols));
            goodVols = setdiff(1:numVols,badVols);
            
            % Check TR of data for regressor definition
            tr = determineTR(rfd,forceTR,defaultTR,fslPrefix);
            
            if mod(tr/2,upsampledTR)~=0
                fprintf('%s\n',['TR/2 is not a multiple of upsampledTR for '...
                    subject ', expt ' expts{i} runSuffix '.  Skipping this run.']);
                continue;
            end;
            
            exptDuration = tr*numVols;
            
            % Load data/mask
            ffd = [prepRegTargDir '/filtered_func_data.nii.gz'];
            funcData = MRIread(ffd);
            dims = size(funcData.vol);
            maskData = MRIread(funcMask);
            maskInd = find(maskData.vol==1);
            funcMat = reshape(funcData.vol,[prod(dims(1:3)) dims(4)])';
            funcMat = funcMat(:,maskInd);
            funcMat = bsxfun(@minus,funcMat,mean(funcMat));
            newData = maskData;
            newData.vol = zeros(size(newData.vol));
            
            cont = 0;
            for iter=0:permIters
                
                % Redefine regressors3Col for new iteration
                regressors3Col = regressors3ColOrig;
                
                if iter>0
                    outputDir = [outputDirOrig '/perms/iter' int2str(iter)];
                    mkdir(outputDir);
                end;
                
                % Add rest blocks to regressors3Col if we want to permute them
                if permuteRest && iter>0
                    allRegr3Col = [];
                    for r=1:numTaskRegr
                        allRegr3Col = [allRegr3Col; regressors3Col{r}];
                    end;
                    regressors3Col{numTaskRegr+1} = defineRestBlocks(allRegr3Col,exptDuration,1);
                end;
                
                % Permute block order
                allRegr3Col = [];
                for r=1:length(regressors3Col)
                    allRegr3Col = [allRegr3Col; regressors3Col{r}];
                end;
                if iter>0
                    blockPerm{iter} = randperm(size(allRegr3Col,1));
                    allRegr3Col = allRegr3Col(blockPerm{iter},:);
                    blockCount = 0;
                    for r=1:length(regressors3Col)
                        regressors3Col{r} = allRegr3Col(blockCount+1:...
                            blockCount+size(regressors3Col{r},1),:);
                        blockCount = blockCount+size(regressors3Col{r},1);
                    end;
                end;
                
                if permuteRest && iter>0
                    regressors3Col = regressors3Col(1:end-1);
                end;
                
                % Define task regressors by convolving boxcar with double-gamma HRF.
                taskRegrMat = [];
                for r=1:numTaskRegr
                    taskRegrMat(:,r)= constructRegressor(regressors3Col{r},...
                        hrfType,numVols,tr,upsampledTR,subtractHalfTR);
                end;
                
                % Remove artifact time points and demean
                taskRegrMat = taskRegrMat(goodVols,:);
                taskRegrMat = bsxfun(@minus,taskRegrMat,mean(taskRegrMat));
                
                % Filter task regressors
                if tempFilt
                    [kernel,~] = fir1(filtOrder,2*tr*filtCutoff,'high');
                    for r = 1:size(taskRegrMat,2)
                        Xtmp = zeros(numVols,1);
                        Xtmp(goodVols) = taskRegrMat(:,r);
                        Xtmp = conv(Xtmp,kernel,'same');
                        taskRegrMat(:,r) = Xtmp(goodVols);
                    end;
                end;
                
                % On first iteration, define nuisance regressors and plot
                % regressor information.
                if iter==0
                    
                    nuisRegrMat = [];
                    
                    if useNuisLinear
                        % Linear trend regressor
                        linRegr = linspace(0,1,numVols);
                        linRegr = linRegr(goodVols);
                        linRegr = linRegr - mean(linRegr);
                        nuisRegrMat = [nuisRegrMat linRegr'];
                        regrNames{end+1} = 'NuisLinear';
                    end;
                    
                    if useNuisMotion
                        mcfFile = dir([prepDir '/mc/prefiltered_func_data*_mcf.par']);
                        moRegr = load([prepDir '/mc/' mcfFile(1).name]);
                        nuisRegrMat = [nuisRegrMat moRegr];
                        for r=1:size(moRegr,2)
                            regrNames{end+1} = ['NuisMotion' int2str(r)];
                        end;
                    end;
                    
                    % Add custom nuisance regressors
                    if ~isempty(customNuisRegr)
                        if size(customNuisRegr,1)~=numVols
                            fprintf('%s\n',['Custom nuisance regressor matrix must '...
                                'have first dimension equal to # of volumes in raw '...
                                'dataset.  Not using these regressors for ' subject ...
                                ', expt ' expts{i} runSuffix inputSuffix outputSuffix '.']);
                        end;
                        
                        nuisRegrMat(:,end+1:end+size(customNuisRegr,2)) = customNuisRegr(goodVols,:);
                        
                        for r=1:size(customNuisRegr,2)
                            regrNames{end+1} = ['NuisCustom' int2str(r)];
                        end;
                    end;
                    
                    % Extract WM/CSF PCA regressors.
                    if useNuisWMPCA || useNuisCSFPCA
                        
                        % Load unsmoothed functional data
                        pfdBET = [prepRegTargDir '/prefiltered_func_data_bet.nii.gz'];
                        funcDataTmp = MRIread(pfdBET);
                        funcMatTmp = reshape(funcDataTmp.vol,[prod(dims(1:3)) dims(4)])';
                        
                        noiseMat = [];  % Time point by voxel matrix, corresponding to noise volume
                        
                        if useNuisWMPCA
                            wmData = MRIread(wmMask);
                            wmInd = find(wmData.vol==1);
                            wmMat = funcMatTmp(:,wmInd);
                            noiseMat = [noiseMat wmMat];
                        end;
                        
                        if useNuisCSFPCA
                            csfData = MRIread(csfMask);
                            csfInd = find(csfData.vol==1);
                            csfMat = funcMatTmp(:,csfInd);
                            noiseMat = [noiseMat csfMat];
                        end;
                        
                        noiseMat = bsxfun(@minus,noiseMat,mean(noiseMat));  % Remove mean at each voxel
                        
                        % Orthogonalize w.r.t. existing nuisance regressors
                        % noiseMat = noiseMat - nuisRegrMat*inv(nuisRegrMat'*nuisRegrMat)*nuisRegrMat'*noiseMat;
                        
                        [pcaRegr,~] = svd(noiseMat,'econ');
                        nuisRegrMat(:,end+1:end+pcaOrder) = pcaRegr(:,1:pcaOrder);
                        
                        for r=1:pcaOrder
                            regrNames{end+1} = ['NuisPCA' int2str(r)];
                        end;
                        
                    end;
                    
                    % Demean nuisance regressors
                    nuisRegrMat = bsxfun(@minus,nuisRegrMat,mean(nuisRegrMat));
                    
                    % Filter nuisance regressors
                    if tempFilt
                        [kernel,~] = fir1(filtOrder,2*tr*filtCutoff,'high');
                        for r = 1:size(nuisRegrMat,2)
                            Xtmp = zeros(numVols,1);
                            Xtmp(goodVols) = nuisRegrMat(:,r);
                            Xtmp = conv(Xtmp,kernel,'same');
                            nuisRegrMat(:,r) = Xtmp(goodVols);
                        end;
                    end;
                    
                    % Define design matrix
                    X = [taskRegrMat nuisRegrMat];
                    
                    % Plot regressor information
                    
                    % 1) Plot all regressors
                    regrFig = figure('units','normalized','position',[.1 .1 .8 .8]);
                    if ~plotResults, set(regrFig,'visible','off'); end;
                    numRegrs = length(regrNames);
                    if numRegrs>15
                        maxRows = 10;
                    else
                        maxRows = 5;
                    end;
                    figRows = min(numRegrs,maxRows);
                    figCols = ceil(numRegrs/maxRows);
                    for r=1:numRegrs
                        subplot(figRows,figCols,r)
                        plot(X(:,r));
                        xlim([1 size(X,1)]);
                        if r==1
                            title([subject ', ' expts{i} runSuffix inputSuffix ...
                                outputSuffix ': ' regrNames{r}],'interpreter','none');
                        else
                            title(regrNames{r},'interpreter','none');
                        end;
                    end;
                    saveas(regrFig,[outputDir '/regr_plots.png']);
                    
                    % 2) Plot regressor correlation matrix, scale [-1 1]
                    corrFig = figure;
                    if ~plotResults, set(corrFig,'visible','off'); end;
                    imagesc(corr(X),[-1 1]); colorbar;
                    set(gca,'XTick',1:numRegrs,'YTick',1:numRegrs,'XTickLabel',...
                        [],'YTickLabel',regrNames);
                    title([subject ', ' expts{i} runSuffix inputSuffix outputSuffix ...
                        ': Regressor Correlations'],'interpreter','none');
                    saveas(corrFig,[outputDir '/regr_correlations.png']);
                    
                    % 3) Plot proportion of task regressor variance removed by nuisance regressors
                    varianceFig = figure;
                    if ~plotResults, set(varianceFig,'visible','off'); end;
                    taskRegrMatProj = taskRegrMat-nuisRegrMat*inv(nuisRegrMat'*nuisRegrMat)*nuisRegrMat'*taskRegrMat;
                    varRemoved = ones(1,numTaskRegr)-sum(taskRegrMatProj.^2)./sum(taskRegrMat.^2);
                    bar(varRemoved);
                    ylim([0 1]);
                    set(gca,'XTick',1:numTaskRegr,'XTickLabel',regrNames(1:numTaskRegr));
                    title([subject ', ' expts{i} runSuffix inputSuffix outputSuffix ...
                        ': Variance Explained by Nuisance Regrs'],'interpreter','none');
                    saveas(varianceFig,[outputDir '/regr_variance_removed.png']);
                    
                    if ~plotResults, close(regrFig,corrFig,varianceFig); end;
                    
                end;
                
                % Design matrix X
                X = [taskRegrMat nuisRegrMat];
                
                if size(X,1)<size(X,2) || rank(X)<size(X,2)
                    cont = 1;
                    errorMsg = ['Design matrix is rank-deficient for ' subject ...
                            ', expt ' expts{i} runSuffix inputSuffix outputSuffix '. '...
                            'Skipping this run.'];
                    break;
                end;
                
                % Compute beta values and residuals
                V = inv(X'*X);
                betas = V*X'*funcMat;
                if iter==0
                    resids = funcMat-X*betas;
                    errorDOF = size(X,1)-size(X,2)-1;
                    errorVar = sum(resids.^2)/errorDOF;
                    rSquareds = ones(1,length(maskInd))-sum(resids.^2)./sum(funcMat.^2);
                end;
                
                % Compute contrasts and statistics
                contrasts = conMat*betas(1:numTaskRegr,:);
                conVarBase = diag(conMat*V(1:numTaskRegr,1:numTaskRegr)*conMat');  % Estimated contrast variance, up to error variance term
                if iter==0
                    conVars = repmat(errorVar,[size(contrasts,1) 1]).*repmat(conVarBase,[1 size(contrasts,2)]);
                    tStats = contrasts./sqrt(conVars);
                    zStats = zeros(size(tStats));
                    zStats(tStats>0) = -icdf('norm',tcdf(-tStats(tStats>0),errorDOF),0,1);
                    zStats(tStats<=0) = icdf('norm',tcdf(tStats(tStats<=0),errorDOF),0,1); % Prevents p-values near 1 from rounding to 1
                end;
                
                % Save regression information
                if iter==0
                    save(outputMat,'regressors3Col','X','V','betas','conMat',...
                        'contrasts','conVarBase','resids','errorDOF','errorVar',...
                        'conVars','tStats','zStats','rSquareds','regrNames','tr');
                else
                    conVarBasePerm{iter} = conVarBase;
                end;
                
                % Write contrast maps
                for c=1:numCons
                    outputPath = [outputDir '/cope' int2str(c) '.nii.gz'];
                    newData.vol(maskInd) = contrasts(c,:)';
                    MRIwrite(newData,outputPath);
                end;
                
                % For unpermuted analysis, write pe/var/zstat/etc images
                if iter==0
                    for r=1:size(betas,1)
                        outputPath = [outputDir '/pe' int2str(r) '.nii.gz'];
                        newData.vol(maskInd) = betas(r,:)';
                        MRIwrite(newData,outputPath);
                    end;
                    for c=1:numCons
                        outputPath = [outputDir '/zstat' int2str(c) '_ols.nii.gz'];
                        newData.vol(maskInd) = zStats(c,:)';
                        MRIwrite(newData,outputPath);
                        outputPath = [outputDir '/varcope' int2str(c) '_ols.nii.gz'];
                        newData.vol(maskInd) = conVars(c,:)';
                        MRIwrite(newData,outputPath);
                    end;
                    
                    outputPath = [outputDir '/sigmasquareds_ols.nii.gz'];
                    newData.vol(maskInd) = errorVar';
                    MRIwrite(newData,outputPath);
                    outputPath = [outputDir '/rSquareds_ols.nii.gz'];
                    newData.vol(maskInd) = rSquareds;
                    MRIwrite(newData,outputPath);
                    
                    if writeResiduals && iter==0
                        newData = funcData;
                        newData.vol = zeros(size(newData.vol));
                        for t=1:size(resids,1)
                            tmp = zeros(dims(1:3));
                            tmp(maskInd) = resids(t,:);
                            newData.vol(:,:,:,t) = tmp;
                        end;
                        outputPath = [outputDir '/res4d.nii.gz'];
                        MRIwrite(newData,outputPath);
                    end;
                    
                    system(['cp ' prepRegDir '/target_func.nii.gz ' outputDir '/target_func.nii.gz']);
                end;
                
                if mod(iter,10)==0
                    fprintf('%s\n',[subject ', ' expts{i} runSuffix...
                        inputSuffix outputSuffix ', iter ' int2str(iter)]);
                end;
                
            end;
            if cont
                fprintf('%s\n\n',errorMsg);
                continue;
            end;
            save(outputMat,'-append','blockPerm','conVarBasePerm');
            
            fprintf('%s\n\n',['Finished 1st-level analysis for subject ' subject ...
                ', expt ' expts{i} runSuffix inputSuffix outputSuffix '.']);
        end;
    end;
end;

end