
% model2ndPerm(studyDir,subjects,varargin)
% 
%%% CURRENTLY UNDERGOING MAJOR OVERHAUL
%%% CONVERTING TO BIDS NAMING CONVENTION, etc
%
%
%
%
% 
% Step 2 of a two-step process (modelPerm, model2ndPerm) to perform a
% General Linear Model based analysis of fMRI data, computing statistics
% using a permutation test.  Step 2 combines results across runs, and
% computes permutation-based statistics at the level of both individual
% voxels and clusters, using a cluster extent threshold to correct for 
% multiple comparisons across voxels.  Step 2 can also be used to compute 
% statistics for individual runs, by specifying one run with the runList 
% argument.  Function can be run multiple times with different voxel- and
% cluster-level thresholds, or different conList arguments.  Should be run
% after modelPerm.  Loops across subjects and runs by default.
%
% Example usage: model2ndPerm('/pathto/studyDirectory','SUB*','expt','FaceLoc')
%
% Critical output files (in *.gperm/cope* directory):
% - cope*.nii.gz: cross-run contrast image
% - cope*_standard.nii.gz: contrast image registered to MNI space, which
%   can be used as input to group analysis
% - zstat*.nii.gz: permutation-based voxelwise statistics, transformed to
%   z-scores, with no cluster correction
% - thresh_zstat*.nii.gz: statistical map after applying cluster extent
%   threshold to correct for multiple comparisons
% - cluster_permutation_data*.mat: information about the cluster extent
%   null distribution for a given voxelwise threshold
% - perms: cross-run contrast and z-statistic maps for each permutation
%
% Required arguments:
% - studyDir (string): path to study directory
% - subjects (string or cell array of strings): list of subjects to
%       analyze. Can use asterisk-based regular expression.  Examples:
%       {'SUB01','SUB02'} or 'SUB*'.
% 
% Variable arguments:
% - overwrite (boolean; default=0): whether to overwrite files that have
%       already been written by this function.
% - expt (string): name of experiment to process.
% - runList (vector of integers in (0,Inf)): list of runs to process.
% - conList (vector of integers in (0,Inf)): list of contrasts to process,
%   indexed by the order of conMat argument to modelPerm.
% - inputSuffix (string): suffix for input .perm directories.
% - outputSuffix (string): suffix for .gperm directory. Can be used to run
%       the script multiple times with different options.
% - voxThresh (scalar in (0,1); default=.001): voxelwise p-threshold, i.e.
%   desired Type I error rate.
% - clustThresh (scalar in (0,1); defualt=.05): clusterwise p-thresold,
%   i.e. desired Type I error rate.
% - clustCorrect (boolean; default=1): whether to run cluster-based
%   thresholding.
% - looRun (integer in [0,Inf); default=0): run to leave out on LOO mode 
%   (0 = no run left out).
% - targetName (string; default='target_func'): name of image in func 
%       directory that was used for registration by preprocMRI, without 
%       file extension. If a custom image was specified for preprocMRI,
%       the same image should be used here.

function model2ndPerm(studyDir,subjects,varargin)

% Check system configuration
fpp.util.checkConfig;

if ~exist(studyDir,'dir')
    fprintf('%s\n',['ERROR: Study directory ' studyDir ' does not exist.']);
    return;
end

overwrite = 0;              % Whether to overwrite output
expt = [];                  % Experiment to analyze.
runList = [];               % Runs to analyze.
conList = [];               % Contrasts to analyze (indexed based on conMat used in modelPerm)
inputSuffix = '';           % Suffix of input .perm dirs
outputSuffix = '';          % New suffix for output dir
voxThresh = .001;           % Voxelwise threshold
clustThresh = .05;          % Clusterwise threshold
clustCorrect = 1;           % Whether to run cluster-based thresholding
looRun = 0;                 % Run to leave out on LOO mode (0 = no run left out)
targetName = 'target_func'; % Name of target_func image used for registration
charLimit = 100000;         % Character limit for bash commands

% Convert list of subjects from "dir" inputs to actual names
subjectsNew = {};
if ~iscell(subjects), subjects = {subjects}; end
for s=1:length(subjects)
    if ~ischar(subjects{s})
        fprintf('%s\n','ERROR: subject argument must be a string or cell array of strings.');
        return;
    end
    if strfind(subjects{s},'*')
        subjectsTmp = dir([studyDir '/' subjects{s}]);
        for j=1:length(subjectsTmp)
            subjectsNew{end+1} = subjectsTmp(j).name;
        end
    else
        subjectsNew{end+1} = subjects{s};
    end
end
subjects = subjectsNew;

% Edit variable arguments.  Note: optInputs checks for proper input.
varArgList = {'overwrite','looRun','clustCorrect','expt','runList','conList',...
    'inputSuffix','outputSuffix','voxThresh','clustThresh'};
for i=1:length(varArgList)
    argVal = fpp.util.optInputs(varargin,varArgList{i});
    if ~isempty(argVal)
        eval([varArgList{i} ' = argVal;']);
    end
end

if ~isempty(runList) && looRun>0 && ~ismember(looRun,runList)
    fprintf('%s\n','ERROR: looRun must be in runList when both are specified.');
    return;
end

for s = 1:length(subjects)
    
    subject = subjects{s};
    subjDir = [studyDir '/' subject];
    analDir = [subjDir '/analysis'];
    regDir = [subjDir '/reg'];
    regTargDir = [subjDir '/reg/' targetName];
    
    scanlogFiles = fpp.util.regExpDir(subjDir,'.*scanlo.*([^~]$)');
    
    if exist(subjDir,'dir')==0
        fprintf('%s\n',['Subject directory does not exist for ' subject '.']);
        continue;
    elseif isempty(scanlogFiles)
        fprintf('%s\n',['Missing scan log file for subject ' subject '.']);
        continue;
    elseif exist(analDir,'dir')==0
        fprintf('%s\n',['Analysis directory does not exist for ' subject '.']);
        continue;
    elseif exist(regDir,'dir')==0
        fprintf('%s\n',['Registration directory does not exist for ' subject '.']);
        continue;
    elseif exist(analDir,'dir')==0
        fprintf('%s\n',['Target registration directory does not exist for ' subject '.']);
        continue;
    end
    
    exptTypes = {};    % Cell array of experiments run
    exptNums = [];     % For each acquisition, which expt # (indexed from exptTypes)
    runNums = [];      % For each acquisition, which run #
    permDirsAll = {};     % For each acquisition, .perm directory path
    
    % Get list of experiments, etc
    cont = 0;
    for scan = 1:length(scanlogFiles)
        fileID = fopen([subjDir '/' scanlogFiles(scan).name]);
        fileData = textscan(fileID,'%d%s%d\n');
        fclose(fileID);
        dataLengths = cellfun(@length,fileData);
        if length(unique(dataLengths))>1 || sum(dataLengths==0)>1
            cont = 1;
            fprintf('%s\n',['Scanlog is invalid for ' subject ', scan ' int2str(scan)...
                '. Skipping this subject.']);
            break;
        end
        [acqs, expts, runs] = deal(fileData{1},fileData{2},fileData{3});
        for i = 1:length(acqs)
            if ~ismember(expts{i},exptTypes)
                exptTypes{end+1} = expts{i};
            end
            exptNums(end+1) = find(ismember(exptTypes,expts{i}));
            runNums(end+1) = runs(i);
            runSuffix = ['-' int2str(runs(i))];
            permDirsAll{end+1} = [analDir '/' expts{i} runSuffix inputSuffix '.perm'];
        end
    end
    if cont, continue; end
    
    for e = 1:length(exptTypes)
        
        exptType = exptTypes{e};
        
        % Skip run if it's anatomical, DTI, or rest.
        if sum(regexpi(exptType,'anat'))>0 || sum(regexpi(exptType,'dti'))>0 ...
                || sum(regexpi(exptType,'rest'))>0
            continue;
        end
        
        % Skip run if it's not the specified experiment/run number.
        if ~isempty(expt) && ~strcmpi(exptType,expt)
            continue;
        end
        
        extantRuns = runNums(exptNums==e);
        if ~isempty(runList) && sum(ismember(extantRuns,runList))~=length(runList)
            fprintf('%s\n',['Subject ' subject ' does not have specified runs of expt '...
                exptType '. Skipping this experiment.']);
            continue;
        elseif looRun>0 && ~ismember(looRun,extantRuns)
            fprintf('%s\n',['Subject ' subject ' does not have looRun for expt '...
                exptType '. Skipping this experiment.']);
            continue;
        elseif length(unique(extantRuns))~=length(extantRuns)
            fprintf('%s\n',['Subject ' subject ' has redundant scanlog run numbers for expt '...
                exptType '. Skipping this experiment.']);
            continue;
        end
        
        permDirs = permDirsAll(exptNums==e);   % List of extant .perm dirs for this expt
        
        runSuffix = '-2ndLevel';
        
        if isempty(runList)
            runList = extantRuns;
        end
        if looRun>0
            runList = setdiff(runList,looRun);
            runSuffix = [runSuffix '-LOO' int2str(looRun)];
        end
        permDirs = permDirs(ismember(extantRuns,runList));
        
        numRuns = length(permDirs);
        
        % Check for existence of input directories.
        cont = 0;
        for d = 1:numRuns
            if ~exist(permDirs{d},'dir')
                fprintf('%s\n',['Subject ' subject ' is missing input directory'...
                    ' for expt ' exptType '. Skipping this experiment.']);
                cont=1;
                break;
            end
        end
        if cont, continue; end
        
        for r=1:numRuns
        	rData{r} = load([permDirs{r} '/regressiondata.mat']);
            permItersVec(r) = length(rData{r}.conVarBasePerm);
        end
        if length(unique(permItersVec))>1
            fprintf('%s\n',['Input perm directories have inconsistent # of iterations '...
                'for ' subject ', ' exptType '. Skipping this experiment.']);
            continue;
        end
        permIters = permItersVec(1);
        
        numCons = length(dir([permDirs{1} '/cope*.nii.gz']));
        
        if isempty(conList)
            conList = 1:numCons;
        else
            conList = conList(ismember(conList,1:numCons));
            if isempty(conList)
                fprintf('%s\n',['Issue with conList for ' subject ', ' ...
                    exptType '. Skipping this experiment..']);
                continue;
            end
        end
        
        gpermDir = [analDir '/' exptType runSuffix inputSuffix outputSuffix '.gperm'];
        if exist(gpermDir,'dir') && overwrite, fpp.util.system(['rm -rf ' gpermDir]); end
        if ~exist(gpermDir,'dir'), mkdir(gpermDir); end
        
        targetNew = [gpermDir '/target_func.nii.gz'];
        if ~exist(targetNew,'file')
            fpp.util.system(['cp ' permDirs{1} '/target_func.nii.gz ' targetNew]);
        end
        
        for c=conList
            
            outputDir = [gpermDir '/cope' int2str(c)];
            permsDir = [outputDir '/perms'];
            if ~exist(permsDir,'dir')
                mkdir(permsDir);
            end
            
            outputZStat = [outputDir '/zstat' int2str(c) '.nii.gz'];
            outputThreshZStat = [outputDir '/thresh_zstat' int2str(c)  '_P' ...
                num2str(voxThresh) 'P' num2str(clustThresh) '.nii.gz'];
            if exist(outputThreshZStat,'file'), continue; end
            
            % Check for existence of input images
            cont = 0;
            for r=1:numRuns
                for iter=1:permIters
                    inputCope = [permDirs{r} '/perms/iter' int2str(iter) ...
                        '/cope' int2str(c) '.nii.gz'];
                    if ~exist(inputCope,'file')
                        fprintf('%s\n',['Not all input cope files exist for ' ...
                            subject ', ' exptType runSuffix inputSuffix ...
                            ', cope ' int2str(c) '. Skipping this experiment.']);
                        cont = 1;
                        break;
                    end
                end
                if cont, break; end
            end
            if cont, continue; end
            
            fprintf('%s\n\n',['Running 2nd-level analysis for subject ' subject ', expt ' ...
                exptType runSuffix inputSuffix outputSuffix ', contrast ' int2str(c) '...']);
            
            % Initial voxelwise permutation test
            if ~exist(outputZStat,'file')
                
                % Average input cope images to compute overall
                % permutation statistic image
                outputCope = [outputDir '/cope' int2str(c) '.nii.gz'];
                if ~exist(outputCope,'file')
                    cmd = ['fslmaths '];
                    sumWeights = 0;    % Sum of weighting factors for cope images
                    for r=1:numRuns
                        inputCope = [permDirs{r} '/cope' int2str(c) '.nii.gz'];
                        tmpPath{r} = [permDirs{r} '/tmp_weighted_cope.nii.gz'];
                        fpp.util.system(['fslmaths ' inputCope ' -mul ' ...
                            num2str(1/rData{r}.conVarBase(c)) ' ' tmpPath{r}]);
                        sumWeights = sumWeights+1/rData{r}.conVarBase(c);
                        if r==1
                            cmd = [cmd tmpPath{r} ' '];
                        else
                            cmd = [cmd '-add ' tmpPath{r} ' '];
                        end
                    end
                    cmd = [cmd '-div ' num2str(sumWeights) ' ' outputCope];
                    fpp.util.system(cmd);
                    for r=1:numRuns
                        fpp.util.system(['rm -rf ' tmpPath{r}]);
                    end
                end
                
                % Register cope image to standard space
                target2Highres = [regTargDir '/target_func2highres.mat'];
                copeStd = [outputDir '/cope' int2str(c) '_standard.nii.gz'];
                standard = [getenv('FSL_DIR') '/data/standard/MNI152_T1_2mm_brain.nii.gz'];
                highres2StdWarp = [regDir '/highres2standard_warp.nii.gz'];
                fpp.util.system(['applywarp --ref=' standard ' --in=' outputCope ' --out=' ...
                    copeStd ' --warp=' highres2StdWarp ' --premat=' target2Highres]);
                
                % Do this for many iterations of lists of input permutations
                for iter=1:permIters
                    outputCopePerm = [permsDir '/cope' int2str(c) '_iter' int2str(iter) '.nii.gz'];
                    if ~exist(outputCopePerm,'file')
                        cmd = ['fslmaths '];
                        sumWeights = 0;    % Sum of weighting factors for cope images
                        for r=1:numRuns
                            inputCope = [permDirs{r} '/perms/iter' int2str(iter) '/cope' int2str(c) '.nii.gz'];
                            tmpPath{r} = [permDirs{r} '/perms/iter' int2str(iter) '/tmp_weighted_cope.nii.gz'];
                            fpp.util.system(['fslmaths ' inputCope ' -mul ' ...
                                num2str(1/rData{r}.conVarBasePerm{iter}(c)) ' ' tmpPath{r}]);
                            sumWeights = sumWeights+1/rData{r}.conVarBasePerm{iter}(c);
                            if r==1
                                cmd = [cmd tmpPath{r} ' '];
                            else
                                cmd = [cmd '-add ' tmpPath{r} ' '];
                            end
                        end
                        cmd = [cmd '-div ' num2str(sumWeights) ' ' outputCopePerm];
                        fpp.util.system(cmd);
                        for r=1:numRuns
                            fpp.util.system(['rm -rf ' tmpPath{r}]);
                        end
                    end
                    
                    if mod(iter,10)==0
                        fprintf('%s\n',['Voxelwise, ' subject ', ' exptType runSuffix...
                            inputSuffix outputSuffix ', con ' int2str(c) ', iter ' int2str(iter)]);
                    end
                end
                
                % Merge permstat files
                copeConcat = [permsDir '/cope' int2str(c) '_concat.nii.gz'];
                mergeCmd = ['fslmerge -t ' copeConcat ' '];
                
                % Split up merge command to avoid Linux character limits
                cLength = length([permsDir '/iter' int2str(permIters) ...
                    '/cope' int2str(c) '.nii.gz ']);
                mLength = length(mergeCmd);
                copesPerMerge = floor((charLimit-mLength)/cLength);
                mergeIters = ceil(permIters/copesPerMerge);
                
                for c2=1:mergeIters
                    copeConcats{c2} = [permsDir '/cope' int2str(c) ...
                        '_concat_miter' int2str(c2) '.nii.gz'];
                    mergeCmd = [mergeCmd ' ' copeConcats{c2}];
                    mergeCmd2 = ['fslmerge -t ' copeConcats{c2} ' '];
                    for j=(1+(c2-1)*copesPerMerge):min(permIters,c2*copesPerMerge)
                        outputCopePerm = [permsDir '/cope' int2str(c) '_iter' int2str(j) '.nii.gz'];
                        mergeCmd2 = [mergeCmd2 ' ' outputCopePerm];
                    end
                    fpp.util.system(mergeCmd2);
                end
                mergeCmd = [mergeCmd '; rm -rf ' permsDir '/cope' ...
                    int2str(c) '_concat_miter*'];
                fpp.util.system(mergeCmd);
                
                % Compute zstat file from perm stat files, by fitting a Gaussian
                copeMean = [permsDir '/cope' int2str(c) '_mean.nii.gz'];
                copeStd = [permsDir '/cope' int2str(c) '_std.nii.gz'];
                fpp.util.system(['fslmaths ' copeConcat ' -Tmean ' copeMean]);
                fpp.util.system(['fslmaths ' copeConcat ' -Tstd ' copeStd]);
                fpp.util.system(['fslmaths ' outputCope ' -sub ' copeMean ...
                    ' -div ' copeStd ' ' outputZStat]);
                
            end
            
            if ~clustCorrect
                fprintf('%s\n\n',['Finished 2nd-level analysis for subject ' subject ', expt ' ...
                    exptType runSuffix inputSuffix outputSuffix ', contrast ' int2str(c) '.']);
                continue;
            end
            
            zCutoff = icdf('norm',1-voxThresh,0,1);
            cDataPath = [outputDir '/cluster_permutation_data_P' num2str(voxThresh) '.mat'];
            testFiles = dir([outputDir '/thresh_zstat' int2str(c) '_P' num2str(voxThresh) '*.nii.gz']);
            copeMean = [permsDir '/cope' int2str(c) '_mean.nii.gz'];
            copeStd = [permsDir '/cope' int2str(c) '_std.nii.gz'];
            
            if isempty(testFiles) || ~exist(cDataPath,'file')
                
                % List of cluster sizes
                clustSizes = [];
                
                for iter=1:permIters
                    outputCopePerm = [permsDir '/cope' int2str(c) '_iter' int2str(iter) '.nii.gz'];
                    zStatPermPath = [permsDir '/zstat' int2str(c) '_iter' int2str(iter) '.nii.gz'];
                    fpp.util.system(['fslmaths ' outputCopePerm ' -sub ' copeMean ...
                        ' -div ' copeStd ' ' zStatPermPath]);
                    [~,clustInfo] = fpp.util.system(['cluster -i ' zStatPermPath ' -t ' num2str(zCutoff)]);
                    clustInfo = regexp(clustInfo, '[\f\n\r]', 'split');
                    clustInfo{2} = strread(clustInfo{2});
                    
                    if isempty(clustInfo{2})
                        clustSizes = [clustSizes 0];
                    else
                        clustSizes = [clustSizes clustInfo{2}(2)];
                    end
                    
                    if mod(iter,10)==0
                        fprintf('%s\n',['Clusterwise, ' subject ', ' exptType runSuffix...
                            inputSuffix outputSuffix ', con ' int2str(c) ', iter ' int2str(iter)]);
                    end
                end
                
                % Determine maximum cluster size to consider, based on mean cope
                [~,maxClustSize] = fpp.util.system(['fslstats ' copeMean ' -V']);
                maxClustSize = str2num(strtrim(maxClustSize));
                maxClustSize = maxClustSize(1);
                
                % Compute cluster distribution and CDF
                clustHist = hist(clustSizes,1:maxClustSize);
                clustDist = clustHist/sum(clustHist);
                clustCDF = cumsum(clustDist);
                
                save(cDataPath,'clustSizes','clustHist','clustDist','clustCDF');
                
            else
                load(cDataPath);
            end
            
            % Keep clusters with size (voxel count) >=clustCutoff
            clustCutoff = find(clustCDF > 1-clustThresh,1,'first')+1;
            fprintf('\n%s\n',['Cluster size threshold: ' int2str(clustCutoff) ' voxels']);
            
            zCutoff = icdf('norm',1-voxThresh,0,1);
            tmpPath{1} = [outputDir '/tmp_cluster_thresholding.nii.gz'];
            fpp.util.system(['cluster -i ' outputZStat ' -t ' num2str(zCutoff) ' --othresh=' ...
                outputThreshZStat ' --osize=' tmpPath{1} ' --no_table']);
            fpp.util.system(['fslmaths ' tmpPath{1} ' -thr ' int2str(clustCutoff) ' -bin ' tmpPath{1}]);
            fpp.util.system(['fslmaths ' outputThreshZStat ' -mul ' tmpPath{1} ' ' outputThreshZStat]);
            fpp.util.system(['rm -rf ' tmpPath{1}]);
            
            fprintf('%s\n\n',['Finished 2nd-level analysis for subject ' subject ', expt ' ...
                exptType runSuffix inputSuffix outputSuffix ', contrast ' int2str(c) '.']);
            
        end
    end
end

end
