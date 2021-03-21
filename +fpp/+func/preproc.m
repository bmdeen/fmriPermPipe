%
% fpp.func.preproc(inputPaths,outputDir,varargin)
% 
% Preprocesses a single (multi-echo) fMRI dataset, including including
% motion parameter estimation, despiking, slice timing correction, one-shot
% motion and distoration correction and template registration, TEDANA 
% multi-echo ICA denoising and optimal echo combination, intensity
% normalization, and optional spatial/temporal filtering. Should be run
% after anatomical scripts, func.defineTemplate, and func.register.
% 
% Example usage:
% fpp.func.preproc({'/pathToData/sub-01_task-faceloc_run-01_echo-1_bold.nii.gz',...
%       '/pathToData/sub-01_task-faceloc_run-01_echo-2_bold.nii.gz'},'/bidsDir/sub-01/');
% 
% Required arguments:
% - inputPaths (string or cell array of strings): paths to input fMRI data,
%       one path for each echo.
% - outputDir (string): path to output directory (subject/session dir, in
%       BIDS filesystem)
% 
% Variable arguments:
% - overwrite (boolean; default=0): whether to overwrite files that have
%       already been written by this function.
% - fwhm (scalar in [0,Inf); default=2.4): FWHM of spatial smoothing kernel 
%       (mm).
% - smThresh (scalar in (0,Inf)): brightness threshold for SUSAN-based
%       smoothing.
% - disdaqs = 3 (integer in [0,Inf); default=3): number of disdaq volumes
%       (at beginning of run) to remove. [NOT CURRENTLY IMPLEMENTED]
% - transCutoff (scalar in [0,Inf]; default=.5): artifact detection cutoff
%       for total translation (mm) between volumes.
% - rotCutoff (scalar in [0,Inf]; default=.5): artifact detection cutoff
%   for total rotation (degrees) between volumes.
% - tptsAfter (integer in [0,Inf); default=0): remove this many time points
%       after pairs of volumes with motion.
% - teVals (vector of values in (0,Inf)): TE values (ms) for multi-echo 
%       data (default: read from json)
% - sliceTimes (vector of values in (0,Inf)): vector of slice acquisition
%       times (s) relative to volume acq onset. (default: read from json)
% - tempFilt (boolean; default=0): whether to perform high-pass temporal
%       filtering.
% - filtCutoff (scalar in (0,Inf); default=.01: highpass filter cutoff (Hz)
% - filtOrder (integer in (0,Inf); default=60): Hamming-window FIR filter
%       order.
% - filtType (string): high, low, or bandpass (default: high if filtCutoff
%       has one element, bandpass if it has two)
% - genTSNR (boolean; default=0): whether to generate TSNR maps.
% - plotResults (boolean; default=1): whether to display result plots.
% - useSTC (boolean; default=1): whether to correct for slice timing.
% - useTedana (boolean; default=1): whether to use TEDANA denoising, for
%       multi-echo data.
% - funcTemplatePath (string): path to functional template image to
%       register to (default: determined from input name)
% - funcTemplateSpace (string): space of functional template, needed if
%       using non-standard funcTemplatePath that lacks a BIDS space entity
% - useTaskTemplate (boolean): whether to use a task-specific template
% - undistort (boolean): Whether to distortion-correct functional data
%       using blip-up/blip-down method
% - spinEchoPath (string): path to spin echo image matched in phase-encode
%       direction to input fMRI data
% - topupWarpPath (string): path to undistortion warp image produced by
%       topup
% - topupJacobianPath (string): path to undistortion warp jacobian image
%       produced by topup
% - deleteMidprep (boolean; default=1): whether to delete midprep images
%

function preproc(inputPaths,outputDir,varargin)

% TODO IMMEDIATELY:
% - Add JSON files for additional TEDANA/ts2map outputs
% - At the end of the script, resample to several spaces if desired:
%   - 2mm-res individual space
%   - 2mm-res MNI152Nlin6Asym space
%   - CIFTI: 32k fsLR surface, 2mm individual subcortical
%   - CIFTI: 32k fsLR surface, 2mm MNI152Nlin6Asym subcortical space
% - Run steps 8-10 on optcomb (non-TEDANA) output as well
%
% TODO NEXT:
% - Add spatial smoothing (susan- or wb-command-based)
% - If deleteMidprep==1, remove "Sources" field of final outputs
% - Add suffix option for desc. Need to edit inputNameGeneric, and edit
%   desc of output TEDANA folder (defitions both within preproc script, and
%   tedana script - new tedanaDir optional argument for tedana script?).
%   Also need to edit any other name that has desc removed - e.g. confounds
%   file (fpp.func.preproc.estimateHeadMotion and preproc step 1); xfm
%   matrices (xfmMocoTarget2FuncTemplate and inverse in step 5)
%
% TODO EVENTUALLY:
% - Generate figures on preproc results
% - Method to deal with disdaqs, taking into account BIDS info, adding them
%   as artifact time points, and NonSteadyStateOutlier00, in confound
%   regressors. Based on NumberOfVolumesDiscardedByUser
%   and NumberOfVolumesDiscardedByScanner json fields.
% 
%
%
% BIDS reference info:
%
% SpatialRef options:
%   Standard spaces:
%   fsaverage, fsLR, MNI152NLin6ASym (FSL/HCP), MNI152NLin2009cAsym
%   (fMRIPrep)
%   Individual spaces: fsnative, individual, session, task
%
% Potential BIDS confound TSV files:
% - desc-confounds_regressors.tsv - all regressors!
% - motion.tsv - just motion regressors
% - physreg.tsv - physiological regressors
% - outliers.tsv - artifact TPs (alt: ArtifactTimePoint00 in regressors.tsv)
% - mixing.tsv / components.tsv

% Check system configuration
fpp.util.checkConfig;

% Basic parameters
overwrite = 0;                  % Whether to overwrite output

% General preproc parameters
fwhm = 0;                       % FWHM of spatial smoothing kernel (mm)
smThresh = [];                  % SUSAN brightness threshold
funcTemplatePath = '';          % Functional template image to register to
funcTemplateSpace = '';         % Space of func template, required if using custom funcTemplatePath with BIDs space entity
useTaskTemplate = 0;            % Use a task/acq-specific functional template, rather than a session template
newMedian = 10000;              % Target value for intensity normalization (median of data in brain mask is set to this)
useTedana = 1;                  % Whether to use tedana-based denoising for multi-echo data
useSTC = 1;                     % Whether to use slice timing correction
useDespike = 1;                 % Whether to apply AFNI's 3dDepike before STC
sliceTimes = [];                % Slice time (s) relative to onset of volume acquisition (Default: read from JSON file)
deleteMidprep = 1;              % Whether to delete midprep images

% Undistortion parameters
undistort = 1;                  % Whether to distortion-correct functional data using blip-up/blip-down method
spinEchoPath = '';              % Spin echo path with PE dir matched to functional, if not using default (determined from json metadata)
topupWarpPath = '';             % Undistortion warp path, if not using default (determined from json metadata)
topupJacobianPath = '';         % Undistortion warp jacobian path, if not using default (determined from json metadata)

% Artifact detection parameters
disdaqs = 0;                    % Number of disdaq volumes (at beginning of run) to remove. [NOT CURRENTLY IMPLEMENTED]
transCutoff = .5;               % Artifact detection cutoff: total translation (mm)
rotCutoff = .5;                 % Artifact detection cutoff: total rotation (degrees)
tptsAfter = 0;                  % Remove this many time points after pairs of volumes with motion

% Multi-echo parameters
multiEcho = 0;                  % Whether data is multi-echo
teVals = [];                    % TE values (ms; Default: read from JSON file)
echoForMoCorr = [];             % Which echo to use to estimate motion parameters

% Temporal filtering parameters
tempFilt = 0;                   % Whether to perform temporal filtering
filtCutoff = .01;               % Highpass filter cutoff (Hz)
filtOrder = 60;                 % Hanning-window FIR filter order
filtType = [];                  % Filter type (high, low, bandpass, or bandstop)

% Data generation parameters
genTSNR = 1;                    % Whether to generate TSNR maps
plotResults = 1;                % Whether to display result plots

pathsToDelete = {};             % Midprep paths to delete at the end
midprepIntro = 'Partially preprocessed data generated by fmriPermPipe.';    % Description intro sentence for midprep data
preprocIntro = 'Preprocessed data generated by fmriPermPipe.';              % Description intro sentence for preproc data
steps = {};                     % Processing step descriptions

% Edit variable arguments.  Note: optInputs checks for proper input.
varArgList = {'overwrite','funcTemplatePath','funcTemplateSpace','fwhm',...
    'tempFilt','teVals','transCutoff','sliceTimes','useTaskTemplate',...
    'rotCutoff','tptsAfter','stdCutoff','disdaqs','genTSNR','plotResults',...
    'faValue','smThresh','useSTC','useTedana','echoForMoCorr','undistort',...
    'useDespike','topupWarpPath','topupJacobianPath','spinEchoPath',...
    'deleteMidprep'};
for i=1:length(varArgList)
    argVal = fpp.util.optInputs(varargin,varArgList{i});
    if ~isempty(argVal)
        eval([varArgList{i} ' = argVal;']);
    end
end

% Check if input is multi-echo
if ischar(inputPaths), inputPaths = {inputPaths};
elseif length(inputPaths)>1, multiEcho = 1; end

% Check echo time info for multi-echo data
if multiEcho
    if isempty(teVals)
        teVals = fpp.util.checkMRIProperty('TE',inputPaths{1});
        if isempty(teVals)
            fprintf('%s\n','ERROR: TE values are not specified and cannot be read from json metadata.');
            return;
        end
    end
    if size(teVals,1)<size(teVals,2), teVals = teVals'; end
    nEchoes=length(inputPaths);
    if isempty(echoForMoCorr)
        echoForMoCorr = 2;
    end
else
    nEchoes = 1;
    useTedana = 0;
    echoForMoCorr = 1;
end

% Check if inputs exist
for e=1:length(inputPaths)
    if exist(inputPaths{e},'file')==0
        fprintf('%s\n\n',['ERROR: Input path ' inputPaths{e} ' does not exist.']);
        return;
    end
end

% Define input name, output directories
[~,inputName,~] = fpp.util.fileParts(inputPaths{1});
inputNameGeneric = strrep(fpp.bids.changeName(inputName,'echo',[]),'_bold','');
if strcmp(outputDir(end),'/'), outputDir = outputDir(1:end-1); end
funcPreprocDir = [outputDir '/func'];
fmapPreprocDir = [outputDir '/fmap'];

% If funcTemplatePath isn't specified by user, define it. If it is, check if file exists.
if isempty(funcTemplatePath)
    if useTaskTemplate
        funcTemplatePath = [funcPreprocDir '/' fpp.bids.changeName(inputName,...
            {'desc','run','echo','space'},{'template','','1','task'},'sbref','.nii.gz')];
        if ~exist(funcTemplatePath,'file')
            funcTemplatePath = fpp.bids.changeName(funcTemplatePath,[],[],'bold');
        end
        funcTemplateSpace = 'task';
    else
        funcTemplatePath = [funcPreprocDir '/' fpp.bids.changeName(inputName,...
            {'desc','task','acq','run','echo','space'},{'template','','','','1','session'},'sbref','.nii.gz')];
        if ~exist(funcTemplatePath,'file')
            funcTemplatePath = fpp.bids.changeName(funcTemplatePath,[],[],'bold');
        end
        funcTemplateSpace = 'session';
    end
else
    [spaceInd,spaceEndInd] = regexp(funcTemplateSpace,'_space-[a-zA-Z0-9]+_');
    if sum(spaceInd>0)
        funcTemplateSpace = funcTemplateSpace(spaceInd(1)+7:spaceEndInd(1)-1);
    elseif isempty(funcTemplateSpace)
        error('funcTemplateSpace must be specified when funcTemplatePath lacks a BIDS space entity.')
    end
end

if ~exist(funcTemplatePath,'file')
    error(['Registration target ' funcTemplatePath ' does not exist. Run fpp.func.defineTemplate.']);
end
maskPath = fpp.bids.changeName(funcTemplatePath,{'desc','echo'},{'brain',[]},'mask','.nii.gz');
wmPath = fpp.bids.changeName(funcTemplatePath,{'desc','echo'},{'wmero1',[]},'mask','.nii.gz');
csfPath = fpp.bids.changeName(funcTemplatePath,{'desc','echo'},{'csfero1',[]},'mask','.nii.gz');
if ~exist(maskPath,'file')
    error(['Brain mask ' maskPath ' does not exist. Run fpp.func.register.']);
end

% Check if output exists.
finalOutputPath = fpp.bids.changeName([funcPreprocDir '/' inputName '.nii.gz'],...
    {'desc','echo','space'},{'preproc',[],funcTemplateSpace});
if exist(finalOutputPath,'file') && ~overwrite
    return;
end

% If overwriting, delete existing TEDANA output directory
tedanaDir = [funcPreprocDir '/' inputNameGeneric '_tedana'];
if exist(tedanaDir,'dir') && overwrite
    fpp.util.system(['rm -rf ' tedanaDir]);
end

% Generate output directories
if ~exist(funcPreprocDir,'dir'), mkdir(funcPreprocDir); end

% Copy raw functional data and metadata, convert to .nii.gz if necessary
for e=1:nEchoes
    [~,inputNames{e},~] = fpp.util.fileParts(inputPaths{e});
    outputPaths{e} = [funcPreprocDir '/' fpp.bids.changeName(inputNames{e},'desc','midprep0raw') '.nii.gz'];
    fpp.util.copyImageAndJson(inputPaths{e},outputPaths{e},'midprepfmri');
    fpp.bids.jsonChangeValue(outputPaths{e},{'Description','SkullStripped','RawSources','SpatialRef'},...
        {'Raw data copied to derivative directory.',false,fpp.bids.removeBidsDir(inputPaths{e}),'native'});
end
inputPathsRaw = outputPaths;
pathsToDelete = [pathsToDelete outputPaths];

% Define undistortion warp files, if undistorting
if undistort && (isempty(spinEchoPath) || isempty(topupWarpPath) ||  isempty(topupJacobianPath))
    [spinEchoPath,topupWarpPath,topupJacobianPath] = fpp.func.preproc.checkSpinEcho(inputPaths{1},fmapPreprocDir);
end

% Check TR
tr = fpp.util.checkMRIProperty('tr',outputPaths{1});
if isempty(tr)
    error('Could not compute TR of functional data.');
end
% Check # of volumes
vols = fpp.util.checkMRIProperty('vols',outputPaths{1});
if isempty(vols)
    error('Could not compute # of volumes from functional data header.');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 1: Estimate motion parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 1, Estimate motion parameters             - ' inputNameGeneric]);
steps{end+1} = 'motion parameter estimation (FSL''s mcflirt)';
inputPaths = outputPaths;
mcDir = [funcPreprocDir '/' inputNameGeneric '_motion'];
if ~exist(mcDir,'dir'), mkdir(mcDir); end
moCorrTargetVolNum = ceil(vols/2);   % Example func volume #, indexed by 0
motionParams = fpp.func.preproc.estimateHeadMotion(inputPaths{echoForMoCorr},mcDir,moCorrTargetVolNum);
artifactTPs = fpp.func.preproc.defineMotionArtifactTimePoints(motionParams,transCutoff,rotCutoff,tptsAfter);
% If a head movement occured at moCorr target, switch target volume number
if ismember(moCorrTargetVolNum,artifactTPs)
    goodTPs = setdiff(1:vols,artifactTPs);
    moCorrTargetVolNum = goodTPs(ceil(length(goodTPs)/2)+1)-1;
    motionParams = fpp.func.preproc.estimateHeadMotion(inputPaths{echoForMoCorr},mcDir,moCorrTargetVolNum);
    artifactTPs = fpp.func.preproc.defineMotionArtifactTimePoints(motionParams,transCutoff,rotCutoff,tptsAfter);
end
outlierPath = fpp.bids.changeName(inputPaths{1},{'desc','echo'},{[],[]},'outliers','.tsv');
outlierTSV = fpp.func.preproc.outlierTSV(artifactTPs,vols);
bids.util.tsvwrite(outlierPath,outlierTSV);
confoundPath = fpp.bids.changeName(inputPaths{1},{'desc','echo'},{[],[]},'confounds','.tsv');
pathsToDelete = [pathsToDelete mcDir];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 2: 3dDespike
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if useDespike
    fprintf('%s\n',['Step 2, Despike                                - ' inputNameGeneric]);
    steps{end+1} = 'despiking (AFNI 3dDespike)';
    for e=1:nEchoes
        outputPaths{e} = fpp.bids.changeName(outputPaths{e},'desc','midprep1despike');
%         if exist(outputPaths{e},'file'), continue; end      %%% TEMPORARY DEBUGGING HACK
        fpp.func.preproc.despike(inputPaths{e},outputPaths{e});
        fpp.bids.jsonChangeValue(outputPaths{e},{'Description','Sources'},...
            {fpp.func.preproc.description(midprepIntro,steps),fpp.bids.removeBidsDir(inputPaths{e})});
    end
end
pathsToDelete = [pathsToDelete outputPaths];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 3: Slice timing correct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if useSTC
    fprintf('%s\n',['Step 3, Slice timing correct                   - ' inputNameGeneric]);
    steps{end+1} = 'slice timing correction (FSL''s slicetimer)';
    inputPaths = outputPaths;
    [errorMsg,outputPaths] = fpp.func.preproc.sliceTimingCorrect(inputPaths,sliceTimes,tr);
    if ~isempty(errorMsg)
        error(errorMsg);
    end
    for e=1:nEchoes
        fpp.bids.jsonChangeValue(outputPaths{e},{'Description'},...
            {fpp.func.preproc.description(midprepIntro,steps)});
    end
end
pathsToDelete = [pathsToDelete outputPaths];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 4: Generate and unwarp MocoTarget vol
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 4, Generate/unwarp MocoTargetVol          - ' inputNameGeneric]);
mocoTargetPath = fpp.bids.changeName(outputPaths{echoForMoCorr},'desc','MocoTargetVol');
fpp.util.system(['fslroi ' outputPaths{echoForMoCorr} ' ' mocoTargetPath ' ' int2str(moCorrTargetVolNum) ' 1']);
fpp.bids.jsonReconstruct(outputPaths{echoForMoCorr},mocoTargetPath);
fpp.bids.jsonChangeValue(mocoTargetPath,'Description','Target volume for motion correction.');
% Undistorted MocoTarget vol
if undistort
    mocoTargetUndistortedPath = fpp.bids.changeName(outputPaths{echoForMoCorr},'desc','MocoTargetVolUndistorted');
    [xfmMocoTarget2SpinEcho,xfmSpinEcho2MocoTarget] = fpp.func.preproc.undistort(mocoTargetPath,mocoTargetUndistortedPath,...
        spinEchoPath,topupWarpPath,topupJacobianPath);
    fpp.bids.jsonChangeValue(mocoTargetUndistortedPath,'Description','Target volume for motion correction, undistorted.');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 5: Register to functional template
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 5, Register to FuncTemplate               - ' inputNameGeneric]);
if undistort
    mocoTargetToRegisterPath = mocoTargetUndistortedPath;
else
    mocoTargetToRegisterPath = mocoTargetPath;
end
if multiEcho && exist(fpp.bids.changeName(funcTemplatePath,'echo',echoForMoCorr),'file')
    funcTemplatePathRegTarget = fpp.bids.changeName(funcTemplatePath,'echo',echoForMoCorr);
else
    funcTemplatePathRegTarget = funcTemplatePath;
end
xfmMocoTarget2FuncTemplate = fpp.bids.changeName(mocoTargetPath,{'desc','from','to','mode','echo'},...
    {'','native','session','image',[]},'xfm','.mat');
xfmFuncTemplate2MocoTarget = fpp.bids.changeName(mocoTargetPath,{'desc','from','to','mode','echo'},...
    {'','session','native','image',[]},'xfm','.mat');
fpp.fsl.flirt(mocoTargetToRegisterPath,funcTemplatePathRegTarget,xfmMocoTarget2FuncTemplate,[],'cost','corratio',...
    'dof',6,'searchrx',[-180 180],'searchry',[-180 180],'searchrz',[-180 180]);
fpp.fsl.invertXfm(xfmMocoTarget2FuncTemplate,xfmFuncTemplate2MocoTarget);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 6: Apply motion/distortion correction, and registration to FuncTemplate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 6, Motion/distortion-correct and register - ' inputNameGeneric]);
steps{end+1} = ['one-shot motion/distortion correction and functional template '...
    'registration (using FSL''s topup, flirt, and applywarp)'];
inputPaths = outputPaths;
for e=1:nEchoes
    outputPaths{e} = fpp.bids.changeName(outputPaths{e},{'desc','space'},{'midprep3moco',funcTemplateSpace});
end
%if ~exist(outputPaths{end},'file')      % TEMPORARY DEBUGGING HACK
fpp.func.preproc.oneShotMotionDistortionCorrect(inputPaths,outputPaths,funcTemplatePath,...
    funcTemplateSpace,mcDir,topupWarpPath,topupJacobianPath,xfmMocoTarget2FuncTemplate,...
    xfmSpinEcho2MocoTarget,xfmMocoTarget2SpinEcho,echoForMoCorr);
%end
for e=1:nEchoes
    fpp.bids.jsonChangeValue(outputPaths{e},'Description',fpp.func.preproc.description(midprepIntro,steps));
end
pathsToDelete = [pathsToDelete outputPaths];
% Remove voxels with zero value at any time point after registration from brain mask
maskNonZeroPath = fpp.bids.changeName(outputPaths{1},{'echo','desc'},{[],'brainNonZero'},'mask');
fpp.fsl.maths(outputPaths{end},['-Tmin -bin -mul ' maskPath],maskNonZeroPath);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 7: TEDANA / Multi-echo Combine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if useTedana
    fprintf('%s\n',['Step 7, TEDANA Multi-echo ICA denoise          - ' inputNameGeneric]);
    stepsAlt = steps;
    steps(end+1:end+2) = {'brain masking (Freesurfer-based)','multi-echo ICA denoising (tedana)'};
    stepsAlt(end+1:end+2) = {'brain masking (Freesurfer-based)','optimal multi-echo combination (tedana''s t2smap)'};  % For optcomb output
    inputPaths = outputPaths;
    outputPath = fpp.bids.changeName(inputPaths{1},{'echo','desc'},{[],'midprep4tedana'});
    %if ~exist(outputPath,'file') %%% TEMPORARY DEBUGGING HACK
    fpp.func.preproc.tedana(inputPaths,outputPath,maskNonZeroPath,[],1,teVals);
    %end
elseif multiEcho
    fprintf('%s\n',['Step 7, Multi-echo combine                     - ' inputNameGeneric]);
    steps{end+1} = 'optimal multi-echo combination (tedana''s t2smap)';
    inputPaths = outputPaths;
    outputPath = fpp.bids.changeName(inputPaths{1},{'echo','desc'},{[],'midprep4optcomb'});
    fpp.func.preproc.tedana(inputPaths,outputPath,maskNonZeroPath,[],0,teVals);
else
    inputPath = outputPaths{1};
    outputPath = fpp.bids.changeName(inputPath,'desc','midprep4mask');
    steps{end+1} = 'brain masking (Freesurfer-based)';
    fpp.fsl.maths(inputPath,['-mul ' maskNonZeroPath],outputPath);
end
pathsToDelete = [pathsToDelete outputPath];
fpp.bids.jsonChangeValue(outputPath,{'Description','SkullStripped'},...
    {fpp.func.preproc.description(midprepIntro,steps),true});
if useTedana
    fpp.bids.jsonChangeValue(fpp.bids.changeName(outputPath,'desc','midprep4optcomb'),{'Description','SkullStripped'},...
        {fpp.func.preproc.description(midprepIntro,stepsAlt),true});
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 7.5: Temporal filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if tempFilt
    fprintf('%s\n',['Step 7.5, Temporal filter                      - ' inputNameGeneric]);
    if useTedana
        warning('Temporal filtering is not recommended in additional to ME-ICA denoising.');
    end
    switch filtType
        case 'high'
            filtTypeStr = 'highpass';
            filtCutoffStr = [num2str(filtCutoff) 'Hz'];
        case 'low'
            filtTypeStr = 'lowpass';
            filtCutoffStr = [num2str(filtCutoff) 'Hz'];
        case 'bandpass'
            filtTypeStr = 'bandpass';
            filtCutoffStr = [num2str(filtCutoff(1)) '-' num2str(filtCutoff(1)) 'Hz'];
    end
    steps{end+1} = [filtTypeStr ' temporal filtering (' filtCutoffStr ', using MATLAB''s fir1 and filtfilt)'];
    inputPath = outputPath;
    outputPath = fpp.bids.changeName(inputPath,'desc',{[],'midprep5tempfilt'});
    fpp.util.mriFilter(inputPath,outputPath,filtCutoff,filtType,filtOrder,tr);
    fpp.bids.jsonReconstruct(inputPath,outputPath,'midprepfmri');
    fpp.bids.jsonChangeValue(outputPath,{'Description','Sources'},...
        {fpp.func.preproc.description(midprepIntro,steps),fpp.bids.removeBidsDir(inputPath)});
    pathsToDelete = [pathsToDelete outputPath];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 8: Intensity normalize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 8, Intensity normalize                    - ' inputNameGeneric]);
steps{end+1} = 'intensity normalization';
inputPath = outputPath;
outputPath = fpp.bids.changeName(inputPath,'desc','midprep6intnorm');
% Normalize image median (also called "grand mean scaling" in SPM)
[~,funcMedian] = fpp.util.system(['fslstats ' inputPath ' -k ' maskNonZeroPath ' -p 50']);
funcMedian = str2num(strtrim(funcMedian));
fpp.fsl.maths(inputPath,['-mul ' num2str(newMedian/funcMedian)],outputPath);
fpp.bids.jsonReconstruct(inputPath,outputPath,'midprepfmri');
fpp.bids.jsonChangeValue(outputPath,{'Description','Sources'},...
    {fpp.func.preproc.description(midprepIntro,steps),fpp.bids.removeBidsDir(inputPath)});
pathsToDelete = [pathsToDelete outputPath];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 9: Extract nuisance signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 9, Extract nuisance signals               - ' inputNameGeneric]);
confoundTSV = bids.util.tsvread(confoundPath);
dataMat = fpp.util.readDataMatrix(outputPath);
maskMat = fpp.util.readDataMatrix(maskNonZeroPath);
wmMat = fpp.util.readDataMatrix(wmPath);
csfMat = fpp.util.readDataMatrix(csfPath);
confoundTSV.global_signal = mean(dataMat.*maskMat)';
confoundTSV.white_matter = mean(dataMat.*wmMat)';
confoundTSV.csf = mean(dataMat.*csfMat)';
[confoundTSV.dvars,confoundTSV.dvars_std] = fpp.func.preproc.dvars(inputPath,maskNonZeroPath);
bids.util.tsvwrite(confoundPath,confoundTSV);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 10: Resample to cortical surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 10, Resample to cortical surface          - ' inputNameGeneric]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 11: Spatially smooth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if fwhm>0
    fprintf('%s\n',['Step 11, Spatially smooth                      - ' inputNameGeneric]);
    steps{end+1} = 'spatial smoothing';
    inputPath = outputPath;
    outputPath = fpp.bids.changeName(inputPath,'desc',{[],'midprep7smooth'});
    
    % Spatial smoothing here!
    % FSL's method relies on brain-masked data.
    
    
    
    
    
    
    fpp.bids.jsonReconstruct(inputPath,outputPath,'midprepfmri');
    fpp.bids.jsonChangeValue(outputPath,{'Description','Sources'},...
        {fpp.func.preproc.description(midprepIntro,steps),fpp.bids.removeBidsDir(inputPath)});
    pathsToDelete = [pathsToDelete outputPath];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 11.5: Rename output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputPath = outputPath;
outputPath = fpp.bids.changeName(inputPath,'desc','preproc');
fpp.util.system(['mv ' inputPath ' ' outputPath]);
if exist(fpp.bids.jsonPath(inputPath))
    fpp.util.system(['mv ' fpp.bids.jsonPath(inputPath) ' ' fpp.bids.jsonPath(outputPath)]);
end
fpp.bids.jsonReconstruct(outputPath,outputPath,'fmri');
fpp.bids.jsonChangeValue(outputPath,'Description',fpp.func.preproc.description(preprocIntro,steps));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 11.6: Generate carpet plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuisanceSeries = [confoundTSV.global_signal confoundTSV.white_matter confoundTSV.csf...
    confoundTSV.framewise_displacement confoundTSV.dvars_std];
nuisanceNames = {'Global mean','WM mean','CSF mean','Framewise displacement (mm)','DVARS_std'};
segmentMaskPaths = {fpp.bids.changeName(maskPath,'desc','gm'),fpp.bids.changeName(maskPath,'desc','wm'),...
    fpp.bids.changeName(maskPath,'desc','csf')};
carpetPlotPath = fpp.bids.changeName(outputPath,[],[],'carpetplot','.png');
segmentColors = {[0 .5 1],[0 1 0],[1 1 0]};
nuisanceColors = {[61 165 193]/255,[45 167 111]/255,[164 146 43]/255,...
        [246 102 126]/255,[197 111 242]/255};
fpp.util.carpetPlot(outputPath,maskNonZeroPath,segmentMaskPaths,nuisanceSeries,...
    nuisanceNames,carpetPlotPath,segmentColors,nuisanceColors);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 12: Compute tSNR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 12, Compute tSNR                          - ' inputNameGeneric]);
% Multi-echo combine raw input data for tSNR comparison
if multiEcho
    inputPathRaw = fpp.bids.changeName(inputPathsRaw{1},{'echo','desc','space'},{[],'rawoptcomb','native'});
    outputDescription = 'Raw data, optimally combined across echoes using tedana.';
    fpp.func.preproc.tedana(inputPathsRaw,inputPathRaw,maskNonZeroPath,outputDescription,0);
else
    inputPathRaw = inputPathsRaw{1};
end
% Undistort raw data
fpp.func.preproc.undistort(inputPathRaw,inputPathRaw,spinEchoPath,topupWarpPath,...
    topupJacobianPath,xfmMocoTarget2SpinEcho,xfmSpinEcho2MocoTarget);
% Register raw data to funcTemplate
inputPathRaw2FuncTemplate = fpp.bids.changeName(inputPathRaw,'space',funcTemplateSpace);
fpp.fsl.moveImage(inputPathRaw,funcTemplatePath,inputPathRaw2FuncTemplate,xfmMocoTarget2FuncTemplate);
fpp.bids.jsonReconstruct(inputPathRaw,inputPathRaw2FuncTemplate,'fmri');
fpp.bids.jsonChangeValue(inputPathRaw2FuncTemplate,'Description',...
    'Raw data, optimally combined across echoes, undistorted, and registered to functional template.');
pathsToDelete = [pathsToDelete inputPathRaw];
% Compute tSNR maps
if genTSNR
    fpp.util.tsnrMap(inputPathRaw2FuncTemplate);
    fpp.util.tsnrMap(outputPath);
    fpp.util.tsnrMap(fpp.bids.changeName(outputPath,'desc','midprep4optcomb'));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CLEANUP: Delete midprep files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if deleteMidprep
    for i=1:length(pathsToDelete)
        fpp.util.system(['rm -rf ' pathsToDelete{i}]);
        if exist(fpp.bids.jsonPath(pathsToDelete{i}),'file')
            fpp.util.system(['rm -rf ' fpp.bids.jsonPath(pathsToDelete{i})]);
        end
    end
end


fprintf('%s\n\n',['Finished preproc for ' inputName]);

end