%
% fpp.func.preproc(inputPaths,outputDir,varargin)
% 
% Preprocesses a single (multi-echo) fMRI dataset, including including
% motion parameter estimation, despiking, slice timing correction, one-shot
% motion and distoration correction and template registration, TEDANA 
% multi-echo ICA denoising and optimal echo combination, intensity
% normalization, and optional temporal and spatial filtering. Should be run
% after anat/fmap scripts, func.defineTemplate, and func.register.
% 
% Example usage:
% fpp.func.preproc({'/pathToData/sub-01_task-faceloc_run-01_echo-1_bold.nii.gz',...
%       '/pathToData/sub-01_task-faceloc_run-01_echo-2_bold.nii.gz'},'/bidsDir/sub-01/');
% 
% Arguments:
% - inputPaths (string or cell array of strings): paths to input fMRI data,
%       one path for each echo.
% - outputDir (string): path to output directory (subject/session dir, in
%       BIDS filesystem)
% 
% Variable arguments:
% - overwrite (boolean; default=0): whether to overwrite files that have
%       already been written by this function.
% - fwhm (scalar in [0,Inf); default=0): FWHM of spatial smoothing kernel 
%       (mm).
% - disdaqs (integer in [0,Inf); default=0): number of disdaq volumes (at 
%       beginning of run) to remove. [NOT CURRENTLY IMPLEMENTED]
% - fdCutoff (scalar in [0,Inf]; default=.5): artifact detection cutoff for
%       framewise displacement (mm)
% - transCutoff (scalar in [0,Inf]; default=Inf): artifact detection cutoff
%       for framewise total translation (mm)
% - rotCutoff (scalar in [0,Inf]; default=Inf): artifact detection cutoff
%   for framewise total rotation (degrees)
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
% - templateType (string, 'anat' or 'func'): whether to use functional or
%       low-res anatomical template
% - templateRes (string, default='2'): anatomical template resolution
% - anatTemplatePath (string): path to low-res anatomical image to register
%       to (default: determined from input name)
% - funcTemplatePath (string): path to functional template image to
%       register to (default: determined from input name)
% - funcTemplateSpace (string): space of functional template, needed if
%       using non-standard funcTemplatePath that lacks a BIDS space entity
% - useTaskTemplate (boolean, default=0): whether to use a task-specific
%       template
% - undistort (boolean, default=1): Whether to distortion-correct
%       functional data using blip-up/blip-down method
% - spinEchoPath (string): path to spin echo image matched in phase-encode
%       direction to input fMRI data
% - topupWarpPath (string): path to undistortion warp image produced by
%       topup
% - topupJacobianPath (string): path to undistortion warp jacobian image
%       produced by topup
% - outputOptcomb (boolean): whether to output data preprocessed without 
%       Tedana (ME-ICA denoising), just optimal echo combination
% - deleteMidprep (boolean; default=1): whether to delete midprep images
%

function preproc(inputPaths,outputDir,varargin)

% TODO NEXT:
% - Add JSON files for additional TEDANA/ts2map outputs
% - If deleteMidprep==1, remove "Sources" field of final outputs?
% - Add suffix option for desc. Need to edit input/outputNameGeneric, and
%   edit desc of output TEDANA folder (defitions both within preproc
%   script, and tedana script - new tedanaDir optional argument for tedana
%   script?). Also need to edit any other name that has desc removed -
%   e.g. outlier TSV file preproc step 1); xfm matrices
%   (xfmMocoTarget2FuncTemplate and inverse in step 5)
%
% TODO EVENTUALLY:
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
useDespike = 1;                 % Whether to apply AFNI's 3dDepike before STC
useSTC = 1;                     % Whether to use slice timing correction
sliceTimes = [];                % Slice time (s) relative to onset of volume acquisition (Default: read from JSON file)
useTedana = 1;                  % Whether to use tedana-based denoising for multi-echo data
newMedian = 10000;              % Target value for intensity normalization (median of data in brain mask is set to this)
fwhm = 0;                       % FWHM of spatial smoothing kernel (mm)
nComps = 5;                     % Number of PCs of WM and CSF to include in confound file
outputOptcomb = 1;              % Whether to output data preprocessed without Tedana (ME-ICA denoising)
deleteMidprep = 1;              % Whether to delete midprep images

% Template parameters
templateType = 'anat';          % Whether to register to functional or anatomical template
templateRes = '2';              % Resolution of anatomical template
anatTemplatePath = '';
funcTemplatePath = '';          % Functional template image to register to
funcTemplateSpace = '';         % Space of func template, required if using custom funcTemplatePath with BIDs space entity
useTaskTemplate = 0;            % Use a task/acq-specific functional template, rather than a session template

% Undistortion parameters
undistort = 1;                  % Whether to distortion-correct functional data using blip-up/blip-down method
spinEchoPath = '';              % Spin echo path with PE dir matched to functional, if not using default (determined from json metadata)
topupWarpPath = '';             % Undistortion warp path, if not using default (determined from json metadata)
topupJacobianPath = '';         % Undistortion warp jacobian path, if not using default (determined from json metadata)

% Artifact detection parameters
disdaqs = 0;                    % Number of disdaq volumes (at beginning of run) to remove. [NOT CURRENTLY IMPLEMENTED]
fdCutoff = .5;                  % Artifact detection cutoff: framewise displacement (mm)
transCutoff = Inf;              % Artifact detection cutoff: framewise total translation (mm)
rotCutoff = Inf;                % Artifact detection cutoff: framewise total rotation (degrees)
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
tsnrIntro = 'tSNR of preprocessed data generated by fmriPermPipe.';         % Description intro sentence for tSNR maps
steps = {};                     % Processing step descriptions

% Edit variable arguments.  Note: optInputs checks for proper input.
varArgList = {'overwrite','funcTemplatePath','funcTemplateSpace','fwhm',...
    'tempFilt','teVals','transCutoff','sliceTimes','useTaskTemplate',...
    'rotCutoff','tptsAfter','stdCutoff','disdaqs','genTSNR','plotResults',...
    'faValue','useSTC','useTedana','echoForMoCorr','undistort','fdCutoff'...
    'useDespike','topupWarpPath','topupJacobianPath','spinEchoPath',...
    'deleteMidprep','templateType','templateRes','anatTemplatePath',...
    'outputOptcomb','nComps'};
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
        error(['Input path ' inputPaths{e} ' does not exist.']);
    end
end

% Define input name, output directories
[~,inputName,~] = fpp.util.fileParts(inputPaths{1});
if strcmp(outputDir(end),'/'), outputDir = outputDir(1:end-1); end
funcPreprocDir = [outputDir '/func'];
fmapPreprocDir = [outputDir '/fmap'];
anatPreprocDir = [outputDir '/anat'];

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
    funcTemplateSpace = fpp.bids.checkNameValue(funcTemplatePath);
    if isempty(funcTemplateSpace)
        error('funcTemplateSpace must be specified when funcTemplatePath lacks a BIDS space entity.')
    end
end
[~,funcTemplateName,~] = fpp.util.fileParts(funcTemplatePath);

if ~exist(funcTemplatePath,'file')
    error(['Functional template image ' funcTemplatePath ' does not exist. Run fpp.func.defineTemplate.']);
end

% Process template space options
if strcmpi(templateType,'anat')
    templateSpace = 'individual';
    if ~isempty(anatTemplatePath)
        templatePath = anatTemplatePath;
    else
        templatePath = [anatPreprocDir '/' fpp.bids.changeName(inputName,...
            {'desc','task','acq','run','echo','space','res'},{'preproc','','','','',...
            templateSpace,templateRes},'T1w','.nii.gz')];
    end
    if ~exist(templatePath,'file')
        error(['Template image ' templatePath ' does not exist.'...
            ' Run fpp.anat.postproc or specify a different image.']);
    end
elseif strcmpi(templateType,'func')
    templatePath = funcTemplatePath;
    templateSpace = funcTemplateSpace;
else
    error('Template type must be specified as func or anat.');
end
[~,templateName,~] = fpp.util.fileParts(templatePath);

% Default temporal filter type
if isempty(filtType)
    if length(filtCutoff)==1
        filtType = 'high';
    else
        filtType = 'bandpass';
    end
end

% Define brain mask and segment paths
maskPath = [anatPreprocDir '/' fpp.bids.changeName(templateName,{'desc','echo'},{'brain',[]},'mask','.nii.gz')];
gmPath = [anatPreprocDir '/' fpp.bids.changeName(templateName,{'desc','echo'},{'gm',[]},'mask','.nii.gz')];
wmPath = [anatPreprocDir '/' fpp.bids.changeName(templateName,{'desc','echo'},{'wm',[]},'mask','.nii.gz')];
csfPath = [anatPreprocDir '/' fpp.bids.changeName(templateName,{'desc','echo'},{'csf',[]},'mask','.nii.gz')];
wmEroPath = [anatPreprocDir '/' fpp.bids.changeName(templateName,{'desc','echo'},{'wmero1',[]},'mask','.nii.gz')];
csfEroPath = [anatPreprocDir '/' fpp.bids.changeName(templateName,{'desc','echo'},{'csfero1',[]},'mask','.nii.gz')];
if ~exist(maskPath,'file')
    error(['Brain mask ' maskPath ' does not exist. Run fpp.func.register.']);
end

% Define generic labels for input/output imagesName
inputNameGeneric = strrep(fpp.bids.changeName(inputName,'echo',[]),'_bold','');
outputNameGeneric = strrep(fpp.bids.changeName(inputName,{'echo','space'},{[],templateSpace}),'_bold','');

% Check if output exists.
finalOutputPath = fpp.bids.changeName([funcPreprocDir '/' inputName '.nii.gz'],...
    {'desc','echo','space'},{'preproc',[],templateSpace});
if exist(finalOutputPath,'file') && ~overwrite
    return;
end

% If overwriting, delete existing TEDANA output directory
tedanaDir = [funcPreprocDir '/' outputNameGeneric '_tedana'];   % Depends on output space
if exist(tedanaDir,'dir') && overwrite
    fpp.util.system(['rm -rf ' tedanaDir]);
end

% Generate output directory
if ~exist(funcPreprocDir,'dir'), mkdir(funcPreprocDir); end

% Define undistortion warp files, if undistorting
if undistort && (isempty(spinEchoPath) || isempty(topupWarpPath) ||  isempty(topupJacobianPath))
    [spinEchoPath,topupWarpPath,topupJacobianPath] = fpp.func.preproc.checkSpinEcho(inputPaths{1},fmapPreprocDir);
end

% Check TR
tr = fpp.util.checkMRIProperty('tr',inputPaths{1});
if isempty(tr)
    error('Could not compute TR of functional data.');
end
% Check # of volumes
vols = fpp.util.checkMRIProperty('vols',inputPaths{1});
if isempty(vols)
    error('Could not compute # of volumes from functional data header.');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 0.5: Copy raw functional data and metadata
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for e=1:nEchoes
    [~,inputNames{e},~] = fpp.util.fileParts(inputPaths{e});
    outputPaths{e} = [funcPreprocDir '/' fpp.bids.changeName(inputNames{e},'desc','midprep0raw') '.nii.gz'];
    fpp.util.copyImageAndJson(inputPaths{e},outputPaths{e},'midprepfmri');
    fpp.bids.jsonChangeValue(outputPaths{e},{'Description','SkullStripped','RawSources','SpatialRef'},...
        {'Raw data copied to derivative directory.',false,fpp.bids.removeBidsDir(inputPaths{e}),'native'});
end
% inputPathsRaw = outputPaths;
pathsToDelete = [pathsToDelete outputPaths];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 1: Estimate motion parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 1, Estimate motion parameters             - ' outputNameGeneric]);
steps{end+1} = 'motion parameter estimation (FSL''s mcflirt)';
inputPaths = outputPaths;
mcDir = [funcPreprocDir '/' inputNameGeneric '_motion'];
if ~exist(mcDir,'dir'), mkdir(mcDir); end
moCorrTargetVolNum = ceil(vols/2);   % Example func volume #, indexed by 0
confoundTSV = fpp.func.preproc.estimateHeadMotion(inputPaths{echoForMoCorr},mcDir,moCorrTargetVolNum);
artifactTPs = fpp.func.preproc.defineMotionArtifactTimePoints(confoundTSV,fdCutoff,transCutoff,rotCutoff,tptsAfter);
% If a head movement occured at moCorr target, switch target volume number
if ismember(moCorrTargetVolNum,artifactTPs)
    goodTPs = setdiff(1:vols,artifactTPs);
    moCorrTargetVolNum = goodTPs(ceil(length(goodTPs)/2)+1)-1;
    confoundTSV = fpp.func.preproc.estimateHeadMotion(inputPaths{echoForMoCorr},mcDir,moCorrTargetVolNum);
    artifactTPs = fpp.func.preproc.defineMotionArtifactTimePoints(confoundTSV,fdCutoff,transCutoff,rotCutoff,tptsAfter);
end
outlierPath = fpp.bids.changeName(inputPaths{1},{'desc','echo'},{[],[]},'outliers','.tsv');
outlierTSV = fpp.func.preproc.outlierTSV(artifactTPs,vols);
bids.util.tsvwrite(outlierPath,outlierTSV);
pathsToDelete = [pathsToDelete mcDir];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 2: 3dDespike
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if useDespike
    fprintf('%s\n',['Step 2, Despike                                - ' outputNameGeneric]);
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
    fprintf('%s\n',['Step 3, Slice timing correct                   - ' outputNameGeneric]);
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
fprintf('%s\n',['Step 4, Generate/unwarp MocoTargetVol          - ' outputNameGeneric]);
mocoTargetPath = fpp.bids.changeName(outputPaths{echoForMoCorr},'desc','MocoTargetVol');
if ~exist(mocoTargetPath,'file') || overwrite
    fpp.util.system(['fslroi ' outputPaths{echoForMoCorr} ' ' mocoTargetPath ' ' int2str(moCorrTargetVolNum) ' 1']);
    fpp.bids.jsonReconstruct(outputPaths{echoForMoCorr},mocoTargetPath);
    fpp.bids.jsonChangeValue(mocoTargetPath,'Description','Target volume for motion correction.');
else
    warning(['Not overwriting existing moco target image. If motion artifact thresholds '...
        'have changed, set overwrite to 1 to rewrite this image, or add a suffix.']);
end
% Undistort MocoTarget vol
if undistort
    mocoTargetUndistortedPath = fpp.bids.changeName(outputPaths{echoForMoCorr},'desc','MocoTargetVolUndistorted');
    if ~exist(mocoTargetUndistortedPath,'file') || overwrite
        [xfmMocoTarget2SpinEcho,xfmSpinEcho2MocoTarget] = fpp.func.preproc.undistort(mocoTargetPath,mocoTargetUndistortedPath,...
            spinEchoPath,topupWarpPath,topupJacobianPath);
        fpp.bids.jsonChangeValue(mocoTargetUndistortedPath,'Description','Target volume for motion correction, undistorted.');
    else
        xfmMocoTarget2SpinEcho = fpp.bids.changeName(mocoTargetPath,{'desc','from','to','mode','echo'},...
        {'','native','SpinEcho','image',[]},'xfm','.mat');
        xfmSpinEcho2MocoTarget = fpp.bids.changeName(mocoTargetPath,{'desc','from','to','mode','echo'},...
        {'','SpinEcho','native','image',[]},'xfm','.mat');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 5: Register to functional template
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 5, Register to template                   - ' outputNameGeneric]);
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
if ~exist(xfmMocoTarget2FuncTemplate,'file') || ~exist(xfmFuncTemplate2MocoTarget,'file') || overwrite
    fpp.fsl.flirt(mocoTargetToRegisterPath,funcTemplatePathRegTarget,xfmMocoTarget2FuncTemplate,[],'cost','corratio',...
        'dof',6,'searchrx',[-180 180],'searchry',[-180 180],'searchrz',[-180 180]);
    fpp.fsl.invertXfm(xfmMocoTarget2FuncTemplate,xfmFuncTemplate2MocoTarget);
end
if strcmpi(templateType,'anat')
    xfmFunc2AnatTemplate = [anatPreprocDir '/' fpp.bids.changeName(funcTemplateName,...
        {'desc','space','res','echo','from','to','mode'},...
        {[],[],[],[],funcTemplateSpace,'individual','image'},'xfm','.mat')];
    xfmMocoTarget2Template = fpp.bids.changeName(inputPaths{1},...
        {'desc','space','res','echo','from','to','mode'},...
        {[],[],[],[],'native','individual','image'},'xfm','.mat');
    if ~exist(xfmMocoTarget2Template,'file') || overwrite
        fpp.fsl.concatXfm(xfmFunc2AnatTemplate,xfmMocoTarget2FuncTemplate,xfmMocoTarget2Template);
    end
else
    xfmMocoTarget2Template = xfmMocoTarget2FuncTemplate;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 6: Apply motion/distortion correction, and template registration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 6, Motion/distortion-correct and register - ' outputNameGeneric]);
steps{end+1} = ['one-shot motion/distortion correction and template '...
    'registration (using FSL''s topup, flirt, and applywarp)'];
inputPaths = outputPaths;
for e=1:nEchoes
    outputPaths{e} = fpp.bids.changeName(outputPaths{e},{'desc','space'},{'midprep3moco',templateSpace});
end
if ~exist(outputPaths{end},'file')      % TEMPORARY DEBUGGING HACK
fpp.func.preproc.oneShotMotionDistortionCorrect(inputPaths,outputPaths,templatePath,...
    templateSpace,mcDir,topupWarpPath,topupJacobianPath,xfmMocoTarget2Template,...
    xfmSpinEcho2MocoTarget,xfmMocoTarget2SpinEcho,echoForMoCorr);
end
for e=1:nEchoes
    fpp.bids.jsonChangeValue(outputPaths{e},'Description',fpp.func.preproc.description(midprepIntro,steps));
end
pathsToDelete = [pathsToDelete outputPaths];
% Remove voxels with zero value at any time point after registration from brain mask
maskNonZeroPath = fpp.bids.changeName(outputPaths{1},{'echo','desc'},{[],'brainNonZero'},'mask');
fpp.fsl.maths(outputPaths{end},['-Tmin -bin -mul ' maskPath],maskNonZeroPath);
fpp.bids.jsonReconstruct(outputPaths{end},maskNonZeroPath,'mri');
fpp.bids.jsonChangeValue(maskNonZeroPath,{'Sources','Type','Description','SkullStripped'},...
    {{fpp.bids.removeBidsDir(outputPaths{1}),fpp.bids.removeBidsDir(maskPath)},...
    'Brain','Brain mask intersected with mask of nonzero voxels from this task/run.',[]});



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 7: TEDANA / Multi-echo Combine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: after this point, inputPaths, outputPaths, and steps are cell
% arrays indexing to 1) tedana and 2) optcomb outputs (not echo #!)

tmpsteps = steps; clear steps; steps{1} = tmpsteps;
descSuffices = {'','NoTedana'};
if useTedana
    fprintf('%s\n',['Step 7, TEDANA Multi-echo ICA denoise          - ' outputNameGeneric]);
    steps{2} = steps{1};
    steps{1}(end+1:end+2) = {'brain masking (Freesurfer-based)','multi-echo ICA denoising (tedana)'};
    steps{2}(end+1:end+2) = {'brain masking (Freesurfer-based)','optimal multi-echo combination (tedana''s t2smap)'};  % For optcomb output
    inputPaths = outputPaths;
    outputPathsNew{1} = fpp.bids.changeName(inputPaths{1},{'echo','desc'},{[],'midprep4tedana'});
    outputPathsNew{2} = fpp.bids.changeName(outputPathsNew{1},'desc','midprep4optcomb');
    outputPaths = outputPathsNew;
%     if ~exist(outputPaths{1},'file') %%% TEMPORARY DEBUGGING HACK
    fpp.func.preproc.tedana(inputPaths,outputPaths{1},maskNonZeroPath,[],1,teVals);
%     end
elseif multiEcho
    fprintf('%s\n',['Step 7, Multi-echo combine                     - ' outputNameGeneric]);
    steps{1}(end+1:end+2) = {'brain masking (Freesurfer-based)','optimal multi-echo combination (tedana''s t2smap)'};
    inputPaths = outputPaths;
    outputPathsNew{1} = fpp.bids.changeName(inputPaths{1},{'echo','desc'},{[],'midprep4optcomb'});
    outputPaths = outputPathsNew;
    fpp.func.preproc.tedana(inputPaths,outputPaths{1},maskNonZeroPath,[],0,teVals);
else
    fprintf('%s\n',['Step 7, Brain mask                             - ' outputNameGeneric]);
    steps{1}{end+1} = 'brain masking (Freesurfer-based)';
    inputPaths{1} = outputPaths{1};
    outputPathsNew{1} = fpp.bids.changeName(inputPaths{1},'desc','midprep4mask');
    outputPaths = outputPathsNew;
    fpp.fsl.maths(inputPaths{1},['-mul ' maskNonZeroPath],outputPaths{1});
end
outputPaths = outputPathsNew;   clear inputPaths;
for i=1:length(outputPaths)
    fpp.bids.jsonChangeValue(outputPaths{i},{'Description','SkullStripped'},...
        {fpp.func.preproc.description(midprepIntro,steps{i}),true});
    pathsToDelete = [pathsToDelete outputPaths{i}];
end
if useTedana && ~outputOptcomb
    outputPaths = outputPaths(1);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 7.5: Temporal filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if tempFilt
    fprintf('%s\n',['Step 7.5, Temporal filter                      - ' outputNameGeneric]);
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
    for i=1:length(outputPaths)
        steps{i}{end+1} = [filtTypeStr ' temporal filtering (' filtCutoffStr ', using MATLAB''s fir1 and filtfilt)'];
        inputPaths{i} = outputPaths{i};
        outputPaths{i} = fpp.bids.changeName(inputPaths{i},'desc',['midprep5tempfilt' descSuffices{i}]);
        fpp.util.mriFilter(inputPaths{i},outputPaths{i},filtCutoff,filtType,filtOrder,tr);
        fpp.bids.jsonReconstruct(inputPaths{i},outputPaths{i},'midprepfmri');
        fpp.bids.jsonChangeValue(outputPaths{i},{'Description','Sources'},...
            {fpp.func.preproc.description(midprepIntro,steps{i}),fpp.bids.removeBidsDir(inputPaths{i})});
        pathsToDelete = [pathsToDelete outputPaths{i}];
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 8: Intensity normalize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalize image median (also called "grand mean scaling" in SPM)
fprintf('%s\n',['Step 8, Intensity normalize                    - ' outputNameGeneric]);
for i=1:length(outputPaths)
    steps{i}{end+1} = 'intensity normalization';
    inputPaths{i} = outputPaths{i};
    outputPaths{i} = fpp.bids.changeName(inputPaths{i},'desc',['midprep6intnorm' descSuffices{i}]);
    [~,funcMedian] = fpp.util.system(['fslstats ' inputPaths{i} ' -k ' maskNonZeroPath ' -p 50']);
    funcMedian = str2num(strtrim(funcMedian));
    fpp.fsl.maths(inputPaths{i},['-mul ' num2str(newMedian/funcMedian)],outputPaths{i});
    fpp.bids.jsonReconstruct(inputPaths{i},outputPaths{i},'midprepfmri');
    fpp.bids.jsonChangeValue(outputPaths{i},{'Description','Sources'},...
        {fpp.func.preproc.description(midprepIntro,steps{i}),fpp.bids.removeBidsDir(inputPaths{i})});
    pathsToDelete = [pathsToDelete outputPaths{i}];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 9: Extract nuisance signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 9, Extract nuisance signals               - ' outputNameGeneric]);
for i=1:length(outputPaths)
    confoundTSVs{i} = confoundTSV;
    confoundPath = fpp.bids.changeName(outputPaths{i},'desc',descSuffices{i},'confounds','.tsv');
    mask = fpp.util.mriRead(maskNonZeroPath);
    dataMat = fpp.util.readDataMatrix(outputPaths{i},mask.vol);
    maskMat = fpp.util.readDataMatrix(maskNonZeroPath);
    wmMat = fpp.util.readDataMatrix(wmEroPath,mask.vol);
    csfMat = fpp.util.readDataMatrix(csfEroPath,mask.vol);
    confoundTSVs{i}.global_signal = mean(dataMat)';
    confoundTSVs{i}.white_matter = mean(dataMat.*wmMat)';
    confoundTSVs{i}.csf = mean(dataMat.*csfMat)';
    [confoundTSVs{i}.dvars,confoundTSVs{i}.dvars_std] = fpp.func.preproc.dvars(inputPaths{i},maskNonZeroPath);
    % Anatomical CompCorr regressors
    wmData = dataMat(wmMat==1,:);
    wmData = wmData - mean(wmData,2);
    [~,~,V] = svd(wmData,'econ');
    for c=1:nComps
        confoundTSVs{i}.(['acompcorr_wm_pc' int2str(c)]) = V(:,c);
    end
    csfData = dataMat(csfMat==1,:);
    if size(csfData,1)<nComps
        warning(['Eroded CSF ROI is less than ' int2str(nComps) ' - not including CSF aCompCorr regressors.']);
    else
        csfData = csfData - mean(csfData,2);
        [~,~,V] = svd(csfData,'econ');
        for c=1:nComps
            confoundTSVs{i}.(['acompcorr_csf_pc' int2str(c)]) = V(:,c);
        end
    end
    bids.util.tsvwrite(confoundPath,confoundTSVs{i});
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 10: Spatially smooth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if fwhm>0
    fprintf('%s\n',['Step 10, Spatially smooth                      - ' outputNameGeneric]);
    for i=1:length(outputPaths)
        steps{i}{end+1} = ['volumetric spatial smoothing within tissue component (GM, WM, CSF), '...
            num2str(fwhm) 'mm-fwhm Gaussian kernel'];
        inputPaths{i} = outputPaths{i};
        outputPaths{i} = fpp.bids.changeName(inputPaths{i},'desc',{[],['midprep7smooth' descSuffices{i}]});
        segmentPaths = {gmPath,wmPath,csfPath};
        fpp.util.copyImageAndJson(inputPaths{i},outputPaths{i},'midprepfmri');
        for s=1:length(segmentPaths)
            fpp.util.smoothInMask(outputPaths{i},segmentPaths{s},fwhm,outputPaths{i});
        end
        fpp.bids.jsonChangeValue(outputPaths{i},{'Description','Sources'},...
            {fpp.func.preproc.description(midprepIntro,steps{i}),fpp.bids.removeBidsDir(inputPaths{i})});
        pathsToDelete = [pathsToDelete outputPaths{i}];
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 10.5: Rename output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(outputPaths)
    inputPaths{i} = outputPaths{i};
    outputPaths{i} = fpp.bids.changeName(inputPaths{i},'desc',['preproc' descSuffices{i}]);
    fpp.util.system(['mv ' inputPaths{i} ' ' outputPaths{i}]);
    if exist(fpp.bids.jsonPath(inputPaths{i}),'file')
        fpp.util.system(['mv ' fpp.bids.jsonPath(inputPaths{i}) ' ' fpp.bids.jsonPath(outputPaths{i})]);
    end
    fpp.bids.jsonReconstruct(outputPaths{i},outputPaths{i},'fmri');
    fpp.bids.jsonChangeValue(outputPaths{i},'Description',fpp.func.preproc.description(preprocIntro,steps{i}));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 10.6: Generate carpet plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuisanceNames = {'Global mean','WM mean','CSF mean','Framewise displacement (mm)','DVARS_std'};
segmentMaskPaths = {fpp.bids.changeName(maskPath,'desc','gm'),fpp.bids.changeName(maskPath,'desc','wm'),...
    fpp.bids.changeName(maskPath,'desc','csf')};
segmentColors = {[0 .5 1],[0 1 0],[1 1 0]};
nuisanceColors = {[61 165 193]/255,[45 167 111]/255,[164 146 43]/255,...
        [246 102 126]/255,[197 111 242]/255};
for i=1:length(outputPaths)
    nuisanceSeries = [confoundTSVs{i}.global_signal confoundTSVs{i}.white_matter confoundTSVs{i}.csf...
        confoundTSVs{i}.framewise_displacement confoundTSVs{i}.dvars_std];
    carpetPlotPath = fpp.bids.changeName(outputPaths{i},[],[],'carpetplot','.png');
    fpp.util.carpetPlot(outputPaths{i},maskNonZeroPath,segmentMaskPaths,nuisanceSeries,...
        nuisanceNames,carpetPlotPath,segmentColors,nuisanceColors);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 11: Compute tSNR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 11, Compute tSNR                          - ' outputNameGeneric]);
if genTSNR
    for i=1:length(outputPaths)
        fpp.util.tsnrMap(outputPaths{i},[],fpp.func.preproc.description(tsnrIntro,steps{i}));
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CLEANUP: Delete midprep files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if deleteMidprep
    fpp.util.deleteImageAndJson(pathsToDelete);
end


fprintf('%s\n\n',['Finished preproc for ' outputNameGeneric]);

end