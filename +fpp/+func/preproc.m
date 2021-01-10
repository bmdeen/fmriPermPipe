%
% fpp.func.preproc(inputPaths,outputDir,varargin)
% 
% Preprocesses a single (multi-echo) fMRI dataset, including including
% motion parameter estimation, despiking, slice timing correction, one-shot
% motion and distoration correction and functional template registration,
% TEDANA multi-echo ICA denoising and optimal echo combination, intensity
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
% - faValue (scalar in (0,1); default=.3): fractional intensity threshold
%       for BET2-based brain extraction. Smaller value -> larger mask.
% - disdaqs = 3 (integer in [0,Inf); default=3): number of disdaq volumes
%       (at beginning of run) to remove. [NOT CURRENTLY IMPLEMENTED]
% - transCutoff (scalar in [0,Inf]; default=.5): artifact detection cutoff
%       for total translation (mm).
% - rotCutoff (scalar in [0,Inf]; default=.5): artifact detection cutoff
%   for total rotation (degrees).
% - tptsAfter (integer in [0,Inf); default=0): remove this many time points
%       after pairs of volumes with motion.
% - stdCutoff (scalar in (0,Inf); default=3.5): artifact detection cutoff
%       for mean signal z-score
% - teVals (vector of values in (0,Inf)): TE values (ms) for multi-echo 
%       data (default: read from json)
% - sliceTimes (vector of values in (0,Inf)): vector of slice acquisition
%       times (s) relative to volume acq onset. (default: read from json)
% - tempFilt (boolean; default=0): whether to perform high-pass temporal
%       filtering.
% - filtCutoff (scalar in (0,Inf); default=.01: highpass filter cutoff (Hz)
% - filtOrder (integer in (0,Inf); default=60): Hanning-window FIR filter
%       order.
% - genTSNR (boolean; default=0): whether to generate TSNR maps.
% - plotResults (boolean; default=1): whether to display result plots.
% - useSTC (boolean; default=1): whether to correct for slice timing.
% - useTedana (boolean; default=1): whether to use TEDANA denoising, for
%       multi-echo data.
% - funcTemplatePath (string): path to functional template image to
%   register to (default: determined from input name)
% - funcTemplateSpace (string): space of functional template, needed if
%   using non-standard funcTemplatePath that lacks a BIDS space entity
% - useTaskTemplate (boolean): whether to use a task-specific template
% - funcDataPhaseEncodeDirection (string): phase encode direction of
%   functional data, in BIDS format (e.g. i, j-; default = read from json)
% - spinEchoPaths (cell array of strings): paths to two spin echo "field
%   map" images with opposite phase encoding directions, for distortion
%   correction (default: derive from fmap json, IntendedFor field)
% - spinEchoPhaseEncodeDirections (cell array of strings): phase encode
%   directions of spin echo images, in BIDS format (e.g. i, j-)
% - deleteMidprep (boolean; default=1): whether to delete midprep images
%

function preproc(inputPaths,outputDir,varargin)

% SpatialReference info:
%   Standard spaces:
%   freesurfer, fsLR, MNI152NLin6ASym (FSL/HCP), MNI152NLin2009cAsym
%   (fMRIPrep)
%   Individual spaces: fsnative, individual, native, session, task
%
%
%
% TODO IMMEDIATELY:
% - Make description text augment with each step added (w/ details of
%   processing done), add this as an input to functions. Have displayed
%   number augment automatically as well. Need to edit
%   oneShotMotionDistortionCorrect as well.
% - Add option to delete midprep images. Also to delete internal files,
%   e.g. topup stuff, motion directory.
% - Add JSON files for mocotarget, undistorted mocotarget, func template
% - Add JSON files for additional TEDANA/ts2map outputs
%
% TODO NEXT:
% - At the end of the script, resample to several spaces if desired:
%   - 2mm-res individual space
%   - 2mm-res MNI152Nlin6Asym space
%   - CIFTI: 32k fsLR surface, 2mm individual subcortical
%   - CIFTI: 32k fsLR surface, 2mm MNI152Nlin6Asym subcortical space
% - Consider running steps 8-10 on optcomb (non-TEDANA) output as well
% - Shift field map preproc to separate function
% - Add temporal filtering (compare matlab/FSL)
% - Add spatial smoothing (susan-based, with brain mask)
% - Add additional confounds - global signal, DVARS, trans/rot,
% FramewiseDisplacement, aCompCorr, NonSteadyStateOutlier00,
% ArtifactTimePoint00 - to desc-confounds_regressors.tsv file
% - Potential BIDS confound TSV files:
% -- desc-confounds_regressors.tsv - all regressors!
% -- motion.tsv - just motion regressors
% -- physreg.tsv - physiological regressors
% -- outliers.tsv
% -- mixing.tsv / components.tsv
% - Save artifact time point regressors in BIDS format (time points to remove)
% - Extract artifact time points from 3dDespike
% - Method to deal with disdaqs, taking into account BIDS info, adding them
%   as artifact time points. Based on NumberOfVolumesDiscardedByUser
%   and NumberOfVolumesDiscardedByScanner json fields.
% - If deleteMidprep==1, remove "Sources" field of final output
% - Add suffix option for desc
%
% TODO EVENTUALLY:
% - Add functionality for multiple spin echo field maps IntendedFor one
%   functional dataset (averaging maps together first)?
% - Add physiological regressors!
% - Track versions of all software used in json descriptions.
% 

% Load/check config variables.
configError = fpp.util.checkConfig;
if ~isempty(configError)
    fprintf('%s\n',configError);
    return;
end

% Basic parameters
overwrite = 0;                  % Whether to overwrite output

% General preproc parameters
fwhm = 0;                       % FWHM of spatial smoothing kernel (mm)
smThresh = [];                  % SUSAN brightness threshold
faValue = .3;                   % BET2 fractional intensity threshold. Smaller value -> larger mask.
funcTemplatePath = '';          % Functional template image to register to
funcTemplateSpace = '';         % Space of func template, required if using custom funcTemplatePath with BIDs space entity
useTaskTemplate = 0;            % Use a task/acq-specific functional template, rather than a session template
newMedian = 10000;              % Target value for intensity normalization (median of data is set to this)
useTedana = 1;                  % Whether to use tedana-based denoising for multi-echo data
undistort = 1;                  % Whether to distortion-correct functional data using blip-up/blip-down method
spinEchoPaths = {};             % Paths to spin echo images with forward/reversed phase-encode direction
spinEchoPhaseEncodeDirections = {}; % Spin echo phase-encode directions, BIDS format (Default: read from JSON file)
funcDataPhaseEncodeDirection = [];  % Functional data phase-encode direction, BIDS format (Default: read from JSON file)
fieldMapParamPath = '';         % Topup field map parameter file (Default: read from JSON file)
useSTC = 1;                     % Whether to use slice timing correction
useDespike = 1;                 % Whether to apply AFNI's 3dDepike before STC
sliceTimes = [];                % Slice time (s) relative to onset of volume acquisition (Default: read from JSON file)
deleteMidprep = 1;              % Whether to delete midprep images

% Artifact detection parameters
disdaqs = 0;                    % Number of disdaq volumes (at beginning of run) to remove. [NOT CURRENTLY IMPLEMENTED]
transCutoff = .5;               % Artifact detection cutoff: total translation (mm)
rotCutoff = .5;                 % Artifact detection cutoff: total rotation (degrees)
tptsAfter = 0;                  % Remove this many time points after pairs of volumes with motion
stdCutoff = Inf;                % Artifact detection cutoff: mean signal z-score

% Multi-echo parameters
multiEcho = 0;                  % Whether data is multi-echo
teVals = [];                    % TE values (ms; Default: read from JSON file)
echoForMoCorr = [];             % Which echo to use to estimate motion parameters

% Temporal filtering parameters
tempFilt = 0;                   % Whether to perform temporal filtering
filtCutoff = .01;               % Highpass filter cutoff (Hz)
filtOrder = 60;                 % Hanning-window FIR filter order

% Data generation parameters
genTSNR = 1;                    % Whether to generate TSNR maps
plotResults = 1;                % Whether to display result plots

pathsToDelete = {};             % Midprep paths to delete at the end

% Edit variable arguments.  Note: optInputs checks for proper input.
varArgList = {'overwrite','funcTemplatePath','funcTemplateSpace','fwhm',...
    'tempFilt','teVals','transCutoff','sliceTimes','funcDataPhaseEncodeDirection',...
    'rotCutoff','tptsAfter','stdCutoff','disdaqs','genTSNR','plotResults',...
    'faValue','smThresh','useSTC','useTedana','echoForMoCorr','undistort',...
    'spinEchoPaths','spinEchoPhaseEncodeDirections','useTaskTemplate','useDespike'};
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

% Define input/output directories
[inputDir,inputName,~] = fpp.util.fileParts(inputPaths{1});
if isempty(inputDir), inputDir = pwd; end
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
    [spaceInd,spaceEndInd] = regexp(funcTemplateName,'_space-[a-zA-Z0-9]+_');
    if sum(spaceInd>0)
        funcTemplateSpace = funcTemplateName(spaceInd(1)+7:spaceEndInd(1)-1);
    elseif isempty(funcTemplateSpace)
        error('funcTemplateSpace must be specified when funcTemplatePath lacks a BIDS space entity.')
    end
end

if ~exist(funcTemplatePath,'file')
    error(['Registration target ' funcTemplatePath ' does not exist. Run fpp.func.defineTemplate.']);
end
maskPath = fpp.bids.changeName(funcTemplatePath,{'desc','echo'},{'brain',[]},'mask','.nii.gz');
if ~exist(maskPath,'file')
    error(['Brain mask ' maskPath ' does not exist. Run fpp.func.register.']);
end

% Check if output exists.
finalOutputPath = fpp.bids.changeName([funcPreprocDir '/' inputName '.nii.gz'],...
    {'desc','echo','space'},{'preproc',[],funcTemplateSpace});
if exist(finalOutputPath,'file') && ~overwrite
    return;
end

% Generate output directories
if ~exist(funcPreprocDir,'dir'), mkdir(funcPreprocDir); end
if undistort && ~exist(fmapPreprocDir,'dir'), mkdir(fmapPreprocDir); end

% Copy raw functional data and metadata, convert to .nii.gz if necessary
for e=1:nEchoes
    [~,inputNames{e},~] = fpp.util.fileParts(inputPaths{e});
    outputPaths{e} = [funcPreprocDir '/' inputNames{e} '.nii.gz'];
    fpp.util.copyImageAndJson(inputPaths{e},outputPaths{e},'midprepfmri');
    fpp.bids.jsonChangeValue(fpp.bids.jsonPath(outputPaths{e}),{'Description','RawSources'},...
        {'Raw data copied to derivative directory.',fpp.bids.removeBidsDir(inputPaths{e})});
end
inputPathsRaw = outputPaths;
pathsToDelete = [pathsToDelete outputPaths];

% Copy and validate spin echo "field map" images, if undistorting
if undistort
    [errorMsg,spinEchoPaths,fieldMapParamPath] = fpp.func.preproc.copyAndCheckSpinEcho(inputPaths{1},...
        fmapPreprocDir,spinEchoPaths,spinEchoPhaseEncodeDirections,funcDataPhaseEncodeDirection,...
        fieldMapParamPath);
    if ~isempty(errorMsg)
        error(errorMsg);
    end
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
fprintf('%s\n',['Step 1, Estimate motion parameters             - ' inputName]);
inputPaths = outputPaths;
mcDir = [funcPreprocDir '/' strrep(fpp.bids.changeName(inputName,'echo',int2str(echoForMoCorr)),'_bold','') '_motion'];
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
pathsToDelete = [pathsToDelete mcDir];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 2: 3dDespike
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if useDespike
    fprintf('%s\n',['Step 2, Despike                                - ' inputName]);
    for e=1:nEchoes
        outputPaths{e} = fpp.bids.changeName(outputPaths{e},'desc','midprep1despike');
%         if exist(outputPaths{e},'file'), continue; end      %%% TEMPORARY DEBUGGING HACK
        fpp.util.system(['3dDespike -nomask ' inputPaths{e}]);
        fpp.util.system(['3dAFNItoNIFTI ' pwd '/despike+orig.BRIK']);
        fpp.fs.mriConvert([pwd '/despike.nii'],outputPaths{e});
%         fpp.util.system(['mri_convert ' pwd '/despike.nii ' outputPaths{e}]);
        fpp.bids.jsonReconstruct(inputPaths{e},outputPaths{e},'midprepfmri');
        fpp.bids.jsonChangeValue(outputPaths{e},{'Description','Sources','SkullStripped','SpatialReference'},...
            {'Partially preprocessed data generated by fmriPermPipe, saved after despiking step.',...
            fpp.bids.removeBidsDir(inputPaths{e}),false,'native'});
        fpp.util.system(['rm -rf ' pwd '/despike+orig.BRIK ' pwd '/despike+orig.HEAD '...
            pwd '/despike.nii']);
    end
end
pathsToDelete = [pathsToDelete outputPaths];
%%%
%%% Add: extract ART time points from 3dDespike
%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 3: Slice timing correct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if useSTC
    fprintf('%s\n',['Step 3, Slice timing correct                   - ' inputName]);
    inputPaths = outputPaths;
    [errorMsg,outputPaths] = fpp.func.preproc.sliceTimingCorrect(inputPaths,sliceTimes,tr);
    if ~isempty(errorMsg)
        error(errorMsg);
    end
    fpp.bids.jsonChangeValue(outputPaths{e},{'SkullStripped','SpatialReference'},{false,'native'});   % in case despiking wasn't run
end
pathsToDelete = [pathsToDelete outputPaths];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 4: Generate and unwarp MocoTarget vol
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 4, Generate/unwarp MocoTargetVol          - ' inputName]);
mocoTargetPath = fpp.bids.changeName(outputPaths{echoForMoCorr},'desc','MocoTargetVol');
fpp.util.system(['fslroi ' outputPaths{echoForMoCorr} ' ' mocoTargetPath ' ' int2str(moCorrTargetVolNum) ' 1']);
% Undistorted MocoTarget vol
if undistort
    mocoTargetUndistortedPath = fpp.bids.changeName(outputPaths{echoForMoCorr},'desc','MocoTargetVolUndistorted');
    [errorMsg,topupWarpPath,topupJacobianPath,xfmMocoTarget2SpinEcho,xfmSpinEcho2MocoTarget] = ...
        fpp.func.preproc.undistort(mocoTargetPath,mocoTargetUndistortedPath,spinEchoPaths,fmapPreprocDir,fieldMapParamPath);
    if ~isempty(errorMsg)
        error(errorMsg);
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 5: Register to functional template
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 5, Register to FuncTemplate               - ' inputName]);
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
fpp.fsl.flirt(mocoTargetPath,funcTemplatePathRegTarget,xfmMocoTarget2FuncTemplate,[],'cost','corratio',...
    'dof',6,'searchrx',[-180 180],'searchry',[-180 180],'searchrz',[-180 180]);
fpp.fsl.invertXfm(xfmMocoTarget2FuncTemplate,xfmFuncTemplate2MocoTarget);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 6: Apply motion/distortion correction, and registration to FuncTemplate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 6, Motion/distortion-correct and register - ' inputName]);
inputPaths = outputPaths;
for e=1:nEchoes
    outputPaths{e} = fpp.bids.changeName(outputPaths{e},{'desc','space'},{'midprep3moco',funcTemplateSpace});
end
%if ~exist(outputPaths{end},'file')      % TEMPORARY DEBUGGING HACK
fpp.func.preproc.oneShotMotionDistortionCorrect(inputPaths,outputPaths,funcTemplatePath,...
    funcTemplateSpace,mcDir,topupWarpPath,topupJacobianPath,xfmMocoTarget2FuncTemplate,...
    xfmSpinEcho2MocoTarget,xfmMocoTarget2SpinEcho,echoForMoCorr);
%end
pathsToDelete = [pathsToDelete outputPaths];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 7: TEDANA / Multi-echo Combine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if useTedana
    fprintf('%s\n',['Step 7, TEDANA Multi-echo ICA denoise          - ' inputName]);
    inputPaths = outputPaths;
    outputPath = fpp.bids.changeName(inputPaths{1},{'echo','desc'},{[],'midprep4tedana'});
    outputDescription = 'Partially preprocessed data generated by fmriPermPipe, saved after TEDANA denoising step.';
    %if ~exist(outputPath,'file') %%% TEMPORARY DEBUGGING HACK
    fpp.func.preproc.tedana(inputPaths,outputPath,maskPath,outputDescription,1,teVals);
    %end
elseif multiEcho
    fprintf('%s\n',['Step 7, Multi-echo combine                     - ' inputName]);
    inputPaths = outputPaths;
    outputPath = fpp.bids.changeName(inputPaths{1},{'echo','desc'},{[],'midprep4optcomb'});
    outputDescription = 'Partially preprocessed data generated by fmriPermPipe, saved after multi-echo combination step.';
    fpp.func.preproc.tedana(inputPaths,outputPath,maskPath,outputDescription,0,teVals);
else
    outputPath = outputPaths{1};
end
if useTedana, pathsToDelete = [pathsToDelete outputPath]; end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 8: Intensity normalize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 8, Intensity normalize                    - ' inputName]);
inputPath = outputPath;
outputPath = fpp.bids.changeName(inputPath,'desc','midprep5intnorm');
% Normalize image median (also called "grand mean scaling" in SPM)
[~,funcMedian] = fpp.util.system(['fslstats ' inputPath ' -k ' maskPath ' -p 50']);
funcMedian = str2num(strtrim(funcMedian));
fpp.fsl.fslMaths(inputPath,['-mul ' num2str(newMedian/funcMedian)],outputPath);
% fpp.util.system(['fslmaths ' inputPath ' -mul ' num2str(newMedian/funcMedian) ' ' outputPath]);
fpp.bids.jsonReconstruct(inputPath,outputPath,'midprepfmri');
fpp.bids.jsonChangeValue(outputPath,{'Description','Sources'},...
    {'Partially preprocessed data generated by fmriPermPipe, saved after intensity normalization step.',...
    fpp.bids.removeBidsDir(inputPath)});
pathsToDelete = [pathsToDelete outputPath];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 9: Temporal filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 9, Temporal filter                        - ' inputName]);
if tempFilt
    inputPath = outputPath;
    outputPath = fpp.bids.changeName(inputPath,'desc',{[],'midprep6tempfilt'});
    
    % TEMPORAL FILTERING HERE!!!
    % Compare FSL/matlab-based methods
    
    
    
    
    
    
    
    fpp.bids.jsonReconstruct(inputPath,outputPath,'midprepfmri');
    fpp.bids.jsonChangeValue(outputPath,{'Description','Sources'},...
        {'Partially preprocessed data generated by fmriPermPipe, saved after temporal filtering step.',...
        fpp.bids.removeBidsDir(inputPath)});
    pathsToDelete = [pathsToDelete outputPath];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 10: Spatially smooth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 10, Spatially smooth                      - ' inputName]);
if fwhm>0
    inputPath = outputPath;
    outputPath = fpp.bids.changeName(inputPath,'desc',{[],'midprep7smooth'});
    
    % Spatial smoothing here!
    % FSL's method relies on brain-masked data.
    
    
    
    
    
    
    fpp.bids.jsonReconstruct(inputPath,outputPath,'midprepfmri');
    fpp.bids.jsonChangeValue(outputPath,{'Description','Sources'},...
        {'Partially preprocessed data generated by fmriPermPipe, saved after spatial smoothing step.',...
        fpp.bids.removeBidsDir(inputPath)});
    pathsToDelete = [pathsToDelete outputPath];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 10.5: Rename output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputPath = outputPath;
outputPath = fpp.bids.changeName(inputPath,'desc','preproc');
fpp.util.system(['mv ' inputPath ' ' outputPath]);
fpp.util.system(['mv ' fpp.bids.jsonPath(inputPath) ' ' fpp.bids.jsonPath(outputPath)]);
fpp.bids.jsonReconstruct(inputPath,outputPath,'fmri');
fpp.bids.jsonChangeValue(outputPath,'Description','Preprocessed data generated by fmriPermPipe.');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 11: Compute tSNR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 11, Compute tSNR                          - ' inputName]);
% Multi-echo combine raw input data for tSNR comparison
if multiEcho
    inputPathRaw = fpp.bids.changeName(inputPathsRaw{1},{'echo','desc'},{[],'rawOptcomb'});
    outputDescription = 'Raw data, optimally combined across echoes using tedana.';
    fpp.func.preproc.tedana(inputPathsRaw,inputPathRaw,maskPath,outputDescription,0);
else
    inputPathRaw = inputPathsRaw{1};
end
% Register raw data to funcTemplate
inputPathRaw2FuncTemplate = fpp.bids.changeName(inputPathRaw,'space',funcTemplateSpace);
fpp.fsl.moveImage(inputPathRaw,funcTemplatePath,inputPathRaw2FuncTemplate,xfmMocoTarget2FuncTemplate);
fpp.bids.jsonReconstruct(inputPathRaw,inputPathRaw2FuncTemplate,'fmri');
% fpp.bids.jsonChangeValue(inputPathRaw2FuncTemplate,{'Description','Sources','SpatialReference'},...
%     {'Raw data, optimally combined across echoes, registered to functional template.',...
%     fpp.bids.removeBidsDir(inputPathRaw),fpp.bids.removeBidsDir(funcTemplatePath)});
fpp.bids.jsonChangeValue(inputPathRaw2FuncTemplate,'Description',...
    'Raw data, optimally combined across echoes, registered to functional template.');
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