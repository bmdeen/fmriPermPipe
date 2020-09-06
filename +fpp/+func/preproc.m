
% fpp.func.preproc(inputPaths,outputDir,varargin)
% 
% Preprocesses a single fMRI dataset, including including artifact 
% detection/scrubbing, motion correction, skull-stripping, spatial 
% smoothing, intensity normalization, high-pass temporal filtering 
% (optional), and registration to a functional template image. Should be
% run after convertDCM.
% 
% Example usage:
% fpp.func.preproc({'/pathToData/sub-01_task-faceloc_run-01_echo-1_bold.nii.gz',...
%       '/pathToData/sub-01_task-faceloc_run-01_echo-2_bold.nii.gz'},'/pathToOutput/');
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
% - targetName (string; default='target_func'): name of image in func 
%       directory to register example_func to, without file extension.
% - genTarget (boolean; default=1): whether to generate target_func image
%       from example_func of first run analyzed if it doesn't exist.
% - disdaqs = 3 (integer in [0,Inf); default=3): number of disdaq volumes
%       (at beginning of run) to remove.
% - transCutoff (scalar in [0,Inf]; default=.5): artifact detection cutoff
%       for total translation (mm).
% - rotCutoff (scalar in [0,Inf]; default=.5): artifact detection cutoff
%   for total rotation (degrees).
% - tptsAfter (integer in [0,Inf); default=0): remove this many time points
%       after pairs of volumes with motion.
% - stdCutoff (scalar in (0,Inf); default=3.5): artifact detection cutoff
%       for mean signal z-score
% - multiEcho (boolean; default=0): whether data is multi-echo
% - teVals (vector of values in (0,Inf)): TE values (ms) for multi-echo 
%       data.
% - sliceTimes (vector of values in (0,Inf)): vector of slice acquisition
%       times (s) relative to volume acquisition onset.
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
%%% ADD NEW ARGUMENTS HERE!
%
% SpatialReference info:
%   Standard spaces:
%   freesurfer, fsLR, MNI152NLin6ASym (FSL/HCP), MNI152NLin2009cAsym
%   (fMRIPrep)
%   Individual spaces: fsnative, anat/T1w, orig, session, task
%
%
%
% TODO NEXT:
% - Add TEDANA
% - By default, construct session template by averaging registered SBRef
%   images across given session (use func if SBRef doesn't exist). Make
%   sure it's in LAS orientation, and include skull-stripped version.
% - Allow option to register to a task- and seq-specific template,
%   constructed using only that seq/task.
% - Add spatial smoothing
% - Add temporal filtering (compare matlab/FSL)
% - Add option to delete midprep images
% - Add JSON files for mocotarget, undistorted mocotarget, func template
%
% TODO EVENTUALLY:
% - Save artifact time points in BIDS format
% - Extract artifact time points from 3dDespike
% - Ensure that flirt can compute registrations with a reflection between
%   images, to account for potential orientation changes.
% - Conversion to LAS based on MRIRead2 and MRIWrite2.
% - Method to deal with disdaqs, taking into account BIDS info, adding them
%   as artifact time points. Based on NumberOfVolumesDiscardedByUser
%   and NumberOfVolumesDiscardedByScanner json fields
% - Consider changing "Sources/RawSources" functionality to avoid listing
%   mid-preproc images as sources.
% - Consider moving topup to copyAndCheckSpinEcho. All of the relevant info
%   is there!!
% - For initial brain mask, consider using BET?
% - Change "funcTemplatePath" variable name?
% - Add functionality for multiple spin echo field maps IntendedFor one
%   functional dataset (averaging maps together first)?
% 

function preproc(inputPaths,outputDir,varargin)

addpath([strrep(mfilename('fullpath'),mfilename,'') '/utils']);

% Load/check config variables.
[configError, fslPrefix] = checkConfig;
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
funcTemplateName = '';          % Name of functional template image, required if path is specified
useTaskTemplate = 0;            % Use a task/acq-specific functional template, rather than a session template
newMedian = 10000;              % Target value for intensity normalization (median of data is set to this)
useTedana = 1;                  % Whether to use tedana-based denoising for multi-echo data
undistort = 1;                  % Whether to distortion-correct functional data using blip-up/blip-down method
spinEchoPaths = {};             % Paths to spin echo images with forward/reversed phase-encode direction
spinEchoPhaseEncodeDirections = {}; % Spin echo phase-encode directions, BIDS format (Default: read from JSON file)
funcDataPhaseEncodeDirection = [];  % Functional data phase-encode direction, BIDS format (Default: read from JSON file)
fieldMapParamPath = '';         % Topup field map parameter file (Default: read from JSON file)
useSTC = 1;                     % Whether to use slice timing correction
sliceTimes = [];                % Slice time (s) relative to onset of volume acquisition (Default: read from JSON file)

% Artifact detection parameters
disdaqs = 0;                    % Number of disdaq volumes (at beginning of run) to remove.
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

% Edit variable arguments.  Note: optInputs checks for proper input.
varArgList = {'overwrite','funcTemplatePath','funcTemplateName','fwhm',...
    'tempFilt','teVals','transCutoff','sliceTimes','funcDataPhaseEncodeDirection',...
    'rotCutoff','tptsAfter','stdCutoff','disdaqs','genTSNR','plotResults',...
    'faValue','smThresh','useSTC','useTedana','echoForMoCorr','undistort',...
    'spinEchoPaths','spinEchoPhaseEncodeDirections','useTaskTemplate'};
for i=1:length(varArgList)
    argVal = optInputs(varargin,varArgList{i});
    if ~isempty(argVal)
        eval([varArgList{i} ' = argVal;']);
    end
end

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

% Define output directories
if strcmp(outputDir(end),'/'), outputDir = outputDir(1:end-1); end
funcPreprocDir = [outputDir '/func'];
fmapPreprocDir = [outputDir '/fmap'];

% Check if output exists.
[inputDir,inputName,~] = fpp.util.fileParts(inputPaths{1});
if isempty(inputDir), inputDir = pwd; end
finalOutputPath = fpp.bids.changeName([funcPreprocDir '/' inputName '.nii.gz'],{'desc','echo'},{'preproc',[]});
if exist(finalOutputPath,'file') && ~overwrite
    return;
end

% Generate output directories
if ~exist(funcPreprocDir,'dir'), mkdir(funcPreprocDir); end
if undistort && ~exist(fmapPreprocDir,'dir'), mkdir(fmapPreprocDir); end
bidsBaseDir = fpp.bids.checkBidsDir(outputDir);  % BIDS root directory
if isempty(bidsBaseDir), bidsBaseDir = '123059815310419841xyz'; end     % Hack so that strrep commands below work for non-BIDS input

% Copy raw functional data and metadata, convert to .nii.gz if necessary
for e=1:nEchoes
    [~,inputNames{e},~] = fpp.util.fileParts(inputPaths{e});
    outputPaths{e} = [funcPreprocDir '/' inputNames{e} '.nii.gz'];
    fpp.util.copyImageAndJson(inputPaths{e},outputPaths{e});
end

% Copy and validate spin echo "field map" images, if undistorting
if undistort
    [errorMsg,spinEchoPaths,fieldMapParamPath] = fpp.func.preproc.copyAndCheckSpinEcho(inputPaths{1},...
        fmapPreprocDir,spinEchoPaths,spinEchoPhaseEncodeDirections,funcDataPhaseEncodeDirection,...
        fieldMapParamPath);
    if ~isempty(errorMsg)
        fprintf('%s\n\n',errorMsg);
        return;
    end
end

% If funcTemplatePath isn't specified by user, define it. If it is, check if file exists.
if isempty(funcTemplatePath)
    if useTaskTemplate
        funcTemplatePath = [funcPreprocDir '/' fpp.bids.changeName(inputName,...
            {'desc','run','echo'},{'FuncTemplate','',''},[],'.nii.gz')];
        funcTemplateName = 'task';
    else
        funcTemplatePath = [funcPreprocDir '/' fpp.bids.changeName(inputName,...
            {'desc','task','acq','run','echo'},{'FuncTemplate','','','',''},[],'.nii.gz')];
        funcTemplateName = 'session';
    end
elseif ~exist(funcTemplatePath,'file')
    fprintf('%s\n\n',['ERROR: Registration target ' funcTemplatePath ' does not exist.']);
    return;
elseif isempty(funcTemplateName)
    fprintf('%s\n\n','ERROR: Need to specify funcTemplateName if using non-standard template.');
    return;
end

% Check TR
tr = fpp.util.checkMRIProperty('tr',outputPaths{1});
if isempty(tr)
    fprintf('%s\n\n','ERROR: Could not compute TR of functional data.');
    return;
end

% Check # of volumes
[~,vols] = system([fslPrefix 'fslval ' outputPaths{1} ' dim4']);
vols = str2num(strtrim(vols));

% Create a wrapper function for converting '.nii.gz' to '.json'
convertNiiJson = @(x) strrep(x,'.nii.gz','.json');
% Create a wrapper to remove BIDS base directory from path
removeBidsBaseDir = @(x) strrep(x,[bidsBaseDir '/'],'');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 1: Estimate motion parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 1, Estimate motion parameters             - ' inputName]);
inputPaths = outputPaths;
mcDir = [funcPreprocDir '/' strrep(inputName,'_bold','') '_motion'];
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 2: 3dDespike
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 2, Despike                                - ' inputName]);
for e=1:nEchoes
    outputPaths{e} = fpp.bids.changeName(outputPaths{e},'desc','midprep1despike');
    %if exist(outputPaths{e},'file'), continue; end      %%% TEMPORARY DEBUGGING HACK
    [~,~] = system(['3dDespike ' inputPaths{e}]);
    [~,~] = system(['3dAFNItoNIFTI ' pwd '/despike+orig.BRIK']);
    [~,~] = system(['mri_convert ' pwd '/despike.nii ' outputPaths{e}]);
    jsonReconstruct(convertNiiJson(inputPaths{e}),convertNiiJson(outputPaths{e}));
    jsonChangeValue(convertNiiJson(outputPaths{e}),{'Description','Sources','SkullStripped','SpatialReference'},...
        {'Partially preprocessed data generated by fmriPermPipe, saved after despiking step.',...
        removeBidsBaseDir(inputPaths{e}),false,'orig'});
    system(['rm -rf ' pwd '/despike+orig.BRIK ' pwd '/despike+orig.HEAD '...
        pwd '/despike.nii']);
end
%%%
%%% Add: extract ART time points from 3dDespike
%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 3: Slice timing correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if useSTC
    fprintf('%s\n',['Step 3, Slice timing correct                   - ' inputName]);
    inputPaths = outputPaths;
    [errorMsg,outputPaths] = fpp.func.preproc.sliceTimingCorrect(inputPaths,sliceTimes,tr);
    if ~isempty(errorMsg)
        fprintf('%s\n\n',errorMsg);
        return;
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 4: Generate and unwarp MocoTarget vol
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 4, Generate/unwarp MocoTargetVol          - ' inputName]);
mocoTargetPath = fpp.bids.changeName(outputPaths{echoForMoCorr},'desc','MocoTargetVol');
system(['fslroi ' outputPaths{echoForMoCorr} ' ' mocoTargetPath ' ' int2str(moCorrTargetVolNum) ' 1']);
% Undistorted MocoTarget vol
if undistort
    mocoTargetUndistortedPath = fpp.bids.changeName(outputPaths{echoForMoCorr},'desc','MocoTargetVolUndistorted');
    [errorMsg,topupWarpPath,topupJacobianPath,xfmMocoTarget2SpinEcho,xfmSpinEcho2MocoTarget] = ...
        fpp.func.preproc.undistort(mocoTargetPath,mocoTargetUndistortedPath,spinEchoPaths,fmapPreprocDir,fieldMapParamPath);
    if ~isempty(errorMsg)
        fprintf('%s\n\n',errorMsg);
        return;
    end
end

% If functional template doesn't exist, create it from native func vol
% EDIT TO SWITCH TO SBREF (and switch to a function)
% Default: Take an average of all available SBRef images, registered to 
% each other with spline interpolation, undistorted, and converted to LAS
% space.
if ~exist(funcTemplatePath,'file')
    if undistort
        copyImageAndJson(mocoTargetUndistortedPath,funcTemplatePath);
    else
        copyImageAndJson(mocoTargetPath,funcTemplatePath);
    end
    genTemplate = 1;
end
fpp.bids.jsonChangeValue(convertNiiJson(funcTemplatePath),{'Description','TaskName','TaskDescription',...
    'Instructions','RepetitionTime','AcquisitionDuration','VolumeTiming','DelayTime',...
    'NumberOfVolumesDiscardedByScanner','NumberOfVolumesDiscardedByUser','DelayAfterTrigger'},...
    {'Functional template image, to be used as a registration target.',[],[],[],[],[],[],[],[],[],[]});

% HERE: Convert functional template to LAS, if it isn't.

% Compute initial brain mask for funcTemplate using BET
initMaskPath = fpp.bids.changeName(funcTemplatePath,'desc','FuncTemplateInitMask','mask');
system(['bet2 ' funcTemplatePath ' ' strrep(initMaskPath,'.nii.gz','') ' -f ' num2str(faValue) ' -m -n']);
fpp.bids.jsonReconstruct(convertNiiJson(funcTemplatePath),convertNiiJson(initMaskPath),{'SpatialReference'});
fpp.bids.jsonChangeValue(convertNiiJson(initMaskPath),{'Description','Sources','Type'},...
    {'Initial brain mask for functional template, generated by BET.',removeBidsBaseDir(funcTemplatePath),'Brain'});



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 5: Registration to functional template
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 5, Registration to FuncTemplate           - ' inputName]);
if undistort
    mocoTargetToRegisterPath = mocoTargetUndistortedPath;
else
    mocoTargetToRegisterPath = mocoTargetPath;
end
xfmMocoTarget2FuncTemplate = fpp.bids.changeName(mocoTargetPath,{'desc','from','to','mode'},...
    {'','orig','session','image'},'xfm','.mat');
xfmFuncTemplate2MocoTarget = fpp.bids.changeName(mocoTargetPath,{'desc','from','to','mode'},...
    {'','orig','session','image'},'xfm','.mat');
system([fslPrefix 'flirt -in ' mocoTargetPath ' -ref ' funcTemplatePath ' -omat ' xfmMocoTarget2FuncTemplate ...
    ' -cost corratio -dof 6 -searchrx -180 180 -searchry -180 180 -searchrz -180 180']);
system([fslPrefix 'convert_xfm -omat ' xfmFuncTemplate2MocoTarget ' -inverse ' xfmMocoTarget2FuncTemplate]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 6: Apply motion/distortion correction, and registration to FuncTemplate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 6, Motion/distortion-correct and register - ' inputName]);
inputPaths = outputPaths;
for e=1:nEchoes
    outputPaths{e} = fpp.bids.changeName(outputPaths{e},'desc','midprep3moco');
end
fpp.func.preproc.oneShotMotionDistortionCorrect(inputPaths,outputPaths,mcDir,funcTemplatePath,...
    funcTemplateName,topupWarpPath,topupJacobianPath,xfmMocoTarget2FuncTemplate,...
    xfmSpinEcho2MocoTarget,xfmMocoTarget2SpinEcho);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 7: TEDANA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if useTedana
    fprintf('%s\n',['Step 7, TEDANA                                 - ' inputName]);
    
    
    
    
    % TEDANA
    % NOTE: remove the global signal removal option!
    % File: midprep4tedana
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 7.5: Multi-echo combination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 7.5, Multi-echo combination               - ' inputName]);
if multiEcho
    inputPaths = outputPaths;
    outputPath = fpp.bids.changeName(inputPaths{1},{'echo','desc'},{[],'midprep5optcomb'});
    system(['t2smap -d ' join(inputPaths,' ') ' -e ' sprintf('%f ',teVals) ' --out-dir ' funcPreprocDir]);
    system(['mv ' funcPreprocDir '/desc-full_T2starmap.nii.gz ' strrep(outputPath,'_bold.nii.gz','_T2star.nii.gz')]);
    system(['mv ' funcPreprocDir '/desc-full_S0map.nii.gz ' strrep(outputPath,'_bold.nii.gz','_S0map.nii.gz')]);
    system(['mv ' funcPreprocDir '/desc-optcom_bold.nii.gz ' outputPath]);
    fpp.bids.jsonReconstruct(convertNiiJson(inputPaths{1}),convertNiiJson(outputPath));
    fpp.bids.jsonChangeValue(convertNiiJson(outputPath),{'Description','Sources','EchoTime','EchoNumber'},...
        {'Partially preprocessed data generated by fmriPermPipe, saved after multi-echo combination step.',...
        cellfun(removeBidsBaseDir,inputPaths,'UniformOutput',false),teVals,[]});
else
    outputPath = outputPaths{1};
end
 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 8: Intensity normalization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 8, Intensity normalization                - ' inputName]);
inputPath = outputPath;
outputPath = fpp.bids.changeName(inputPath,'desc',{[],'midprep6intnorm'});
% Normalize image median (also called "grand mean scaling" in SPM)
[~,funcMedian] = system([fslPrefix 'fslstats ' inputPath ' -k ' initMaskPath ' -p 50']);
system([fslPrefix 'fslmaths ' inputPath ' -mul ' num2str(newMedian/funcMedian) ' ' outputPath]);
fpp.bids.jsonReconstruct(convertNiiJson(inputPath),convertNiiJson(outputPath));
fpp.bids.jsonChangeValue(convertNiiJson(outputPath),{'Description','Sources'},...
    {'Partially preprocessed data generated by fmriPermPipe, saved after intensity normalization step.',...
    removeBidsBaseDir(inputPath)});



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 9: Temporal filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 9, Temporal filtering                     - ' inputName]);
if tempFilt
    inputPath = outputPath;
    outputPath = fpp.bids.changeName(inputPath,'desc',{[],'midprep7tempfilt'});
    
    % TEMPORAL FILTERING HERE!!!
    % Compare FSL/matlab-based methods
    
    
    
    
    
    
    
    fpp.bids.jsonReconstruct(convertNiiJson(inputPath),convertNiiJson(outputPath));
    fpp.bids.jsonChangeValue(convertNiiJson(outputPath),{'Description','Sources'},...
        {'Partially preprocessed data generated by fmriPermPipe, saved after temporal filtering step.',...
        removeBidsBaseDir(inputPath)});
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 10: Spatial smoothing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 10, Spatial smoothing                     - ' inputName]);
if fwhm>0
    inputPath = outputPath;
    outputPath = fpp.bids.changeName(inputPath,'desc',{[],'midprep8smooth'});
    
    % Spatial smoothing here!
    % FSL's method relies on brain-masked data. Just use fslmaths? Ideally
    % I won't be using this form of smoothing.
    
    
    
    
    
    
    fpp.bids.jsonReconstruct(convertNiiJson(inputPath),convertNiiJson(outputPath));
    fpp.bids.jsonChangeValue(convertNiiJson(outputPath),{'Description','Sources'},...
        {'Partially preprocessed data generated by fmriPermPipe, saved after spatial smoothing step.',...
        removeBidsBaseDir(inputPath)});
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 11: Rename output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputPath = outputPath;
outputPath = fpp.bids.changeName(inputPath,'desc',{[],'preproc'});
system(['mv ' inputPath ' ' outputPath]);
system(['mv ' convertNiiJson(inputPath) ' ' convertNiiJson(outputPath)]);
fpp.bids.jsonChangeValue(convertNiiJson(outputPath),{'Description'},{'Preprocessed data generated by fmriPermPipe.'});
fprintf('%s\n\n',['Finished preproc for ' inputName]);

end