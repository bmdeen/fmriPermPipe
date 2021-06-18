%
% fpp.func.defineTemplate(inputPath,outputDir,varargin)
%
% Defines a functional template from raw bold/sbref data
%
% Arguments:
% - inputPath (string): path to input bold/sbref image (any one echo)
% - outputDir (string): output directory (BIDS subject/session level).
%    Outputs will be written to outputDir/func
%
% Variable arguments:
% - overwrite (boolean; default=0): whether to overwrite files that have
%       already been written by this function.
% - useTaskTemplate (boolean): whether to use a task-specific template
% - undistort (boolean): Whether to distortion-correct functional data
%       using blip-up/blip-down method
% - spinEchoPath (string): path to spin echo image matched in phase-encode
%       direction to input fMRI data
% - topupWarpPath (string): path to undistortion warp image produced by
%       topup
% - topupJacobianPath (string): path to undistortion warp jacobian image
%       produced by topup

function defineTemplate(inputPath,outputDir,varargin)

% Check system configuration
fpp.util.checkConfig;

% Basic parameters
overwrite = 0;                  % Whether to overwrite output
useTaskTemplate = 0;            % Use a task/acq-specific functional template, rather than a session template

% Undistortion parameters
undistort = 1;                  % Whether to distortion-correct functional template using blip-up/blip-down method
spinEchoPath = '';              % Spin echo path with PE dir matched to functional, if not using default (determined from json metadata)
topupWarpPath = '';             % Undistortion warp path, if not using default (determined from json metadata)
topupJacobianPath = '';         % Undistortion warp jacobian path, if not using default (determined from json metadata)

% Edit variable arguments.  Note: optInputs checks for proper input.
varArgList = {'overwrite','useTaskTemplate','undistort','spinEchoPath',...
    'topupWarpPath','topupJacobianPath'};
for i=1:length(varArgList)
    argVal = fpp.util.optInputs(varargin,varArgList{i});
    if ~isempty(argVal)
        eval([varArgList{i} ' = argVal;']);
    end
end

if ~exist(inputPath,'file')
    error('inputPath does not exist.')
end

% Define input/output directories
[inputDir,inputName,~] = fpp.util.fileParts(inputPath);
if isempty(inputDir), inputDir = pwd; end
if strcmp(outputDir(end),'/'), outputDir = outputDir(1:end-1); end
funcPreprocDir = [outputDir '/func'];
fmapPreprocDir = [outputDir '/fmap'];
if ~exist(funcPreprocDir,'dir') mkdir(funcPreprocDir); end

% Check if output exists.
if useTaskTemplate
    funcTemplateSpace = 'task';
    finalOutputPath = fpp.bids.changeName([funcPreprocDir '/' inputName '.nii.gz'],...
        {'run','space','desc'},{[],[],funcTemplateSpace,'template'});
else
    funcTemplateSpace = 'session';
    finalOutputPath = fpp.bids.changeName([funcPreprocDir '/' inputName '.nii.gz'],...
        {'task','run','space','desc'},{[],[],funcTemplateSpace,'template'});
end
if exist(finalOutputPath,'file') && ~overwrite
    return;
end

% Define undistortion warp files, if undistorting
if undistort && (isempty(spinEchoPath) || isempty(topupWarpPath) ||  isempty(topupJacobianPath))
    [spinEchoPath,topupWarpPath,topupJacobianPath] = fpp.func.preproc.checkSpinEcho(inputPath,fmapPreprocDir);
end

% Check if data is multi-echo
echoExp = '_echo-[1-9]_';
if sum(regexpi(inputName,echoExp))>0
    nEchoes = 0;
    for e=1:9   % Determine # of echoes
        if exist(regexprep(inputPath,echoExp,['_echo-' int2str(e) '_']),'file')
            nEchoes = nEchoes+1;
        end
    end
    for e=1:nEchoes
    	inputPaths{e} = regexprep(inputPath,echoExp,['_echo-' int2str(e) '_']);
    end
else
    nEchoes = 1;
    inputPaths{1} = inputPath;
end

% Copy raw data/metadata for each echo, convert to .nii.gz if necessary
for e=1:nEchoes
    [~,inputNames{e},~] = fpp.util.fileParts(inputPaths{e});
    if useTaskTemplate
        outputPaths{e} = [funcPreprocDir '/' fpp.bids.changeName(inputNames{e},{'run','space','desc'},...
            {[],'task','template'},[],'.nii.gz')];
    else
        outputPaths{e} = [funcPreprocDir '/' fpp.bids.changeName(inputNames{e},{'acq','task','run','space','desc'},...
            {[],[],[],'session','template'},[],'.nii.gz')];
    end
    fpp.util.copyImageAndJson(inputPaths{e},outputPaths{e},'mri');
    fpp.bids.jsonChangeValue(outputPaths{e},{'Description','RawSources','SpatialRef'},...
        {'Functional template space.',fpp.bids.removeBidsDir(inputPaths{e}),fpp.bids.removeBidsDir(outputPaths{1})});
end

% For 4D bold images, extract middle image
for e=1:nEchoes
    vols = fpp.util.checkMRIProperty('vols',outputPaths{e});
    if vols>1
        fpp.util.system(['fslroi ' outputPaths{e} ' ' outputPaths{e} ' ' int2str(ceil(vols/2)) ' 1']);
    end
end

% Undistort functional template
if undistort
    xfmFunc2SpinEcho = fpp.bids.changeName(outputPaths{1},{'desc','from','to','mode','echo'},...
        {'',funcTemplateSpace,'SpinEcho','image',[]},'xfm','.mat');
    xfmSpinEcho2Func = fpp.bids.changeName(outputPaths{1},{'desc','from','to','mode','echo'},...
        {'','SpinEcho',funcTemplateSpace,'image',[]},'xfm','.mat');
    fpp.fsl.flirt(outputPaths{1},spinEchoPath,xfmFunc2SpinEcho,[],'cost','corratio',...
        'dof',6,'searchrx',[-90 90],'searchry',[-90 90],'searchrz',[-90 90]);
    fpp.fsl.invertXfm(xfmFunc2SpinEcho,xfmSpinEcho2Func);
    for e=1:nEchoes
        fpp.func.preproc.undistort(outputPaths{e},outputPaths{e},spinEchoPath,...
            topupWarpPath,topupJacobianPath,xfmFunc2SpinEcho,xfmSpinEcho2Func);
    end
    fpp.util.deleteImageAndJson({xfmFunc2SpinEcho,xfmSpinEcho2Func});
end


end