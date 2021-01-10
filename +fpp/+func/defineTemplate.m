%
% fpp.func.defineTemplate(inputPath,outputDir,varargin)
%
% Defines a functional template from raw bold/sbref data
%
% Arguments:
% - inputPath (string): path to input bold/sbref image (any one echo)
% - outputDir (string): output directory (BIDS subject/session level).
%    Outputs will be written to outputDir/func

function defineTemplate(inputPath,outputDir,varargin)

% Load/check config variables.
configError = fpp.util.checkConfig;
if ~isempty(configError)
    fprintf('%s\n',configError);
    return;
end

% Basic parameters
overwrite = 0;                  % Whether to overwrite output
useTaskTemplate = 0;            % Use a task/acq-specific functional template, rather than a session template

% Edit variable arguments.  Note: optInputs checks for proper input.
varArgList = {'overwrite','useTaskTemplate'};
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
if ~exist(funcPreprocDir,'dir') mkdir(funcPreprocDir); end

% Check if output exists.
if useTaskTemplate
    finalOutputPath = fpp.bids.changeName([funcPreprocDir '/' inputName '.nii.gz'],...
        {'run','space','desc'},{[],[],'task','template'});
else
    finalOutputPath = fpp.bids.changeName([funcPreprocDir '/' inputName '.nii.gz'],...
        {'task','run','space','desc'},{[],[],'session','template'});
end
if exist(finalOutputPath,'file') && ~overwrite
    return;
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
    fpp.bids.jsonChangeValue(fpp.bids.jsonPath(outputPaths{e}),{'Description','RawSources','SpatialRef'},...
        {'Functional template space.',fpp.bids.removeBidsDir(inputPaths{e}),fpp.bids.removeBidsDir(outputPaths{1})});
end

% For 4D bold images, extract middle image
for e=1:nEchoes
    vols = fpp.util.checkMRIProperty('vols',inputPaths{e});
    if vols>1
        fpp.util.system(['fslroi ' inputPaths{e} ' ' inputPaths{e} ' ' int2str(ceil(vols/2)) ' 1']);
    end
end

end