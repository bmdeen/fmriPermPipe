%
% fpp.anat.preprocCoronal(inputCoronalPaths,preprocT2Path,varargin)
% 
% Preprocesses high-inplane-resolution coronal T2-weighted anatomical
% images, intended for medial temporal lobe structural analysis. Register
% and average multiple images if they exist, register image to individual
% space. A typical resolution for these images would be .5x.5x1.5mm.
%
% Arguments:
% - inputCoronalPaths (cell array of strings): paths to input coronal images
% - preprocT2Path (string): path to preprocessed T2 image in individual
%   space
% 
% Variable arguments:
% - overwrite (boolean): whether to overwrite existing outputs
%

function preprocCoronal(inputCoronalPaths,preprocT2Path,varargin)
%
% TODO
% - Add registration to functional space, ASHS segmentation, conversion of
%   MTL ROIs across spaces

% Check system configuration
fpp.util.checkConfig;

% Basic parameters
overwrite = 0;                  % Whether to overwrite output

% Edit variable arguments.  Note: optInputs checks for proper input.
varArgList = {'overwrite'};
for i=1:length(varArgList)
    argVal = fpp.util.optInputs(varargin,varArgList{i});
    if ~isempty(argVal)
        eval([varArgList{i} ' = argVal;']);
    end
end

% Define wrapper function for fpp.bids.removeBidsDir (for cellfun functionality)
removeBidsDir = @(x) fpp.bids.removeBidsDir(x);

% Define input/output directories
[inputDir,inputName,~] = fpp.util.fileParts(inputCoronalPaths{1});
inputNameGeneric = strrep(fpp.bids.changeName(inputName,'run',[]),'_T1w','');
if isempty(inputDir), inputDir = pwd; end
[anatPreprocDir,~,~] = fpp.util.fileParts(preprocT2Path);
if isempty(anatPreprocDir), inputDir = pwd; end
if strcmp(anatPreprocDir(end),'/'), anatPreprocDir = anatPreprocDir(1:end-1); end
% anatResolution = fpp.bids.checkNameValue(preprocT2Path,'res');

% Define registration outputs
coronal2IndividualXfm = fpp.bids.changeName(preprocT2Path,{'desc','space','res','from','to','mode'},...
    {[],[],[],'nativeCoronal','individual','image'},'xfm','.mat');
individual2CoronalXfm = fpp.bids.changeName(coronal2IndividualXfm,{'from','to'},{'individual','nativeCoronal'});
preprocT2InCoronalPath = fpp.bids.changeName(preprocT2Path,{'space','res'},{'nativeCoronal',''});
preprocT1Path = fpp.bids.changeName(preprocT2Path,[],[],'T1w');
preprocT1InCoronalPath = fpp.bids.changeName(preprocT1Path,{'space','res'},{'nativeCoronal',''});

% Check if output exists.
if exist(preprocT2InCoronalPath,'file') && ~overwrite
    return;
end

% Copy raw data/metadata, convert to .nii.gz if necessary
for i=1:length(inputCoronalPaths)
    [~,inputCoronalNames{i},~] = fpp.util.fileParts(inputCoronalPaths{i});
    outputCoronalPaths{i} = [anatPreprocDir '/' fpp.bids.changeName(inputCoronalNames{i},'acq','HighResCoronal','T2w','.nii.gz')];
    fpp.util.copyImageAndJson(inputCoronalPaths{i},outputCoronalPaths{i},'mri');
    fpp.bids.jsonChangeValue(outputCoronalPaths{i},{'Description','RawSources'},...
        {'Raw data copied to derivative directory.',fpp.bids.removeBidsDir(inputCoronalPaths{i})});
end
inputCoronalPathsRaw = outputCoronalPaths;

pathsToDelete = inputCoronalPathsRaw;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 1: Register all coronal anatomicals to first run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 1, Register and average anatomicals       - ' inputNameGeneric]);
inputCoronalPaths = outputCoronalPaths;
outputCoronalPath = fpp.bids.changeName(inputCoronalPaths{1},{'run','space'},{[],'nativeCoronal'});
addCmd = '';
for i=2:length(inputCoronalPaths)
    outputCoronalPaths{i} = fpp.bids.changeName(inputCoronalPaths{i},'space','nativeCoronal');
    xfmThisRun2FirstRun = fpp.bids.changeName(inputCoronalPaths{i},{'run','from','to','mode'},...
    {'',['CoronalIm' fpp.util.numPad(i,2)],'nativeCoronal','image'},'xfm','.mat');
    fpp.fsl.flirt(inputCoronalPaths{i},inputCoronalPaths{1},xfmThisRun2FirstRun,outputCoronalPaths{i},...
        'dof',6,'interp','sinc');
    fpp.bids.jsonChangeValue(outputCoronalPaths{i},{'Description','Sources'},...
        {'Raw data, registered to first anatomical with sinc interpolation.',...
        fpp.bids.removeBidsDir(inputCoronalPaths{i})});
    addCmd = [addCmd '-add ' outputCoronalPaths{i} ' '];
    pathsToDelete = [pathsToDelete outputCoronalPaths{i} xfmThisRun2FirstRun];
end
if length(inputCoronalPaths)>1
    addCmd = [addCmd '-div ' int2str(length(inputCoronalPaths))];
    fpp.fsl.maths(inputCoronalPaths{1},addCmd,outputCoronalPath);
    fpp.bids.jsonChangeValue(outputCoronalPath,{'Description','Sources','RawSources'},...
        {'Raw data, averaged across runs.',cellfun(removeBidsDir,outputCoronalPaths,'UniformOutput',false),...
        cellfun(removeBidsDir,inputCoronalPathsRaw,'UniformOutput',false)});
else
    fpp.util.copyImageAndJson(inputCoronalPaths{1},outputCoronalPath,'mri');
end



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% STEP 2: Registration to individual space
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 2, Register to individual space           - ' inputNameGeneric]);
fpp.fsl.flirt(outputCoronalPath,preprocT2Path,coronal2IndividualXfm,[],'dof',6);   % Initial FLIRT-based xfm
fpp.fsl.invertXfm(coronal2IndividualXfm,individual2CoronalXfm);
fpp.fsl.moveImage(preprocT2Path,outputCoronalPath,preprocT2InCoronalPath,individual2CoronalXfm,'interp','sinc');
fpp.fsl.moveImage(preprocT1Path,outputCoronalPath,preprocT1InCoronalPath,individual2CoronalXfm,'interp','sinc');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CLEANUP: Delete unneeded intermediate images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(pathsToDelete)
    fpp.util.system(['rm -rf ' pathsToDelete{i}]);
    if exist(fpp.bids.jsonPath(pathsToDelete{i}),'file')
        fpp.util.system(['rm -rf ' fpp.bids.jsonPath(pathsToDelete{i})]);
    end
end

end