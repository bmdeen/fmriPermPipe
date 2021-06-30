%
% fpp.anat.preprocCoronal(inputCoronalPaths,preprocT2Path,varargin)
% 
% Preprocesses high-inplane-resolution coronal T2-weighted anatomical
% images, intended for medial temporal lobe structural analysis. Register
% and average multiple images if they exist, register image to individual
% space, and generate ASHS-based hippocampal/MTL segmentation. A typical
% resolution for input images would be .5x.5x1.5mm.
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

% Check system configuration
fpp.util.checkConfig;

% Constants: ASHS segmentation info
segInd = {1,2:4,8,10,11:12,13};
segNames = {'CA1','CA23DG','Sub','ERC','PRC','PHC'};

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

% Define ASHS template directory
[fppFuncDir,~,~]		= fileparts(mfilename('fullpath'));			% path to the directory containing this script
tmp = dir([fppFuncDir '/../../data']);
dataDir = tmp(1).folder;
ashsDir = [dataDir '/ashs_atlas_upennpmc_20170810'];

% Define input/output directories
[inputDir,inputName,~] = fpp.util.fileParts(inputCoronalPaths{1});
inputNameGeneric = strrep(fpp.bids.changeName(inputName,'run',[]),'_T1w','');
if isempty(inputDir), inputDir = pwd; end
[anatPreprocDir,~,~] = fpp.util.fileParts(preprocT2Path);
if isempty(anatPreprocDir), inputDir = pwd; end
if strcmp(anatPreprocDir(end),'/'), anatPreprocDir = anatPreprocDir(1:end-1); end
subjID = fpp.bids.checkNameValue(inputName,'sub');
anatResolution = fpp.bids.checkNameValue(preprocT2Path,'res');

% Define registration outputs
coronal2IndividualXfm = fpp.bids.changeName(preprocT2Path,{'desc','space','res','from','to','mode'},...
    {[],[],[],'nativeCoronal','individual','image'},'xfm','.mat');
individual2CoronalXfm = fpp.bids.changeName(coronal2IndividualXfm,{'from','to'},{'individual','nativeCoronal'});
individual2CoronalXfmITK = fpp.bids.changeName(individual2CoronalXfm,'desc','ITKFormat');
preprocT2InCoronalPath = fpp.bids.changeName(preprocT2Path,{'space','res'},{'nativeCoronal',''});
preprocT1Path = fpp.bids.changeName(preprocT2Path,[],[],'T1w');
preprocT1InCoronalPath = fpp.bids.changeName(preprocT1Path,{'space','res'},{'nativeCoronal',''});

% Check if outputs exist.
outputSegIndivPath = fpp.bids.changeName(preprocT2Path,'desc','ashs','dseg');
if exist(preprocT2InCoronalPath,'file') && exist(outputSegIndivPath,'file') && ~overwrite
    return;
end

% Copy raw data/metadata, convert to .nii.gz if necessary
for i=1:length(inputCoronalPaths)
    [~,inputCoronalNames{i},~] = fpp.util.fileParts(inputCoronalPaths{i});
    outputCoronalPaths{i} = [anatPreprocDir '/' fpp.bids.changeName(inputCoronalNames{i},...
        'acq','HighResCoronal','T2w','.nii.gz')];
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
    fpp.bids.jsonChangeValue(outputCoronalPath,{'Description','Sources','RawSources','SpatialRef'},...
        {'Raw data, averaged across runs.',cellfun(removeBidsDir,outputCoronalPaths,'UniformOutput',false),...
        cellfun(removeBidsDir,inputCoronalPathsRaw,'UniformOutput',false),removeBidsDir(outputCoronalPath)});
else
    fpp.util.copyImageAndJson(inputCoronalPaths{1},outputCoronalPath,'mri');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 2: Registration to individual space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 2, Register to individual space           - ' inputNameGeneric]);
fpp.fsl.flirt(outputCoronalPath,preprocT2Path,coronal2IndividualXfm,[],'dof',6);   % Initial FLIRT-based xfm
fpp.fsl.invertXfm(coronal2IndividualXfm,individual2CoronalXfm);
% fpp.fsl.moveImage(preprocT2Path,outputCoronalPath,preprocT2InCoronalPath,individual2CoronalXfm,'interp','sinc');
fpp.fsl.moveImage(preprocT1Path,outputCoronalPath,preprocT1InCoronalPath,individual2CoronalXfm,'interp','sinc');
fpp.bids.jsonChangeValue(preprocT2InCoronalPath,'Description','Raw data, averaged, in NativeCoronal space.');
fpp.bids.jsonChangeValue(preprocT1InCoronalPath,'Description','Raw data, averaged, in NativeCoronal space.');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 3: ASHS-based MTL segmentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist(ashsDir,'dir')
    fprintf('%s\n',['Step 3, Medial temporal lobe segmentation      - ' inputNameGeneric]);
    hemis = {'L','R',''};
    hemiPrefix = {'L_','R_',''};
    outputSegDir = fpp.bids.changeName(outputCoronalPath,'acq','','ashs','');
    outputSegPath = fpp.bids.changeName(outputCoronalPath,{'desc','acq'},{'ashs',''},'dseg');
    % Note: ASHS requires native-space T1 image, not MNI-aligned
    nativeT1Path = fpp.bids.changeName(preprocT1Path,{'desc','res','space'},{'','','nativeT1w'});
    % Perform ASHS segmentation
    if ~exist(outputSegDir,'dir')       %%% TEMPORARY FOR DEBUGGING
        fpp.util.system(['$ASHS_ROOT/bin/ashs_main.sh -I sub-' subjID ' -a ' ashsDir...
            ' -g ' nativeT1Path ' -f ' outputCoronalPath ' -w ' outputSegDir]);
    end
    % Copy and rename output segmentation
    tmpPath = fpp.bids.changeName(outputSegPath,'desc','ashsTmpPreprocCoronal20920398542');
    fpp.fsl.maths([outputSegDir '/final/sub-' subjID '_right_lfseg_corr_nogray.nii.gz'],...
        '-bin -mul 100',tmpPath);
    fpp.fsl.maths([outputSegDir '/final/sub-' subjID '_right_lfseg_corr_nogray.nii.gz'],...
        ['-add ' tmpPath ' -add ' outputSegDir '/final/sub-' subjID '_left_lfseg_corr_nogray.nii.gz'],...
        outputSegPath);
    fpp.util.deleteImageAndJson(tmpPath);
    fpp.util.system(['cp ' fpp.bids.jsonPath(outputCoronalPath) ' ' fpp.bids.jsonPath(outputSegPath)]);
    fpp.bids.jsonChangeValue(outputSegPath,{'Description','Sources'},{['Medial temporal lobe segmentation file, '...
        'generated by Automatic Segmentation of Hippocampal Subfields.'],removeBidsDir(outputCoronalPath)});
    fpp.wb.command('volume-label-import',outputSegPath,[dataDir '/desc-ashs_lut.txt'],...
        outputSegPath,'-drop-unused-labels');
    % Convert dseg to subregion masks
    for i=1:length(segInd)
        for h=1:length(hemis)
            switch h
                case 1
                    ind = segInd{i};
                case 2
                    ind = segInd{i}+100;
                case 3
                    ind = [segInd{i} segInd{i}+100];
            end
            outputMaskPath = fpp.bids.changeName(outputSegPath,'desc',['ashs' hemis{h} segNames{i}],'mask');
            fpp.util.label2ROI(outputSegPath,ind,outputMaskPath);
            fpp.bids.jsonChangeValue(outputMaskPath,{'Type','Description','Sources'},...
                {'ROI',[hemiPrefix{h} segNames{i} ' mask generated by Automatic Segmentation'...
                ' of Hippocampal Subfields'],outputSegPath});
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 4: Move ASHS ROIs to individual space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist(ashsDir,'dir')
    fprintf('%s\n',['Step 4, Move MTL segments to individual space  - ' inputNameGeneric]);
    fpp.fsl.moveImage(outputSegPath,preprocT2Path,outputSegIndivPath,...
        coronal2IndividualXfm,'interp','nn','isLabel',1);
    % Convert dseg to subregion masks
    for i=1:length(segInd)
        for h=1:length(hemis)
            switch h
                case 1
                    ind = segInd{i};
                case 2
                    ind = segInd{i}+100;
                case 3
                    ind = [segInd{i} segInd{i}+100];
            end
            outputMaskPath = fpp.bids.changeName(outputSegIndivPath,'desc',['ashs' hemis{h} segNames{i}],'mask');
            fpp.util.label2ROI(outputSegIndivPath,ind,outputMaskPath);
            fpp.bids.jsonChangeValue(outputMaskPath,{'Type','Description','Sources'},...
                {'ROI',[hemis{h} '_' segNames{i} ' mask generated by Automatic Segmentation'...
                ' of Hippocampal Subfields'],outputSegIndivPath});
        end
    end
end


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