
% Function to generate ROI image from statistical map and search space
% mask, by choosing the top N% of voxels (and optionally applying a
% statistical threshold). Takes NIFTI or CIFTI dscalar inputs.
%
% fpp.util.defineROI(statPath,searchPath,outputPath,varargin)
%
% Arguments:
% - statPath (string or cell array of strings): path to statistical map. If
%       cell array, maps will be averaged first.
% - searchPath (string): path to search space mask image
% - outputPath (string): path to output image
%
% Variable arguments:
% - roiSize (scalar): size of ROI, in % or # of voxels in a search space
% - sizeType ('pct' or 'num'): whether size is measured in % or # of
%       coordinates (voxels or grayordinates)
% - statThresh (scalar): threshold for statistical map, applied after size
% - invertStats (boolean, default = 0): whether to invert statistical map
% - maskPath (string): path to brain mask, to intersect with search space
% - parcPath (string): path to parcellation. If parcPath and parcInds are
%       specified, the search space is defined by the parcellation, using
%       fpp.util.label2ROI.
% - parcInds (numeric vector): indices of parcels to include in the search
%       space.

function defineROI(statPath,searchPath,outputPath,varargin)

roiSize = 5;
sizeType = 'pct';
statThresh = -Inf;
invertStats = 0;
maskPath = '';
parcPath = '';
parcInds = [];

% Define wrapper function for fpp.bids.removeBidsDir (for cellfun functionality)
removeBidsDir = @(x) fpp.bids.removeBidsDir(x);

% Edit variable arguments.  Note: optInputs checks for proper input.
varArgList = {'roiSize','sizeType','statThresh','invertStats','maskPath','parcPath','parcInds'};
for i=1:length(varArgList)
    argVal = fpp.util.optInputs(varargin,varArgList{i});
    if ~isempty(argVal)
        eval([varArgList{i} ' = argVal;']);
    end
end

if ~(strcmpi(sizeType,'pct') || strcmpi(sizeType,'num'))
    error('sizeType must be set to pct or num');
end
if strcmpi(sizeType,'pct') && roiSize>100
    roiSize = 100;
end
if strcmpi(sizeType,'num')
    roiSize = round(roiSize);
end

% Check whether input is NIFTI or CIFTI
[~,~,inputExt] = fpp.util.fileParts(searchPath);
if ~ismember(lower(inputExt),{'.nii.gz','.nii','.dscalar.nii'})
    error('Inputs must be a NIFTI or CIFTI dscalar files.');
end
isCifti = 0;
if strcmpi(inputExt,'.dscalar.nii')
    isCifti = 1;
end
if isCifti
    coord = 'grayordinate';
    imageType = 'cifti';
else
    coord = 'voxel';
    imageType = 'volume';
end

% Generate search space from parcellation, if parcPath/parcInds specified
if ~isempty(parcPath) && ~isempty(parcInds)
    fpp.util.label2ROI(parcPath,parcInds,searchPath)
end

% Average statistical maps, if statPath is a cell array
statAvg = 0;
if iscell(statPath)
    statAvg = 1;
    statPaths = statPath;
    [statDir,statName,~] = fpp.util.fileParts(statPaths{1});
    statPath = [statDir '/' statName '_tmpDefineROI2098315310985' inputExt];
    nStats = length(statPaths);
    
    % Average maps using wb_command volume-math/cifti-math
    weightEquation = '(';
    flagText = '';
    for s=1:nStats
        statStr = ['s' int2str(s)];
        flagText = [flagText ' -var ' statStr ' ' statPaths{s}];
        if s==1
            weightEquation = [weightEquation statStr];
        else
            weightEquation = [weightEquation '+' statStr];
        end
    end
    weightEquation = [weightEquation ')/' num2str(nStats)];
    fpp.wb.command([imageType '-math'],[],weightEquation,statPath,flagText);
end

% Load statistical map and search space
if ~isempty(maskPath)
    maskData = fpp.util.mriRead(maskPath);
    [statVec,~] = fpp.util.readDataMatrix(statPath,maskData.vol);
    [searchVec,hdr] = fpp.util.readDataMatrix(searchPath,maskData.vol);
else
    [statVec,~] = fpp.util.readDataMatrix(statPath);
    [searchVec,hdr] = fpp.util.readDataMatrix(searchPath);
end

searchInd = find(searchVec==1);
searchSize = size(searchInd);
statVec = statVec(searchInd);
if invertStats, statVec = -statVec; end
roiInSearchVec = zeros(size(statVec));    % ROI indices, within search space

% Find top N %/# of voxels in search space
[~,sortInd] = sort(statVec,'descend');
switch lower(sizeType)
    case 'pct'
        roiInSearchVec(sortInd(1:roiSize*searchSize/100)) = 1;
    case 'vox'
        roiInSearchVec(sortInd(1:roiSize)) = 1;
end

% Apply statistical threshold
roiInSearchVec(statVec<statThresh) = 0;

% Define and write final ROI
roiVec = zeros(size(searchVec));
roiVec(searchInd) = roiInSearchVec;
if ~isempty(maskPath)
    fpp.util.writeDataMatrix(roiVec,hdr,outputPath,maskData.vol);
else
    fpp.util.writeDataMatrix(roiVec,hdr,outputPath);
end

% Check if ROI is empty, delete and send warning if so
[~,roiOutputSize] = fpp.util.system(['fslstats ' outputPath ' -V']);
roiOutputSize = str2num(roiOutputSize); roiOutputSize = roiOutputSize(1);
if roiOutputSize==0
    fpp.util.system(['rm -rf ' outputPath]);
    warning(['Output ROI had zero ' coord 's: ' outputPath]);
    return;
elseif roiOutputSize<10
    warning(['Output ROI < 10 ' coord 's: ' outputPath]);
end

% Reconstruct json metadata
if exist(fpp.bids.jsonPath(searchPath),'file')
    fpp.bids.jsonReconstruct(searchPath,outputPath);
    if invertStats, invertStr = ' (inverted)'; else invertStr = ''; end
    if strcmpi(sizeType,'pct'), sizeTypeStr = ['% of ' coord 's'];
    else sizeTypeStr = [' ' coord 's']; end
    outputDescription = ['Region-of-interest file generated by fpp.util.defineROI, using statistical map '...
        removeBidsDir(statPath) invertStr ' and search space ' removeBidsDir(searchPath) '. ROI defined as the top '...
        num2str(roiSize) sizeTypeStr ' in the search space.'];
    if statThresh>-Inf
        outputDescription = [outputDescription(1:end-1) ', followed by a statistical threshold of ' num2str(statThresh) '.'];
    end
    outputSources = cellfun(removeBidsDir,{statPath,searchPath},'UniformOutput',false);
    fpp.bids.jsonChangeValue(outputPath,{'Description','Sources','RawSources'},{outputDescription,outputSources,[]});
end

% Delete temporary, averaged statistic map, if necessary
if statAvg
    system(['rm -rf ' statPath]);
    if exist(fpp.bids.jsonPath(statPath),'file')
        system(['rm -rf ' fpp.bids.jsonPath(statPath)]);
    end
end

end