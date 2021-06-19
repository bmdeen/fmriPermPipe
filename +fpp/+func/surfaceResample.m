
% fpp.func.surfaceResample(inputNiftiPath,inputSurfacePaths,surfaceROIPaths,outputCiftiPath,varargin)
%
% Function to resample a NIFTI 3D or 4D volumetric image file to the
% cortical surface, generating a CIFTI file.
%
% By default, resamples data to fsnative resolution. If sphereRegFsLRPaths
% and midthickFsLRPaths are defined, data are resampled to fsLR 32k space.
%
% If subcortical (volumetric) information is included in the CIFTI file, it
% can be kept in high-resolution individual space, or can be transformed to
% low-resolution individual space by specifying referenceFuncResPath, or
% warped to a standard space by specifying warp and referenceNonlinPath.
%
% The script assumes 3D inputs are shape files, and 4D inputs are
% functional data. Use isShape and isLabel to change this.
%
% Uses wb_command -volume-to-surface-mapping with -ribbon-constrained
% option.
%
% Arguments:
% - inputNiftiPath (string): path to 3D/4D NIFTI volume to resample
% - inputSurfacePaths (cell array): paths to input surf.gii files.
%     inputSurfacePaths{1} = {hemi-L_midthickness.surf.gii, hemi-L_white.surf.gii, hemi-L_pial.surf.gii}
%     inputSurfacePaths{2} = {hemi-R_midthickness.surf.gii, hemi-R_white.surf.gii, hemi-R_pial.surf.gii}
% - surfaceROIPaths (cell array): paths to cortex ROI files, in space of
%       output data (fsnative or fsLR)
%     surfaceROIPaths{1} = hemi-L_desc-cortexAtlas_mask.shape.gii
%     surfaceROIPaths{2} = hemi-R_desc-cortexAtlas_mask.shape.gii
% - outputCiftiPath (string): path to output CIFTI file
%
% Variable arguments:
% - isLabel (boolean, default false): whether data should be treated as a
%       label file (by default, images are treated as metric data).
% - isShape (boolean): whether data should be treated by a shape file.
%       This is default for 3D images; func is default for 4D.
% - subcortSegPath (string): path to subcortical segmentation label image,
%       in individual or standard space. If specified, CIFTI will include
%       subcortical volumetric in addition to surface data.
% - premat (string): FSL-style linear transformation to apply to volumetric
%       data before surface sampling. Should be used if data are not in
%       individual space.
% - referencePath (string): path to reference image in individual space,
%       target of premat registration. Required if premat is specified.
% - referenceFuncResPath (string): Low resolution individual space image.
%       Volumetric data will be downsampled to this space, while surface
%       data will be resampled directly from the high-resolution reference.
%       Mutually exclusive with referencePathNonlinPath.
% - outputNiftiPath (string): path to volumetric output image.
% - warp (string): path to warp file specifying transformation from
%       individual to standard image
% - referenceNonlinPath (string): path to standard image to warp to.
%       Mutually exclusive with referenceFuncResPath.
% - sphereRegFsLRPaths (cell array of strings): paths to spherical
%       registration files (input coordinates registered to fsLR 32k).
%     sphereRegFsLRPaths{1} = hemi-L_desc-reg2fsLR_sphere.surf.gii
%     sphereRegFsLRPaths{2} = hemi-R_desc-reg2fsLR_sphere.surf.gii
% - midthickFsLRPaths (cell array of strings): paths to midthickness
%       surfaces resampled to fsLR 32k.
% - surfDilation (scalar >= 0): distance (mm) to dilate CIFTI surface data.
%       Intended to remove zero values in output.
% - volDilation (scalar >= 0): distance (mm) to dilate CIFTI volume data.
% - fwhm (scalar): FWHM of Gaussian kernel used for CIFTI smoothing. Set to
%       0 for no smoothing.
% - mapNames (cell array of strings): map names for output label or scalar
%      CIFTI file

function surfaceResample(inputNiftiPath,inputSurfacePaths,surfaceROIPaths,outputCiftiPath,varargin)

% TO ADD:
% - Deal with volume sampling for data resampled to fsLR. Option to
%       downsample data to 2mm, or nonlinear registration to MNI.
% - Option to exclude local CoefVar outliers from 4D data (a la HCP pipeline)?
% - Related: input cortical GM ROI option
% - Resample subcortical to CIFTI properly. Need to modify labels of input
%       segmentation.

% Constants
hemis = {'L','R'};
structures = {'CORTEX_LEFT','CORTEX_RIGHT'};

% FPP data directory
[fppFuncDir,~,~]		= fileparts(mfilename('fullpath'));			% path to the directory containing this script
tmp = dir([fppFuncDir '/../../data']);
dataDir = tmp(1).folder;

% Variable arguments
isLabel = 0;
isShape = 1;
subcortSegPath = [];
premat = [];
referencePath = [];
fwhm = 0;
surfDilation = 0;
volDilation = 0;
outputNiftiPath = '';
referenceFuncResPath = [];
warp = '';
referenceNonlinPath = [];
sphereRegFsLRPaths = {};
midthickFsLRPaths = {};
mapNames = {};

% Edit variable arguments.  Note: optInputs checks for proper input.
varArgList = {'subcortSegPath','premat','referencePath','fwhm','isLabel',...
    'surfDilation','volDilation','outputNiftiPath','sphereRegFsLRPaths',...
    'midthickFsLRPaths','referenceFuncResPath','warp','referenceNonlinPath',...
    'mapNames','isShape'};
for i=1:length(varArgList)
    argVal = fpp.util.optInputs(varargin,varArgList{i});
    if ~isempty(argVal)
        eval([varArgList{i} ' = argVal;']);
    end
end

if ~isempty(premat) && isempty(referencePath)
    error('If premat is specified, referencePath (registration target image) must also be specified.');
elseif isempty(premat) && ~isempty(referencePath)
    error('If referencePath is specified, premat (registration xfm.mat file) must also be specified.');
end
if ~isempty(warp) && isempty(referenceNonlinPath)
    error('If warp is specified, referenceNonlinPath (registration target image) must also be specified.');
elseif isempty(warp) && ~isempty(referenceNonlinPath)
    error('If referenceNonlinPath is specified, warp (registration warp.nii.gz file) must also be specified.');
end

% Load fsLR sphere files, if needed
regFsLR = 0;
if ~isempty(sphereRegFsLRPaths) && ~isempty(midthickFsLRPaths)
    regFsLR = 1;
    [fppFuncDir,~,~]		= fileparts(mfilename('fullpath'));			% path to the directory containing this script
    tmp = dir([fppFuncDir '/../../data']);
    dataDir = tmp(1).folder;
    for h=1:2
        fsLRSpherePaths{h} = [dataDir '/hemi-' hemis{h} '_space-fsLR_den-32k_sphere.surf.gii'];
        outputMidthickPaths{h} = midthickFsLRPaths{h};
    end
else
    for h=1:2
        outputMidthickPaths{h} = inputSurfacePaths{h}{1};
    end
end

[outputDir,~,~] = fpp.util.fileParts(outputCiftiPath);
[~,inputName,~] = fpp.util.fileParts(inputNiftiPath);

% Check if input is 3D or 4D
giftiType = 'shape';
dims = fpp.util.checkMRIProperty('dims',inputNiftiPath);
tr = fpp.util.checkMRIProperty('tr',inputNiftiPath);
if length(dims)>3 && ~isempty(tr) && ~isShape
    giftiType = 'func';
end
if isLabel
    giftiType = 'label';
    fwhm = 0;
    interpStr = 'nn';
else
    interpStr = 'trilinear';
end
if length(dims)==3, dims(4) = 1; end

% Check if input is a t- or z-statistical map
isStat = 0;
if contains(inputNiftiPath,{'tstat.nii','zstat.nii'}) && ~isLabel, isStat = 1; end

% Register image to high-res individual space, if necessary
tmpNiftiPath = [outputDir '/' inputName '_tmpSurfaceResample21093520813502.nii.gz'];
if ~isempty(premat)
    fpp.fsl.moveImage(inputNiftiPath,referencePath,tmpNiftiPath,premat,'interp',interpStr);
    if ~isempty(outputNiftiPath) && isempty(referenceFuncResPath) && isempty(referenceNonlinPath)
        fpp.util.copyImageAndJson(tmpNiftiPath,outputNiftiPath);
    end
else
    fpp.util.copyImageAndJson(inputNiftiPath,tmpNiftiPath);
end
if isLabel  % For labels, re-import label table text file to resampled output
    tmpLUTPath = [outputDir '/' inputName '_tmpSurfaceResample21093520813502_lut.txt'];
    fpp.wb.command('volume-label-export-table',inputNiftiPath,'1',tmpLUTPath);
    fpp.wb.command('volume-label-import',tmpNiftiPath,tmpLUTPath,tmpNiftiPath);
end

% Convert to low-res individual space, if necessary
tmpNiftiVolumePath = [outputDir '/' inputName '_tmpSurfaceResample21093520813502_volume.nii.gz'];
if ~isempty(referenceFuncResPath)
    if isempty(premat)
        premat = [dataDir '/eye.mat'];  % Use identity matrix if converting from high-res to low-res individual
    end
    fpp.fsl.moveImage(inputNiftiPath,referenceFuncResPath,tmpNiftiVolumePath,premat,'interp',interpStr);
elseif ~isempty(referenceNonlinPath) && ~isempty(warp)
    fpp.fsl.moveImage(inputNiftiPath,referenceNonlinPath,tmpNiftiVolumePath,...
        premat,'interp',interpStr,'warp',warp);
else
    tmpNiftiVolumePath = tmpNiftiPath;
end

% Copy volumetric output to outputNiftiPath
if ~isempty(outputNiftiPath)
    fpp.util.copyImageAndJson(tmpNiftiVolumePath,outputNiftiPath);
end

% Resample volume (high-res individual) to surface
for h=1:2
    tmpGiftiPaths{h} = [outputDir '/' inputName '_hemi-' hemis{h}...
        '_tmpSurfaceResample21093520813502.' giftiType '.gii'];
    if isLabel
        fpp.wb.command('volume-label-to-surface-mapping',tmpNiftiPath,inputSurfacePaths{h}{1},tmpGiftiPaths{h},...
            ['-ribbon-constrained ' inputSurfacePaths{h}{2} ' ' inputSurfacePaths{h}{3}]);
    else
        fpp.wb.command('volume-to-surface-mapping',tmpNiftiPath,inputSurfacePaths{h}{1},tmpGiftiPaths{h},...
            ['-ribbon-constrained ' inputSurfacePaths{h}{2} ' ' inputSurfacePaths{h}{3}]);
    end
    fpp.wb.command('set-structure',tmpGiftiPaths{h},structures{h});
end

% Resample to fsLR, if specified
if regFsLR
    for h=1:2
        if isLabel
            fpp.wb.command('label-resample',tmpGiftiPaths{h},[sphereRegFsLRPaths{h} ' '...
                fsLRSpherePaths{h} ' BARYCENTRIC'],tmpGiftiPaths{h},'-largest');
        else
            fpp.wb.command('metric-resample',tmpGiftiPaths{h},[sphereRegFsLRPaths{h} ' ' fsLRSpherePaths{h}...
                ' ADAP_BARY_AREA'],tmpGiftiPaths{h},['-area-surfs ' inputSurfacePaths{h}{1} ' ' midthickFsLRPaths{h}]);
        end
    end
end

% Combine GIFTI hemispheres and NIFTI volume into CIFTI file
volumeStr = '';
if ~isempty(subcortSegPath)
    volumeStr = [' -volume ' tmpNiftiVolumePath ' ' subcortSegPath];
end
if strcmp(giftiType,'shape')
    fpp.wb.command('cifti-create-dense-scalar',[],[],outputCiftiPath,...
        ['-left-metric ' tmpGiftiPaths{1} ' -roi-left ' surfaceROIPaths{1}...
        ' -right-metric ' tmpGiftiPaths{2} ' -roi-right ' surfaceROIPaths{2} volumeStr]);
    if isStat
        fpp.wb.command('cifti-palette',outputCiftiPath,'MODE_USER_SCALE',outputCiftiPath,...
            '-pos-user 2.5 6 -neg-user -2.5 -6 -palette-name FSL -disp-pos true');
    end
elseif strcmp(giftiType,'label')
    fpp.wb.command('cifti-create-label',[],[],outputCiftiPath,...
        ['-left-label ' tmpGiftiPaths{1} ' -roi-left ' surfaceROIPaths{1}...
        ' -right-label ' tmpGiftiPaths{2} ' -roi-right ' surfaceROIPaths{2} volumeStr]);
else
    fpp.wb.command('cifti-create-dense-timeseries',[],[],outputCiftiPath,...
        ['-left-metric ' tmpGiftiPaths{1} ' -roi-left ' surfaceROIPaths{1}...
        ' -right-metric ' tmpGiftiPaths{2} ' -roi-right ' surfaceROIPaths{2}...
        volumeStr ' -timestep ' num2str(tr)]);
end
if ismember(giftiType,{'shape','label'})
    if isempty(mapNames)
        for m=1:dims(4)
            mapNames{m} = [fpp.bids.changeName(inputName,{'space','den','res'},{[],[],[]},[],'') '-' int2str(m)];
        end
    end
    for m=1:dims(4)
        fpp.wb.command('set-map-names',outputCiftiPath,[],[],['-map ' int2str(m) ' ' mapNames{m}]);
    end
end

% Dilate, if specified
if surfDilation>0 || volDilation>0
    fpp.wb.command('cifti-dilate',outputCiftiPath,['COLUMN ' num2str(surfDilation)...
        ' ' num2str(volDilation)],outputCiftiPath,['-left-surface ' outputMidthickPaths{1}...
        ' -right-surface ' outputMidthickPaths{2} ' -nearest']);
end

% CIFTI smoothing
if fwhm>0
    sigmaSm = fwhm/2.355; % Standard deviation of Gaussian kernel
    outputDesc = fpp.bids.checkNameValue(outputCiftiPath,'desc');
    fwhmStr = ['Sm' strrep(num2str(fwhm),'.','p')];
    outputCiftiSmPath = fpp.bids.changeName(outputCiftiPath,'desc',[outputDesc fwhmStr]);
    fpp.wb.command('cifti-smoothing',outputCiftiPath,[num2str(sigmaSm) ' '...
        num2str(sigmaSm) ' COLUMN'],outputCiftiSmPath,['-left-surface '...
        outputMidthickPaths{1} ' -right-surface ' outputMidthickPaths{2}]);
end

% Delete temporary paths
fpp.util.deleteImageAndJson(tmpNiftiPath);
for h=1:2, fpp.util.deleteImageAndJson(tmpGiftiPaths{h}); end
if isLabel && exist(tmpLUTPath), fpp.util.system(['rm -rf ' tmpLUTPath]); end

end