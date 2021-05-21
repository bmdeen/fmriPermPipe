
% Function to resample a NIFTI 3D or 4D volumetric image file to the
% cortical surface.
%
% By default, resamples data to fsnative resolution. If sphereRegFsLRPaths
% and midthickFsLRPaths are defined, data are resampled to fsLR 32k space.
%
% Uses wb_command -volume-to-surface-mapping with -ribbon-constrained
% option.
%
% fpp.func.surfaceResample(inputNiftiPath,inputSurfacePaths,surfaceROIPaths,outputCiftiPath,varargin)
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
%       label file (by default, 3D images are treated as metric data)
% - subcortSegPath (string): path to subcortical segmentation label image,
%       in individual space. If specified, CIFTI will include subcortical
%       volumetric in addition to surface data.
% - premat (string): FSL-style linear transformation to apply to volumetric
%       data before surface sampling. Should be used if data are not in
%       individual space.
% - referencePath (string): path to reference image in individual space,
%       target ofpremat registration. Required if premat is specified.
% - outputNiftiPath (string): path to volumetric output image.
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

% Variable arguments
isLabel = 0;
subcortSegPath = [];
premat = [];
referencePath = [];
fwhm = 0;
surfDilation = 0;
volDilation = 0;
outputNiftiPath = '';
sphereRegFsLRPaths = {};
midthickFsLRPaths = {};

% Edit variable arguments.  Note: optInputs checks for proper input.
varArgList = {'subcortSegPath','premat','referencePath','fwhm','isLabel',...
    'surfDilation','volDilation','outputNiftiPath','sphereRegFsLRPaths',...
    'midthickFsLRPaths'};
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
if length(dims)>3
    giftiType = 'func';
    tr = fpp.util.checkMRIProperty('tr',inputNiftiPath);
end
if isLabel
    if length(dims)>3
        error('fpp.func.surfaceResample can only handle 3D label files, not 4D.');
        % Note: this is due to set-map-names command.
    end
    giftiType = 'label';
    fwhm = 0;
    interpStr = 'nn';
else
    interpStr = 'trilinear';
end

% Check if input is a t- or z-statistical map
isStat = 0;
if contains(inputNiftiPath,{'tstat.nii','zstat.nii'}) && ~isLabel, isStat = 1; end

% Register image to individual space, if necessary
tmpNiftiPath = [outputDir '/' inputName '_tmpSurfaceResample21093520813502.nii.gz'];
if ~isempty(premat)
    fpp.fsl.moveImage(inputNiftiPath,referencePath,tmpNiftiPath,premat,'interp',interpStr);
    if ~isempty(outputNiftiPath)
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

% Resample volume to surface
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

% Combine GIFTI hemispheres into CIFTI file
volumeStr = '';
if ~isempty(subcortSegPath)
    volumeStr = [' -volume ' tmpNiftiPath ' ' subcortSegPath];
end
if strcmp(giftiType,'shape')
    fpp.wb.command('cifti-create-dense-scalar',[],[],outputCiftiPath,...
        ['-left-metric ' tmpGiftiPaths{1} ' -roi-left ' surfaceROIPaths{1}...
        ' -right-metric ' tmpGiftiPaths{2} ' -roi-right ' surfaceROIPaths{2} volumeStr]);
    mapName = fpp.bids.changeName(inputName,{'space','den','res'},{[],[],[]},[],'');
    fpp.wb.command('set-map-names',outputCiftiPath,[],[],['-map 1 ' mapName]);
    if isStat
        fpp.wb.command('cifti-palette',outputCiftiPath,'MODE_USER_SCALE',outputCiftiPath,...
            '-pos-user 2.5 6 -neg-user -2.5 -6 -palette-name FSL -disp-pos true');
    end
elseif strcmp(giftiType,'label')
    fpp.wb.command('cifti-create-label',[],[],outputCiftiPath,...
        ['-left-label ' tmpGiftiPaths{1} ' -roi-left ' surfaceROIPaths{1}...
        ' -right-label ' tmpGiftiPaths{2} ' -roi-right ' surfaceROIPaths{2} volumeStr]);
    mapName = fpp.bids.changeName(inputName,{'space','den','res'},{[],[],[]},[],'');
    fpp.wb.command('set-map-names',outputCiftiPath,[],[],['-map 1 ' mapName]);
else
    fpp.wb.command('cifti-create-dense-timeseries',[],[],outputCiftiPath,...
        ['-left-metric ' tmpGiftiPaths{1} ' -roi-left ' surfaceROIPaths{1}...
        ' -right-metric ' tmpGiftiPaths{2} ' -roi-right ' surfaceROIPaths{2}...
        volumeStr ' -timestep ' num2str(tr)]);
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