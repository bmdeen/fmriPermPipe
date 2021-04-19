
% Function to resample a NIFTI 3D or 4D volumetric image file to the
% cortical surface.
%
% Currently assumes that surfaces are at native (freesurfer) resolution.
%
% Uses wb_command -volume-to-surface-mapping with -ribbon-constrained
% option.
%
% fpp.func.surfaceResample(inputNiftiPath,inputSurfacePaths,inputSurfaceROIPaths,outputCiftiPath,varargin)
%
% Arguments:
% - inputNiftiPath (string): path to 3D/4D NIFTI volume to resample
% - inputSurfacePaths (cell array): paths to input surf.gii files.
%       inputSurfacePaths{1} = {hemi-L_midthickness.surf.gii, hemi-L_white.surf.gii, hemi-L_pial.surf.gii}
%       inputSurfacePaths{2} = {hemi-R_midthickness.surf.gii, hemi-R_white.surf.gii, hemi-R_pial.surf.gii}
% - inputSurfaceROIPaths (cell array): paths to input cortex ROI files
%       inputSurfaceROIPaths{1} = hemi-L_desc-cortexAtlas_mask.shape.gii
%       inputSurfaceROIPaths{2} = hemi-R_desc-cortexAtlas_mask.shape.gii
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
% - fwhm (scalar): FWHM of Gaussian kernel used for CIFTI smoothing. Set to
%       0 for no smoothing.

function surfaceResample(inputNiftiPath,inputSurfacePaths,inputSurfaceROIPaths,outputCiftiPath,varargin)

% TO ADD:
% - Option to exclude local CoefVar outliers from 4D data (a la HCP pipeline)
% - Related: input cortical GM ROI option
% - Resample subcortical to CIFTI
% - Resample to fsLR, not just native space. Need spherical coordinate
%       systems and registration as variable arguments.
% - Consider adding surface and/or volume dilation, as in HCP pipeline?

% Constants
hemis = {'L','R'};
structures = {'CORTEX_LEFT','CORTEX_RIGHT'};

% Variable arguments
isLabel = 0;
subcortSegPath = [];
premat = [];
referencePath = [];
fwhm = 0;

% Edit variable arguments.  Note: optInputs checks for proper input.
varArgList = {'subcortSegPath','premat','referencePath','fwhm'};
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
else
    fpp.util.copyImageAndJson(inputNiftiPath,tmpNiftiPath);
end

for h=1:2
    tmpGiftiPaths{h} = [outputDir '/' inputName '_hemi-' hemis{h}...
        '_tmpSurfaceResample21093520813502.' giftiType '.gii'];
    if isLabel
        fpp.wb.command('volume-label-to-surface-mapping',tmpNiftiPath,inputSurfacePaths{h}{1},tmpGiftiPaths{h},...
            ['-ribbon-constrained ' inputSurfacePaths{h}{2} ' ' inputSurfacePaths{h}{3}]);
    else
        fpp.wb.command('volume-to-surface-mapping',tmpNiftiPath,inputSurfacePaths{h}{1},tmpGiftiPaths{h},...
            ['-ribbon-constrained ' inputSurfacePaths{h}{2} ' ' inputSurfacePaths{h}{3}]);
%         if isStat     % Not needed if GIFTI files are deleted
%             fpp.wb.command('metric-palette',tmpGiftiPaths{h},'MODE_USER_SCALE',[],...
%                 '-pos-user 2.5 6 -neg-user -2.5 -6 -palette-name FSL -disp-pos true');
%         end
    end
    fpp.wb.command('set-structure',tmpGiftiPaths{h},structures{h});

end

% Combine GIFTI hemispheres into CIFTI file
volumeStr = '';
if ~isempty(subcortSegPath)
    volumeStr = [' -volume ' tmpNiftiPath ' ' subcortSegPath];
end
if strcmp(giftiType,'shape')
    fpp.wb.command('cifti-create-dense-scalar',[],[],outputCiftiPath,...
        ['-left-metric ' tmpGiftiPaths{1} ' -roi-left ' inputSurfaceROIPaths{1}...
        ' -right-metric ' tmpGiftiPaths{2} ' -roi-right ' inputSurfaceROIPaths{2} volumeStr]);
    mapName = fpp.bids.changeName(inputName,{'space','den','res'},{[],[],[]},[],'');
    fpp.wb.command('set-map-names',outputCiftiPath,[],[],['-map 1 ' mapName]);
    if isStat
        fpp.wb.command('cifti-palette',outputCiftiPath,'MODE_USER_SCALE',outputCiftiPath,...
            '-pos-user 2.5 6 -neg-user -2.5 -6 -palette-name FSL -disp-pos true');
    end
elseif strcmp(giftiType,'label')
    fpp.wb.command('cifti-create-label',[],[],outputCiftiPath,...
        ['-left-label ' tmpGiftiPaths{1} ' -roi-left ' inputSurfaceROIPaths{1}...
        ' -right-label ' tmpGiftiPaths{2} ' -roi-right ' inputSurfaceROIPaths{2} volumeStr]);
    mapName = fpp.bids.changeName(inputName,{'space','den','res'},{[],[],[]},[],'');
    fpp.wb.command('set-map-names',outputCiftiPath,[],[],['-map 1 ' mapName]);
else
    fpp.wb.command('cifti-create-dense-timeseries',[],[],outputCiftiPath,...
        ['-left-metric ' tmpGiftiPaths{1} ' -roi-left ' inputSurfaceROIPaths{1}...
        ' -right-metric ' tmpGiftiPaths{2} ' -roi-right ' inputSurfaceROIPaths{2}...
        volumeStr ' -timestep ' num2str(tr)]);
end

% CIFTI smoothing
if fwhm>0
    sigmaSm = fwhm/2.355; % Standard deviation of Gaussian kernel
    outputDesc = fpp.bids.checkNameValue(outputCiftiPath,'desc');
    fwhmStr = ['Sm' strrep(num2str(fwhm),'.','p')];
    outputCiftiSmPath = fpp.bids.changeName(outputCiftiPath,'desc',[outputDesc fwhmStr]);
    fpp.wb.command('cifti-smoothing',outputCiftiPath,[num2str(sigmaSm) ' '...
        num2str(sigmaSm) ' COLUMN'],outputCiftiSmPath,['-left-surface '...
        inputSurfacePaths{1}{1} ' -right-surface ' inputSurfacePaths{2}{1}]);
end

% Delete temporary paths
fpp.util.deleteImageAndJson(tmpNiftiPath);
for h=1:2, fpp.util.deleteImageAndJson(tmpGiftiPaths{h}); end

end