
% fpp.util.sphericalROI(brainPath,outPath,coords,r,maskPath)
%
% Function to create a spherical ROI around a coordinate.
%
% Arguments:
%   brainPath = path to brain image defining coordinate space (leave empty
%       to use 2mm MNI space)
%   outPath = path to output ROI image
%   coords = image coordinates of ROI center (zero-indexed)
%   r = sphere radius in mm (optional, default 7.5mm)
%   maskPath = path to mask image to constrain ROI (optional)

function sphericalROI(brainPath,outPath,coords,r,maskPath)

radius = 7.5;
tmpName = 'tmp163140398143';

if exist('r','var'), radius = r; end
if exist('maskPath','var') && ~exist(maskPath,'file'), clear maskPath; end
if isempty(brainPath), brainPath = [getenv('FSLDIR') '/data/standard/MNI152_T1_2mm_brain_mask.nii.gz']; end

[outDir,~,~] = fileparts(outPath);
if ~exist(outDir,'dir'), mkdir(outDir); end
tmpPath = [outDir '/' tmpName '.nii.gz'];

fpp.util.system(['fslmaths ' brainPath ' -Tmean -mul 0 -add 1 -roi ' int2str(coords(1)) ' 1 ' int2str(coords(2)) ...
    ' 1 ' int2str(coords(3)) ' 1 0 1 ' tmpPath ' -odt float']);
fpp.util.system(['fslmaths ' tmpPath ' -kernel sphere ' num2str(radius) ...
    ' -fmean -thr .0000001 -bin ' outPath]);
if exist('maskPath','var')
    fpp.util.system(['fslmaths ' outPath ' -mul ' maskPath ' ' outPath]);
end
fpp.util.system(['rm -rf ' tmpPath]);


end