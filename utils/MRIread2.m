
% Wrapper for Freesurfer's MRIread function, which results in the
% appropriate x/y coordinate labeling

function mri = MRIread2(fstring,headeronly)

if(exist('headeronly')~=1) headeronly = 0; end

mri = MRIread(fstring,headeronly);
mri.vol = permute(mri.vol,[2 1 3:length(size(mri.vol))]);

end