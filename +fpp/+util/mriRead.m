
% Wrapper for Freesurfer's MRIread function, which swaps the first two
% dimensions so that first matlab dimension corresponds to first image
% dimension.
%
% mri = fpp.util.mriRead(fstring,headeronly)

function mri = mriRead(fstring,headeronly)

if(exist('headeronly')~=1) headeronly = 0; end

mri = MRIread(fstring,headeronly);
mri.vol = permute(mri.vol,[2 1 3:length(size(mri.vol))]);

end