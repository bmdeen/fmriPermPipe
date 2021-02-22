
% Wrapper for Freesurfer's MRIwrite function, which swaps the first two
% dimensions so that first matlab dimension corresponds to first image
% dimension.
%
% err = fpp.util.mriWrite(mri,fstring,datatype)

function err = mriWrite(mri,fstring,datatype)

if(nargin == 2), datatype = 'float'; end

mri.vol = permute(mri.vol,[2 1 3:length(size(mri.vol))]);
err = MRIwrite(mri,fstring,datatype);

end