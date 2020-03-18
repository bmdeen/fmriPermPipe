
% Wrapper for Freesurfer's MRIwrite function, which results in the
% appropriate x/y coordinate labeling

function err = MRIwrite2(mri,fstring,datatype)

if(nargin == 2) datatype = 'float'; end

mri.vol = permute(mri.vol,[2 1 3:length(size(mri.vol))]);
MRIwrite(mri,fstring,datatype);

end