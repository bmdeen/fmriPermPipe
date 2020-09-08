
% Function to get the orientation (e.g. RAS, LAS) of an MRI image

function orientString = getImageOrientation(inputPath)

[~,fslValOutput] = fpp.util.system(['fslval ' inputPath ' qform_xorient']);
orientString(1) = fslValOutput(findstr(fslValOutput,'-to-')+4);
[~,fslValOutput] = fpp.util.system(['fslval ' inputPath ' qform_yorient']);
orientString(2) = fslValOutput(findstr(fslValOutput,'-to-')+4);
[~,fslValOutput] = fpp.util.system(['fslval ' inputPath ' qform_zorient']);
orientString(3) = fslValOutput(findstr(fslValOutput,'-to-')+4);

end