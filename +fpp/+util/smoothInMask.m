
% Function to smooth a volumetric brain image within a mask, to avoid
% smoothing values from external structures into the mask.
%
% fpp.util.smoothInMask(inputPath,maskPath,fwhm,outputPath)
%
% Arguments:
% - inputPath (string): path to input dataset
% - maskPath (string): path to mask to smooth within
% - fwhm (scalar): full width at half max of Gaussian smoothing kernel
% - outputPath (string): path to output dataset

function smoothInMask(inputPath,maskPath,fwhm,outputPath)

sigma = fwhm/2.355; % Standard deviation of Gaussian kernel
tmpPath1 = fpp.bids.changeName(inputPath,'desc','tmpInputMaskedSmoothed1309851857428');	% Input image, masked and smoothed
tmpPath2 = fpp.bids.changeName(inputPath,'desc','tmpMaskSmoothed1309851857428');        % Mask image, smoothed
fpp.fsl.maths(inputPath,['-mul ' maskPath ' -s ' num2str(sigma)],tmpPath1);
fpp.fsl.maths(maskPath,['-s ' num2str(sigma)],tmpPath2);
fpp.fsl.maths(tmpPath1,['-div ' tmpPath2 ' -mul ' maskPath],tmpPath1);
fpp.fsl.maths(maskPath,['-mul -1 -add 1 -mul ' inputPath ' -add ' tmpPath1],outputPath);
fpp.util.system(['rm -rf ' tmpPath1 ' ' tmpPath2]);
if exist(fpp.bids.jsonPath(tmpPath1),'file')
    fpp.util.system(['rm -rf ' fpp.bids.jsonPath(tmpPath1)]);
end
if exist(fpp.bids.jsonPath(tmpPath2),'file')
    fpp.util.system(['rm -rf ' fpp.bids.jsonPath(tmpPath2)]);
end

if ~isempty(fpp.bids.getMetadata(inputPath)) && ~strcmp(inputPath,outputPath)
    fpp.bids.jsonReconstruct(inputPath,outoutPath);
end

end