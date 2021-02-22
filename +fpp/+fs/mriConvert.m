
% Wrapper for Freesurfer's mri_convert, which copies input metadata when
% filename is changed.
%
% fpp.fs.mriConvert(inputPath,outputPath[,flagText])
%
% Runs command: mri_convert flagText inputPath outputPath

function mriConvert(inputPath,outputPath,flagText)

if ~exist('flagText','var')
    flagText = '';
end

[~,inputName,~] = fpp.util.fileParts(inputPath);
[~,outputName,~] = fpp.util.fileParts(outputPath);

cmd = ['mri_convert ' flagText ' ' inputPath ' ' outputPath];
fpp.util.system(cmd);

if ~strcmp(inputName,outputName) && ~isempty(fpp.bids.getMetadata(inputPath))
    fpp.bids.jsonReconstruct(inputPath,outputPath);
end


end