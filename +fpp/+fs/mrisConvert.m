
% Wrapper for Freesurfer's mris_convert, which copies input metadata when
% filename is changed.
%
% fpp.fs.mrisConvert(inputPath,outputPath[,flagText])
%
% Runs command: mris_convert flagText inputPath outputPath

function mrisConvert(inputPath,outputPath,flagText)

if ~exist('flagText','var') || isempty(flagText)
    flagText = '';
end

[~,inputName,~] = fpp.util.fileParts(inputPath);
[~,outputName,~] = fpp.util.fileParts(outputPath);

cmd = ['mris_convert ' flagText ' ' inputPath ' ' outputPath];
fpp.util.system(cmd);

if ~strcmp(inputName,outputName) && ~isempty(fpp.bids.getMetadata(inputPath))
    fpp.bids.jsonReconstruct(inputPath,outputPath);
end


end