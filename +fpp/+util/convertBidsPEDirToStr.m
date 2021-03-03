
% Function to convert BIDS-format phase-encode direction (e.g. i, -j) to
% string format (e.g. LR, PA).
%
% convertBidsPEDirToStr(inputPath,phaseEncodeDir)
%
% Arguments:
% - inputPath (string): input image
% - phaseEncodeDir (string, optional): BIDS-format phase-encode direction
%   of the image, needed if PE dir isn't specified by JSON metadata

function phaseEncodeDirStr = convertBidsPEDirToStr(inputPath,phaseEncodeDir)

if ~exist('phaseEncodeDir','var') || isempty(phaseEncodeDir)
    jsonData = fpp.bids.getMetadata(inputPath);
    if isfield(jsonData,'PhaseEncodingDirection')
        phaseEncodeDir = jsonData.PhaseEncodingDirection;
    else
        error('phaseEncodeDir must be specified if PE dir is not specified by JSON metadata.');
    end
end

imageOrientation = fpp.util.getImageOrientation(inputPath);
orientationLabels = {'L','R','A','P','S','I'};
orientationLabelsInverted = {'R','L','P','A','I','S'};
switch phaseEncodeDir(1)
    case 'i'
        peDir = 1;
    case 'j'
        peDir = 2;
    case 'k'
        peDir = 3;
end
phaseEncodeDirStr = [orientationLabelsInverted{strcmp(imageOrientation(peDir),...
    orientationLabels)} imageOrientation(peDir)];
if length(phaseEncodeDir)>1
    phaseEncodeDirStr = fliplr(phaseEncodeDirStr);
end

end
