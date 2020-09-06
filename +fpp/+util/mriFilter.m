
% Function to filter fMRI dataset, using a Hamming-window based FIR filter
% and MATLAB's filtfilt.
%
% inputPath (string): path to input image
% outputPath (string): path to output image
% filtCutoff (one or two element vector): cutoff frequency(s) in Hz
% tr (scalar): repetition time in seconds
% type (string, optional): low, high, bandpass, or bandstop (default: low
%   for one frequency, bandpass for two)
% filtOrder (scalar, optional): FIR filter order (default: 60)

function mriFilter(inputPath,outputPath,filtCutoff,tr,type,filtOrder)

if ~exist('type','var')
    type = [];
end
if ~exist('filtOrder','var')
    filtOrder = 60;
end

inputData = fpp.util.mriRead(inputPath);
outputData = inputData;

dims = size(inputData.vol);
inputMat = reshape(inputData.vol,[prod(dims(1:3)) dims(4)])';
outputMat = fpp.util.firFilter(inputMat,filtCutoff,1/tr,type,filtOrder) + repmat(mean(inputMat),[size(inputMat,1) 1]);

outputData.vol = reshape(outputMat',dims);
fpp.util.mriWrite(outputData,outputPath);

end