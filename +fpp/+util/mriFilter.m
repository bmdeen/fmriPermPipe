
% Function to filter fMRI dataset (NIFTI or CIFTI), using a Hamming-window
% based FIR filter and MATLAB's filtfilt.
%
% fpp.util.mriFilter(inputPath,outputPath,filtCutoff,filtType,filtOrder,tr)
%
% Arguments:
% - inputPath (string): path to input image
% - outputPath (string): path to output image
% - filtCutoff (one or two element vector): cutoff frequencies in Hz
% - filtType (string, optional): high, low, bandpass, or bandstop (default: 
%   high for one frequency, bandpass for two)
% - filtOrder (scalar, optional): FIR filter order (default: 60)
% - tr (scalar, optional): repetition time in seconds (default: determine
%   from data)

function mriFilter(inputPath,outputPath,filtCutoff,filtType,filtOrder,tr)

if ~exist('tr','var') || isempty(tr)
    tr = fpp.util.checkMRIProperty('tr',inputPath);
    if isempty(tr)
        error('TR of input data could not be determined.');
    end
end
if ~exist('filtType','var') || isempty(filtType)
    filtType = [];
end
if ~exist('filtOrder','var') || isempty(filtOrder)
    filtOrder = 60;
end

[inputMat,hdr] = fpp.util.readDataMatrix(inputPath);
inputMat = inputMat';
outputMat = fpp.util.firFilter(inputMat,1/tr,filtCutoff,filtType,filtOrder) + repmat(mean(inputMat),[size(inputMat,1) 1]);
outputMat = outputMat';
fpp.util.writeDataMatrix(outputMat,hdr,outputPath);

end