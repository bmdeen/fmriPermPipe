
% Function to filter time series x, using a Hamming-window based FIR filter
% and MATLAB's filtfilt. NOTE: removes mean from time series.
%
% y = fpp.util.firFilter(x,filtCutoff,sr,filtType,filtOrder)
%
% Arguments:
% - x (numeric matrix): time series input. If this is an N-D matrix, time
%   series must be in the first dimension.
% - sr (scalar): sampling rate in same unit as filtCutoff
% - filtCutoff (one or two element vector): cutoff frequencies
% - filtType (string, optional): high, low, bandpass, or bandstop (default: 
%   high for one frequency, bandpass for two)
% - filtOrder (scalar, optional): FIR filter order (default: 60)

function y = firFilter(x,sr,filtCutoff,filtType,filtOrder)

if ~exist('filtType','var') || isempty(filtType)
    if length(filtCutoff)==1
        filtType = 'high';
    else
        filtType = 'bandpass';
    end
end
if ~exist('filtOrder','var') || isempty(filtOrder)
    filtOrder = 60;
end

% Demean and zero-pad input
padSize = filtOrder*3;  % Size of zero-pad: minimum size required by filtfilt
dims = size(x);
x = bsxfun(@minus, x, mean(x));
x = [x; zeros([padSize dims(2:end)])];

% Define Hamming-window-based FIR filter
[kernel,~] = fir1(filtOrder,2/sr*filtCutoff,filtType);

% Perform zero-phase filtering
yPadded = filtfilt(kernel,1,x);

% Remove zero-padding
S.subs = repmat({':'},1,ndims(x));
S.subs{1} = dims(1)+1:dims(1)+padSize;
S.type = '()';
y = subsasgn(yPadded,S,[]);

end