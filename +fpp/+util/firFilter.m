
% Function to filter time series x, using a Hamming-window based FIR filter
% and MATLAB's filtfilt. NOTE: removes mean from time series.
%
% x (matrix): time series input
% filtCutoff (one or two element vector): cutoff frequency(s)
% sr (scalar): sampling rate in same unit as filtCutoff
% type (string, optional): low, high, bandpass, or bandstop
% filtOrder (scalar, optional): FIR filter order

function y = firFilter(x,filtCutoff,sr,type,filtOrder)

if ~exist('type','var') || isempty(type)
    if length(filtCutoff)==1
        type = 'low';
    else
        type = 'bandpass';
    end
end
if ~exist('filtOrder','var')
    filtOrder = 60;
end

[kernel,~] = fir1(filtOrder,2/sr*filtCutoff,type);
%y = conv(x,kernel,'same');
y = filtfilt(kernel,1,x);

end