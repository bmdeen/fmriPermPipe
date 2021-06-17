
% thisTimeOut = fpp.util.convertSiemensTimestamp(thisTime)
%
% Function to convert Siemens Vb timestamp to raw seconds (since midnight).
%
% Example input: 101609.885000, as in 10 hrs, 16 mins, 9s, and 885ms.

function thisTimeOut = convertSiemensTimestamp(thisTime)

if ischar(thisTime), thisTime = str2num(thisTime); end

thisTimeMS = mod(thisTime,1);
thisTimeSecs = mod(thisTime,100)-thisTimeMS;
thisTimeMins = (mod(thisTime,10000)-thisTimeMS-thisTimeSecs)/100;
thisTimeHrs = (thisTime-(thisTimeMins*100+thisTimeSecs+thisTimeMS))/10000;
thisTimeOut = thisTimeHrs*3600+thisTimeMins*60+thisTimeSecs+thisTimeMS;   % Time in s

end