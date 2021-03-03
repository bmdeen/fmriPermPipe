
% Function to compute DVARS for an fMRI dataset

function [output,outputStd] = dvars(inputPath)

dataMat = fpp.util.readDataMatrix(inputPath);
dataDiff = zeros(size(dataMat));
dataDiff(:,2:end) = diff(dataMat,1,2);
output = rms(dataDiff);
outputStd = output/mean(std(dataDiff,0,2));

end