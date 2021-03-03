
% Function to generate TSV file of outlier/artifact time points, from a set
% of time indices.
%
% outMat = fpp.func.preproc.outlierTSV(ind,vols)

function tsv = outlierTSV(ind,vols)

tsv = struct();

for i=1:length(ind)
    tmp = zeros(vols,1);
    tmp(ind(i),end) = 1;
    tsv.(['motion_outlier' fpp.util.numPad(i-1,2)]) = tmp;
end

end