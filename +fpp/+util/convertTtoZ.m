
% fpp.util.convertTtoZ(inputPath,outputPath,dof)
%
% Function to convert a t-statistic map to a z-statistic map

function convertTtoZ(inputPath,outputPath,dof)

[tVec,hdr] = fpp.util.readDataMatrix(inputPath);

zVec = zeros(size(tVec));
zVec(tVec>0) = -icdf('norm',tcdf(-double(tVec(tVec>0)),dof),0,1);
zVec(tVec<0) = icdf('norm',tcdf(double(tVec(tVec<0)),dof),0,1); % Prevents p-values near 1 from rounding to 1

if sum(zVec==Inf)>0
    warning(['Z-stat map was written with infinite values. Cannot convert such '
        'high t-statistics to z-statistics given MATLAB precision limits - ' outputPath]);
end

fpp.util.writeDataMatrix(zVec,hdr,outputPath);

if ~isempty(fpp.bids.getMetadata(inputPath))
    fpp.bids.jsonReconstruct(inputPath,outputPath);
end

end