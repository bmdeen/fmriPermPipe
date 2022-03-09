
% fpp.util.convertTtoZ(inputPath,outputPath,dof)
%
% Function to convert a t-statistic map to a z-statistic map

function convertTtoZ(inputPath,outputPath,dof)

[tVec,hdr] = fpp.util.readDataMatrix(inputPath);

zVec = zeros(size(tVec));
zVec(tVec>0) = -icdf('norm',tcdf(-double(tVec(tVec>0)),dof),0,1);
zVec(tVec<0) = icdf('norm',tcdf(double(tVec(tVec<0)),dof),0,1); % Prevents p-values near 1 from rounding to 1

if sum(zVec==Inf | zVec==-Inf)>0
%     warning(['Z-stat map was written with infinite values. Cannot convert such '...
%         'high t-statistics to z-statistics given MATLAB precision limits - ' outputPath]);
    zVec(zVec==Inf) = max(zVec(zVec~=Inf));
    zVec(zVec==-Inf) = min(zVec(zVec~=-Inf));
    warning(['WARNING: Some of the t-statistic values in map ' inputPath ' were too large'...
        ' to convert to z-statistics given MATLAB precision limits. These values have been'...
        ' set to the extrema of the remaining z-statistic values.']);
end

fpp.util.writeDataMatrix(zVec,hdr,outputPath);

if ~isempty(fpp.bids.getMetadata(inputPath))
    fpp.bids.jsonReconstruct(inputPath,outputPath);
end

end