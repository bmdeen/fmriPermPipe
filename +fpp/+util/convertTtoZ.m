
% Function to convert a t-statistic map to a z-statistic map

function convertTtoZ(inputPath,outputPath,dof)

[tVec,hdr] = fpp.util.readDataMatrix(inputPath);

zVec = zeros(size(tVec));
zVec(tVec>0) = -icdf('norm',tcdf(-tVec(tVec>0),dof),0,1);
zVec(tVec<0) = icdf('norm',tcdf(tVec(tVec<0),dof),0,1); % Prevents p-values near 1 from rounding to 1

fpp.util.writeDataMatrix(zVec,hdr,outputPath);

if ~isempty(fpp.bids.getMetadata(inputPath))
    fpp.bids.jsonReconstruct(inputPath,outputPath);
end

end