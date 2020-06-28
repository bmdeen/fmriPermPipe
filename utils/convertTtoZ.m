
% Function to convert an t-statistic map to a z-statistic map

function convertTtoZ(inputPath,outputPath,dof)

tMap = MRIread2(inputPath);
tVec = reshape(tMap.vol,[prod(size(tMap.vol)) 1]);

zVec = zeros(size(tVec));
zVec(tVec>0) = -icdf('norm',tcdf(-tVec(tVec>0),dof),0,1);
zVec(tVec<0) = icdf('norm',tcdf(tVec(tVec<0),dof),0,1); % Prevents p-values near 1 from rounding to 1

zMap = tMap;
zMap.vol = reshape(zVec,size(zMap.vol));
MRIwrite2(zMap,outputPath);

end