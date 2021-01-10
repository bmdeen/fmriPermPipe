
% Function to FDR-correct z-statistic image, output thresholded image.
%
% fpp.func.analysis.fdrCorrect(inputPath,outputPath,maskPath,qThresh,tails)

function fdrCorrect(inputPath,outputPath,maskPath,qThresh,tails)

method = 'pdep';    % pdep = Benjamini & Hochberg, dep = Benjamini & Yekutieli

if ~exist('qThresh','var') || ~(isnumeric(qThresh) && isscalar(qThresh) && ...
        qThresh<1 && qThresh>0)
    qThresh = .05;
end

if ~exist('tails','var') || ~(isnumeric(tails) && isscalar(tails) && ismember(tails,[1 2]))
    tails = 2;
end

zMap = fpp.util.mriRead(inputPath);
mask = fpp.util.mriRead(maskPath);

zVec = zMap.vol(mask.vol(:)==1);
pVec = zeros(size(zVec));
if tails==1
    pVec = normcdf(-zVec);
else
    pVec(zVec>0) = 2*normcdf(-zVec(zVec>0));
    pVec(zVec<0) = 2*normcdf(zVec(zVec<0));
end

[h,critP,~,adjP] = fpp.func.analysis.fdrBH(pVec,qThresh,method);
disp(['Critical P-value: ' num2str(critP)]);

zVec(h==0) = 0;
zMap.vol(mask.vol(:)==1) = zVec;
fpp.util.mriWrite(zMap,outputPath);


end