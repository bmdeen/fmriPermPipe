
% Function to FDR-correct z-statistical image, output thresholded image.
%
% fdrCorrect(inputPath,outputPath,maskPath,qThresh,tails)

function fdrCorrect(inputPath,outputPath,maskPath,qThresh,tails)

method = 'pdep';    % pdep = Benjamini & Hochberg, dep = Benjamini & Yekutieli

if ~exist('qThresh','var') || ~(isnumeric(qThresh) && isscalar(qThresh) && ...
        qThresh<1 && qThresh>0)
    qThresh = .05;
end

if ~exist('tails','var') || ~(isnumeric(tails) && isscalar(tails) && ismember(tails,[1 2]))
    tails = 2;
end

zMap = MRIread2(inputPath);
mask = MRIread2(maskPath);

zVec = zMap.vol(mask.vol(:)==1);
pVec = zeros(size(zVec));
if tails==1
    pVec = normcdf(-zVec);
else
    pVec(zVec>0) = 2*normcdf(-zVec(zVec>0));
    pVec(zVec<0) = 2*normcdf(zVec(zVec<0));
end

[h,critP,~,adjP] = fdr_bh(pVec,qThresh,method);
disp(['Critical P-value: ' num2str(critP)]);

zVec(h==0) = 0;
zMap.vol(mask.vol(:)==1) = zVec;
MRIwrite2(zMap,outputPath);


end