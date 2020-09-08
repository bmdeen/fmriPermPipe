
% multiEchoCombine(dataPath,teVals,varargin)
%
% Function to take multi-echo data, output T2* and R2* estimates, as well
% as weighted average data.
%
% To add:
% - Brain mask option
%
%
% NEXT MAJOR CHANGE:
% - Switch to ts2map functionality.
% - Bids output naming.
% - Don't split data

function multiEchoCombine(inputPath,teVals,varargin)

% Basic parameters
overwrite = 0;                  % Whether to overwrite output
t2StarPath = '';                % Path of output T2* image
r2StarPath = '';                % Path of output R2* image
outputPath = '';                % Path of output combined data

if ~contains(inputPath,'.nii.gz')
    disp('ERROR: Input image to multiEchoCombine must have .nii.gz extension.');
    return;
end

% Edit variable arguments.  Note: optInputs checks for proper input.
varArgList = {'overwrite','t2StarPath','r2StarPath','outputPath'};
for i=1:length(varArgList)
    argVal = optInputs(varargin,varArgList{i});
    if ~isempty(argVal)
        eval([varArgList{i} ' = argVal;']);
    end
end

[outputDir,inputName,inputExt] = fpp.util.fileParts(inputPath);

if isempty(outputDir)
    disp('ERROR: Please provide full path for input data to multiEchoCombine.');
    return;
end

inputName = [inputName inputExt];
tmpDir = [outputDir '/tmp10913245098143940814'];
mkdir(tmpDir);

if isempty(outputPath)
    outputPath = [outputDir '/' strrep(inputName,'.nii.gz','_MEcomb.nii.gz')];
end
if isempty(t2StarPath)
    t2StarPath = [outputDir '/' strrep(inputName,'.nii.gz','_T2StarEst.nii.gz')];
end
if isempty(r2StarPath)
    r2StarPath = [outputDir '/' strrep(inputName,'.nii.gz','_R2StarEst.nii.gz')];
end

if size(teVals,1)<size(teVals,2), teVals = teVals'; end
numTEs=length(teVals);
for i=1:numTEs
    teSuffices{i} = ['_TE' int2str(i)];
end

% Separate data from different numTEs
[~,vols] = system(['fslval ' inputPath ' dim4']);
vols = str2num(strtrim(vols));
volsPerTE = vols/numTEs;
tr = fpp.util.checkMRIProperty('tr',inputPath);
system(['fslsplit ' inputPath ' ' tmpDir '/split -t']);
for e=1:numTEs
    sepTEPaths{e} = [tmpDir '/rawdata' teSuffices{e} '.nii.gz'];
    mergeCmd{e} = ['fslmerge -tr ' sepTEPaths{e}];
end
for t=1:vols
    echoNum = mod(t,numTEs); if echoNum==0, echoNum=numTEs; end
    mergeCmd{echoNum} = [mergeCmd{echoNum} ' ' tmpDir '/split' fpp.util.numPad(t-1,4) '.nii.gz'];
end
for e=1:numTEs
    mergeCmd{e} = [mergeCmd{e} ' ' num2str(tr)];
    system(mergeCmd{e});
end


for e=1:numTEs
    meanFuncSep{e} = [tmpDir '/mean_func' teSuffices{e} '.nii.gz'];
    system(['fslmaths ' sepTEPaths{e} ' -Tmean ' meanFuncSep{e}]);
end
meanFuncCombined = [tmpDir '/mean_func_TEsep.nii.gz'];
mergeCmd{1} = ['fslmerge -tr ' meanFuncCombined];
for e=1:numTEs
    mergeCmd{1} = [mergeCmd{1} ' ' meanFuncSep{e}];
end
mergeCmd{1} = [mergeCmd{1} ' ' num2str(tr)];
system(mergeCmd{1});

funcData = fpp.util.mriRead(meanFuncCombined);
dims = size(funcData.vol);

logFuncMat = reshape(log(funcData.vol),[prod(dims(1:3)) numTEs]);
X = [ones(numTEs,1) teVals];
betas = inv(X'*X)*X'*logFuncMat';
r2Vec = -betas(2,:)';      % Estimated T2* (ms)
r2Vec(funcData.vol(:,:,:,1)==0) = 0;	% Apply simple brain mask
r2Vol = reshape(r2Vec,dims(1:3));
t2Vec = -betas(2,:)'.^-1;      % Estimated T2* (ms)
t2Vec(funcData.vol(:,:,:,1)==0) = 0;    % Apply simple brain mask
t2Vol = reshape(t2Vec,dims(1:3));

newData = funcData;
newData.vol = r2Vol;
fpp.util.mriWrite(newData,r2StarPath);
newData.vol = t2Vol;
fpp.util.mriWrite(newData,t2StarPath);

% Weight by estimated T2*
% For cases where estimated T2*<=0, use equal weighting of all three TEs.
for e=1:numTEs
    funcSepData{e} = fpp.util.mriRead(sepTEPaths{e});
end
weightVol = zeros([dims(1:3) numTEs]);
for e=1:numTEs
    weightVol(:,:,:,e) = teVals(e)*(t2Vol>0).*exp(-teVals(e)./t2Vol) + (t2Vol<=0);
end
weightVol = weightVol./repmat(sum(weightVol,4),[1 1 1 numTEs]);
weightedFuncVol = zeros([dims(1:3) volsPerTE]);
for e=1:numTEs
    weightedFuncVol = weightedFuncVol + funcSepData{e}.vol.*repmat(weightVol(:,:,:,e),[1 1 1 volsPerTE]);
end

newData.vol = weightedFuncVol;
fpp.util.mriWrite(newData,outputPath,'short');

system(['rm -rf ' tmpDir]);

end
