
% Function to separate multi-echo data stored into individual echo series.
% - Assumes input has echoes concatenated across time:
%       T1E1, T1E2, ..., T1EN, T2E1, T2E2, ...
% - Uses BIDS naming convention: "_echo-N_"

function outputPaths = multiEchoSeparate(inputPath,nEchoes,tr)

[outputDir,~,~] = fileparts(inputPath);
randStem = 'split123140975165-';

for e=1:nEchoes
    outputPaths{e} = fpp.bids.changeName(inputPath,'echo',int2str(e));
    mergeCmd{e} = ['fslmerge -tr ' outputPaths{e}];
end

[~,vols] = fpp.util.system(['fslval ' inputPath ' dim4']);
vols = str2num(strtrim(vols));
fpp.util.system(['fslsplit ' inputPath ' ' outputDir '/' randStem ' -t']);

for t=1:vols
    tReal = ceil(t/nEchoes);
    echoNum = mod(t,nEchoes); if echoNum==0, echoNum=nEchoes; end
    mergeCmd{echoNum} = [mergeCmd{echoNum} ' ' outputDir '/' randStem fpp.util.numPad(t-1,4) '.nii.gz'];
end
for e=1:nEchoes
    mergeCmd{e} = [mergeCmd{e} ' ' num2str(tr)];
    fpp.util.system(mergeCmd{e});
end
fpp.util.system(['rm -rf ' outputDir '/' randStem '*']);

end