
% Function to separate multi-echo data into individual echo series, using
% BIDS naming convention: "_echo-N_"

function outputPaths = multiEchoSeparate(inputPath,nEchoes,tr)

[outputDir,inputName,inputExt] = filepartsGZ(inputPath);
randStem = 'split123140975165-';

for e=1:nEchoes
    outputPaths{e} = bidsChangeEntity(inputPath,'echo',int2str(e));
    mergeCmd{e} = ['fslmerge -tr ' outputPaths{e}];
end

[~,vols] = system(['fslval ' inputPath ' dim4']);
vols = str2num(strtrim(vols));
system(['fslsplit ' inputPath ' ' outputDir '/' randStem ' -t']);

for t=1:vols
    tReal = ceil(t/nEchoes);
    echoNum = mod(t,nEchoes); if echoNum==0, echoNum=nEchoes; end
    mergeCmd{echoNum} = [mergeCmd{echoNum} ' ' outputDir '/' randStem numPad(t-1,4) '.nii.gz'];
end
for e=1:nEchoes
    mergeCmd{e} = [mergeCmd{e} ' ' num2str(tr)];
    system(mergeCmd{e});
end
system(['rm -rf ' outputDir '/' randStem '*']);

end