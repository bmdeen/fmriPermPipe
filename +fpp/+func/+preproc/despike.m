
% Wrapper for AFNI's 3dDespike script. Converts output to .nii.gz, and
% reconstructs json metadata.
%
% fpp.func.preproc.despike(inputPath,outputPath)

function despike(inputPath,outputPath)

[inputDir,inputName,~] = fpp.util.fileParts(inputPath);
tmpDir = [inputDir '/' inputName '_despike10398513402'];
curDir = pwd;
mkdir(tmpDir);  cd(tmpDir);

fpp.util.system(['3dDespike -nomask ' inputPath]);
fpp.util.system(['3dAFNItoNIFTI ' pwd '/despike+orig.BRIK']);
fpp.fs.mriConvert([pwd '/despike.nii'],outputPath);
if ~strcmp(inputPath,outputPath)
    fpp.bids.jsonReconstruct(inputPath,outputPath,'midprepfmri');
end

cd(curDir);
fpp.util.system(['rm -rf ' tmpDir]);

end