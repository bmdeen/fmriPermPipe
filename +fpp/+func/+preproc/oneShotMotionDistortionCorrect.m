
% Function to apply motion and distortion correct, and registration to
% template space with a single interpolation step.
%
% fpp.func.preproc.oneShotMotionDistortionCorrect(inputPaths,outputPaths,templatePath,...
%   templateSpace,mcDir,topupWarpPath,topupJacobianPath,xfmNativeFunc2Template,...
%   xfmSpinEcho2NativeFunc,xfmNativeFunc2SpinEcho,echoForMoCorr)
%
% Assumes matched naming between input and mcDir.
%

function oneShotMotionDistortionCorrect(inputPaths,outputPaths,templatePath,templateSpace,mcDir,topupWarpPath,...
    topupJacobianPath,xfmNativeFunc2Template,xfmSpinEcho2NativeFunc,xfmNativeFunc2SpinEcho,echoForMoCorr)

xfmSpinEcho2Template = fpp.bids.changeName(topupJacobianPath,{'desc','from','to','mode'},...
    {'','native',templateSpace,'image'},'xfm','.mat');
topupJacobian2TemplatePath = fpp.bids.changeName(topupJacobianPath,'space',templateSpace,'jacobian');
[~,mcName,~] = fileparts(mcDir);

% Move warp Jacobian to template space
if ~exist(xfmSpinEcho2Template,'file')
    fpp.fsl.concatXfm(xfmNativeFunc2Template,xfmSpinEcho2NativeFunc,xfmSpinEcho2Template);
end
if ~exist(topupJacobian2TemplatePath,'file')
    fpp.fsl.moveImage(topupJacobianPath,templatePath,topupJacobian2TemplatePath,xfmSpinEcho2Template);
end

% Copy warp/postmat to avoid applywarp segmentation faults when running multiple runs in parallel
[~,xfmSE2TempName,~] = fileparts(xfmSpinEcho2Template);
xfmSpinEcho2TemplateCopy = [mcDir '/' xfmSE2TempName '.mat'];
fpp.util.system(['cp ' xfmSpinEcho2Template ' ' xfmSpinEcho2TemplateCopy]);
[~,topupWarpName,topupWarpExt] = fpp.util.fileParts(topupWarpPath);
topupWarpPathCopy = [mcDir '/' topupWarpName topupWarpExt];
fpp.util.system(['cp ' topupWarpPath ' ' topupWarpPathCopy]);

% Check # of volumes
vols = fpp.util.checkMRIProperty('vols',inputPaths{1});
tr = fpp.util.checkMRIProperty('tr',inputPaths{1});

for e=1:length(outputPaths)
    mergeCmd = ['fslmerge -tr ' outputPaths{e}];
    inputSplitStem = [strrep(inputPaths{e},'.nii.gz','') '_SplitForMoco'];
    outputSplitStem = [strrep(fpp.bids.changeName(inputPaths{e},'space',templateSpace),'.nii.gz','') '_SplitForMoco'];
    fpp.util.system(['fslsplit ' inputPaths{e} ' ' inputSplitStem ' -t']);
    for t=0:vols-1
        inputVolPath = [inputSplitStem fpp.util.numPad(t,4) '.nii.gz'];
        outputVolPath = [outputSplitStem fpp.util.numPad(t,4) '.nii.gz'];
        xfmInputVol2NativeFunc = [mcDir '/' fpp.bids.changeName(mcName,{'echo','desc','from','to','mode'},...
            {int2str(echoForMoCorr),'',['native' fpp.util.numPad(t,4)],'native','image'},'xfm','.mat')];
        xfmInputVol2SpinEcho = [mcDir '/' fpp.bids.changeName(mcName,{'echo','desc','from','to','mode'},...
            {int2str(echoForMoCorr),'',['native' fpp.util.numPad(t,4)],'SpinEcho','image'},'xfm','.mat')];
        fpp.fsl.concatXfm(xfmNativeFunc2SpinEcho,xfmInputVol2NativeFunc,xfmInputVol2SpinEcho);
        for i=1:10      % Bug fix, catches random segmentation faults by repeating the applywarp command up to 10 times
            try
                fpp.fsl.moveImage(inputVolPath,templatePath,outputVolPath,xfmInputVol2SpinEcho,...
                    'warp',topupWarpPathCopy,'postmat',xfmSpinEcho2TemplateCopy);
                break;
            catch exception
                if i<10
                    warning(exception.message);
                    warning('applywarp errored out, trying again');
                else
                    error('applywarp errored out multiple times.');
                end
            end
        end
        fpp.util.system(['fslmaths ' outputVolPath ' -mul ' topupJacobian2TemplatePath ' ' outputVolPath]);
        mergeCmd = [mergeCmd ' ' outputVolPath];
        if mod(t+1,10)==0
           fprintf('\t%s\n',['Finished warping echo ' int2str(e) ', volume ' int2str(t+1)]);
        end
    end
    
    % Merge output into single time series file, delete individual time points
    mergeCmd = [mergeCmd ' ' num2str(tr)];
    fpp.util.system(mergeCmd);
    fpp.util.system(['rm -rf ' inputSplitStem '*.* ' outputSplitStem '*.*']);
    
    % Generate output JSON file
    fpp.bids.jsonReconstruct(inputPaths{e},outputPaths{e},'midprepfmri');
    fpp.bids.jsonChangeValue(outputPaths{e},{'Sources','SpatialRef'},...
        {fpp.bids.removeBidsDir(inputPaths{e}),fpp.bids.removeBidsDir(templatePath)});
end

fpp.util.system(['rm -rf ' xfmSpinEcho2TemplateCopy ' ' topupWarpPathCopy]);

end