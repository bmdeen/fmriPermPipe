
% Function to undistort a functional volume using blip-up/blip-down spin
% echo measurements.
%
% Arguments:
% inputPath (string): path to input distorted functional
% outputPath (string): path to output undistorted functional
% spinEchoPaths (cell array of strings): paths spin echo maps with opposite
%   phase encode directions (first must be matched in PE dir to functional)
% fmapOutputDir (string): directory in which to output field map-related
%   derivatives
% fieldMapParamPath: text file with field map parameters for topup
%
%
% TODO:
% - Topup: accomodate odd slice number by adding dummy slice in this case
% - Generate json files for produced images, with Description, Sources,
%   SkullStripped, and SpatialReference fields.

function [errorMsg,topupWarpPath,topupJacobianPath,xfmFunc2SpinEcho,xfmSpinEcho2Func] = ...
    undistort(inputFuncPath,inputFuncUndistortedPath,spinEchoPaths,fmapPreprocDir,fieldMapParamPath)

errorMsg = [];

% Topup output files
[~,spinEchoName,~] = fpp.util.fileParts(spinEchoPaths{1});
topupOutputStem = [fmapPreprocDir '/' strrep(fpp.bids.changeName(spinEchoName,'dir',''),'_epi','')];
topupWarpPath = [fpp.bids.changeName(topupOutputStem,'desc','UndistortionWarp') '_xfm.nii.gz'];
topupJacobianPath = [fpp.bids.changeName(topupOutputStem,'desc','UndistortionWarp') '_jacobian.nii.gz'];
topupJacobian2FuncPath = fpp.bids.changeName(inputFuncPath,{'desc','echo'},{'UndistortionWarp',[]},'jacobian','.nii.gz');

% Check topup properties, write to field map param file
if ~exist('fieldMapParamPath','var') || isempty(fieldMapParamPath) || ~exist(fieldMapParamPath,'file')
    fmapProperties = [];
    for f=1:length(spinEchoPaths)
        fmapProperties(end+1,:) = fpp.util.checkMRIProperty('Topup',spinEchoPaths{f});
    end
    if isempty(fmapProperties)
        errorMsg = 'ERROR: Could not determine spin echo topup parameters.';
        return;
    end
    fieldMapParamPath = fpp.bids.changeName(spinEchoPaths{3},[],[],'topupparams','.txt');
    fid = fopen(fieldMapParamPath,'w');
    fprintf(fid,'%d %d %d %f\n',fmapProperties');
    fclose(fid);
end

% Run topup. NOTE: assuming that first spin echo file is matched in
% phase encode direction to functional data.
%%% TO ADD: ACCOMODATE ODD-SLICE-NUMBER BY ADDING DUMMY SLICE
%%% Generate json files for warp/jacobian
if ~exist(topupWarpPath,'file') || ~exist(topupJacobianPath,'file')
    fpp.util.system(['topup --imain=' spinEchoPaths{end} ' --datain=' fieldMapParamPath ' --config=b02b0.cnf --out=' ...
        topupOutputStem ' --iout=' fpp.bids.changeName(topupOutputStem,'desc','Undistorted') '_epi --fout=' ...
        topupOutputStem '_fieldmapHz --dfout=' topupOutputStem '_warp --jacout=' topupOutputStem '_jacobian']);
    fpp.util.system(['mv ' topupOutputStem '_jacobian_01.nii.gz ' topupJacobianPath]);
    fpp.util.system(['rm -rf ' topupOutputStem '_jacobian_02.nii.gz']);
    fpp.util.system(['mv ' topupOutputStem '_warp_01.nii.gz ' topupWarpPath]);
    fpp.util.system(['rm -rf ' topupOutputStem '_warp_02.nii.gz']);
end

% Register input func volume to spin echo image
xfmFunc2SpinEcho = fpp.bids.changeName(inputFuncPath,{'desc','from','to','mode'},...
    {'','orig','SpinEcho','image'},'xfm','.mat');
xfmSpinEcho2Func = fpp.bids.changeName(inputFuncPath,{'desc','from','to','mode'},...
    {'','SpinEcho','orig','image'},'xfm','.mat');
fpp.util.system(['flirt -in ' inputFuncPath ' -ref ' spinEchoPaths{1} ' -omat ' xfmFunc2SpinEcho ...
    ' -cost corratio -dof 6 -searchrx -90 90 -searchry -90 90 -searchrz -90 90']);
fpp.util.system(['convert_xfm -omat ' xfmSpinEcho2Func ' -inverse ' xfmFunc2SpinEcho]);

% Undistort, by applying warp and multiplying by Jacobian.
fpp.util.system(['applywarp --in=' inputFuncPath ' --ref=' inputFuncPath ' --out=' inputFuncUndistortedPath ...
    ' --warp=' topupWarpPath ' --premat=' xfmFunc2SpinEcho ' --postmat=' xfmSpinEcho2Func]);
fpp.util.system(['flirt -in ' topupJacobianPath ' -ref ' inputFuncPath ' -out ' topupJacobian2FuncPath ...
    ' -applyxfm -init ' xfmSpinEcho2Func]);
fpp.util.system(['fslmaths ' inputFuncUndistortedPath ' -mul ' topupJacobian2FuncPath ' ' inputFuncUndistortedPath]);
fpp.util.system(['rm -rf ' topupJacobian2FuncPath]);

% Generate output JSON file
if exist(fpp.bids.jsonPath(inputFuncPath))
    fpp.bids.jsonReconstruct(inputFuncPath,inputUndistortedFuncPath);
end

end