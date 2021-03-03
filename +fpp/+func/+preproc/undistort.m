
% Function to undistort a functional volume using blip-up/blip-down spin
% echo measurements.
%
% [xfmFunc2SpinEcho,xfmSpinEcho2Func] = fpp.func.preproc.undistort(...
%   inputPath,outputPath,topupWarpPath,topupJacobianPath)
%
% Arguments:
% - inputPath (string): path to input distorted functional
% - outputPath (string): path to output undistorted functional
% - spinEchoPath (string): path to spin-echo image with phase-encode
%   direction matched to inputPath
% - topupWarpPath (string): path to topup-based undistortion warp
% - topupJacobianPath (string): path to topup-based warp jacobian
%

function [xfmFunc2SpinEcho,xfmSpinEcho2Func] = undistort(inputPath,outputPath,...
    spinEchoPath,topupWarpPath,topupJacobianPath)

topupJacobian2FuncPath = fpp.bids.changeName(inputPath,{'desc','echo'},{'UndistortionWarp',[]},'jacobian','.nii.gz');

% Register input func volume to spin echo image
xfmFunc2SpinEcho = fpp.bids.changeName(inputPath,{'desc','from','to','mode','echo'},...
    {'','native','SpinEcho','image',[]},'xfm','.mat');
xfmSpinEcho2Func = fpp.bids.changeName(inputPath,{'desc','from','to','mode','echo'},...
    {'','SpinEcho','native','image',[]},'xfm','.mat');
fpp.fsl.flirt(inputPath,spinEchoPath,xfmFunc2SpinEcho,[],'cost','corratio',...
    'dof',6,'searchrx',[-90 90],'searchry',[-90 90],'searchrz',[-90 90]);
fpp.fsl.invertXfm(xfmFunc2SpinEcho,xfmSpinEcho2Func);

% Undistort, by applying warp and multiplying by Jacobian.
fpp.fsl.moveImage(inputPath,inputPath,outputPath,xfmFunc2SpinEcho,...
    'warp',topupWarpPath,'postmat',xfmSpinEcho2Func);
fpp.fsl.moveImage(topupJacobianPath,inputPath,topupJacobian2FuncPath,xfmSpinEcho2Func);
fpp.fsl.maths(outputPath,['-mul ' topupJacobian2FuncPath],outputPath);
fpp.util.system(['rm -rf ' topupJacobian2FuncPath]);

% Generate output JSON file
fpp.bids.jsonReconstruct(inputPath,outputPath,'midprepfmri');

end