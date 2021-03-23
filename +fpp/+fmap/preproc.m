% 
% fpp.fmap.preproc(inputPaths,outputDir,varargin)
% 
% Preprocesses a pair of opposing phase-encode-direction spin-echo
% acquisitions, including copying from raw data to derivatives, and running
% topup to compute an undistortion warp field for functional data.
% 
% Example usage:
% fpp.func.preproc({'/pathToData/sub-01_dir-AP_run-01_epi.nii.gz',...
%       '/pathToData/sub-01_dir-PA_run-01_epi.nii.gz'},'/bidsDir/sub-01/');
% 
% Required arguments:
% - inputPaths (cell array of strings): paths to input spin-echo data,
%       one path for each of two opposing phase-encode directions. The
%       first must functional data in phase-encode direction.
% - outputDir (string): path to output directory (subject/session dir, in
%       BIDS filesystem)
% 
% Variable arguments:
% - overwrite (boolean; default=0): whether to overwrite files that have
%       already been written by this function.
% - phaseEncodeDirections (cell array of strings): phase-encode
%       direction of each spin echo map, in BIDS format (e.g. i, -j)
% - fieldMapParamPath (string): path to field map info file (topup input)
% 

function preproc(inputPaths,outputDir,varargin)

% Check system configuration
fpp.util.checkConfig;

jsonOpts.indent = '\t';     % Use tab indentation for JSON outputs

% Define wrapper function for fpp.bids.removeBidsDir (for cellfun functionality)
removeBidsDir = @(x) fpp.bids.removeBidsDir(x);

% Basic parameters
overwrite = 0;                      % Whether to overwrite output
phaseEncodeDirections = {}; % Spin echo map phase encode directions (overrides metadata info)

% Edit variable arguments.  Note: optInputs checks for proper input.
varArgList = {'overwrite','phaseEncodeDirections','fieldMapParamPath'};
for i=1:length(varArgList)
    argVal = fpp.util.optInputs(varargin,varArgList{i});
    if ~isempty(argVal)
        eval([varArgList{i} ' = argVal;']);
    end
end

% Check spin echo field map PE directions, if not specified
if ~exist('phaseEncodeDirections','var') || isempty(phaseEncodeDirections)
    for f=1:2
        tmp = fpp.util.checkMRIProperty('PEDir',inputPaths{f});
        if ~isempty(tmp)
            phaseEncodeDirections{f} = tmp;
        else
            error(['Could not determine phase-encode direction from JSON metadata;'...
                   'phaseEncodeDirections must be specified.']);
        end
    end
elseif size(phaseEncodeDirections,1)>size(phaseEncodeDirections,2)
    phaseEncodeDirections = phaseEncodeDirections'; % Ensure row array
end
% Convert PE directions to string output format (e.g. 'AP')
for f=1:2
    phaseEncodeDirectionsStr{f} = fpp.util.convertBidsPEDirToStr(inputPaths{f},phaseEncodeDirections{f});
end

% Define input name, output directories
[~,inputName,~] = fpp.util.fileParts(inputPaths{1});
inputNameGeneric = strrep(fpp.bids.changeName(inputName,'dir',[]),'_epi','');
if strcmp(outputDir(end),'/'), outputDir = outputDir(1:end-1); end
fmapPreprocDir = [outputDir '/fmap'];
if ~exist(fmapPreprocDir,'dir'), mkdir(fmapPreprocDir); end

% Topup output files
topupOutputStem = [fmapPreprocDir '/' inputNameGeneric];
topupWarpPath = [fpp.bids.changeName(topupOutputStem,{'from','to','mode','desc'},...
    {'native','Undistorted','image',['UndistortionWarp' phaseEncodeDirectionsStr{1}]}) '_xfm.nii.gz'];
topupJacobianPath = [fpp.bids.changeName(topupOutputStem,{'from','to','mode','desc'},...
    {'native','Undistorted','image',['UndistortionWarp' phaseEncodeDirectionsStr{1}]}) '_jacobian.nii.gz'];
topupFieldMapPath = [topupOutputStem '_fieldmapHz.nii.gz'];
topupFieldCoefPath = [topupOutputStem '_fieldcoef.nii.gz'];
topupUndistortedPath = [fpp.bids.changeName(topupOutputStem,'desc','Undistorted') '_epi.nii.gz'];

% Check if output exists.
if exist(topupUndistortedPath,'file') && ~overwrite
    return;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 1: Copy and concatenate images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 1, Copy and concatenate images            - ' inputNameGeneric]);
for f=1:2
    [~,spinEchoName,~] = fpp.util.fileParts(inputPaths{f});
    outputPaths{f} = [fmapPreprocDir '/' fpp.bids.changeName(spinEchoName,'dir',phaseEncodeDirectionsStr{f},'epi') '.nii.gz'];
    fpp.util.copyImageAndJson(inputPaths{f},outputPaths{f},'fmap');
    fpp.bids.jsonChangeValue(outputPaths{f},'RawSources',inputPaths{f});
end
rawPaths = inputPaths;
inputPaths = outputPaths;
inputPaths{3} = fpp.bids.changeName(inputPaths{1},'dir','Both');
fpp.util.system(['fslmerge -t ' inputPaths{3} ' ' inputPaths{1} ' ' inputPaths{2}]);
fpp.bids.jsonReconstruct(inputPaths{1},inputPaths{3},'fmap');
fpp.bids.jsonChangeValue(inputPaths{3},{'PhaseEncodingDirection','RawSources'},...
    {[],cellfun(removeBidsDir,rawPaths,'UniformOutput',false)});
topupLogPath = [strrep(inputPaths{3},'.nii.gz','') '.topup_log'];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 2: Write field map parameter file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 2, Write field map parameter file         - ' inputNameGeneric]);
% Check topup properties, write to field map param file
if ~exist('fieldMapParamPath','var') || isempty(fieldMapParamPath) || ~exist(fieldMapParamPath,'file')
    fmapProperties = [];
    for f=1:2
        fmapProperties(end+1,:) = fpp.util.checkMRIProperty('Topup',inputPaths{f});
    end
    if isempty(fmapProperties)
        error('Could not determine spin echo topup parameters; fieldMapParamPath must be specified.');
    end
    fieldMapParamPath = strrep(inputPaths{3},'_epi.nii.gz','_topupparams.txt');
    fid = fopen(fieldMapParamPath,'w');
    fprintf(fid,'%d %d %d %f\n',fmapProperties');
    fclose(fid);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STEP 3: Run topup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',['Step 3, Compute undistortion warp (topup)      - ' inputNameGeneric]);
% Determine configuration file based on field map dimensions
dims = fpp.util.checkMRIProperty('dims',inputPaths{1});
if sum(mod(dims,4))==0
    configName = 'b02b0_4.cnf';
elseif sum(mod(dims,2))==0
    configName = 'b02b0_2.cnf';
else
    configName = 'b02b0_1.cnf';
end
cmd = ['topup --imain=' inputPaths{end} ' --datain=' fieldMapParamPath ' --config=' configName ' --out=' ...
    topupOutputStem ' --iout=' fpp.bids.changeName(topupOutputStem,'desc','Undistorted') '_epi --fout=' ...
    topupOutputStem '_fieldmapHz --dfout=' topupOutputStem '_warp --jacout=' topupOutputStem '_jacobian'];
fpp.util.system(cmd);
fpp.util.system(['mv ' topupOutputStem '_jacobian_01.nii.gz ' topupJacobianPath]);
fpp.util.system(['rm -rf ' topupOutputStem '_jacobian_02.nii.gz']);
fpp.util.system(['mv ' topupOutputStem '_warp_01.nii.gz ' topupWarpPath]);
fpp.util.system(['rm -rf ' topupOutputStem '_warp_02.nii.gz']);
% Define json info
jsonData.Sources = cellfun(removeBidsDir,inputPaths(1:2),'UniformOutput',false);
inputJsonData = fpp.bids.getMetadata(inputPaths{1});
if ~isempty(inputJsonData) && isfield(inputJsonData,'RawSources')
    jsonData.RawSources = inputJsonData.RawSources;
end
% For fieldcoef
jsonData.Description = 'Topup-derived field map spline coefficients, to be used with applytopup.';
bids.util.jsonencode(fpp.bids.jsonPath(topupFieldCoefPath),jsonData,jsonOpts);
% For fieldmap
jsonData.SpatialRef = removeBidsDir(topupUndistortedPath);
jsonData.Description = 'Topup-derived field map (Hz).';
bids.util.jsonencode(fpp.bids.jsonPath(topupFieldMapPath),jsonData,jsonOpts);
% For undistorted
jsonData.Description = 'Undistored spin-echo images, generated by topup.';
bids.util.jsonencode(fpp.bids.jsonPath(topupUndistortedPath),jsonData,jsonOpts);
% For warp
rmfield(jsonData,'SpatialRef');
jsonData.Type = 'Nonlinear';
jsonData.Software = 'topup';
jsonData.Invertible = true;
jsonData.FromFile = inputPaths{1};
jsonData.ToFile = topupUndistortedPath;
jsonData.CommandLine = cmd;
jsonData.Description = 'Warp coefficient file for undistortion, generated by topup.';
bids.util.jsonencode(fpp.bids.jsonPath(topupWarpPath),jsonData,jsonOpts);
% For jacobian
jsonData.Description = 'Jacobian file for undistortion warp, generated by topup.';
bids.util.jsonencode(fpp.bids.jsonPath(topupJacobianPath),jsonData,jsonOpts);

end