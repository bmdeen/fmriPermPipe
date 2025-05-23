
% Wrapper for FSL's fnirt, a nonlinear registration tool. Generates .json
% metadata files for output images and transform files.
%
% fpp.fsl.fnirt(inputPath,referencePath,varargin)
%
% Arguments:
% - inputPath (string): path to input data
% - referencePath (string): path to target image to warp to
%
% Variable arguments:
% - aff, inwarp, intin, cout, iout, fout, jout, logout, config, refmask,
%   inmask, minmet, miter, subsamp, warpres, infwhm, reffwhm, regmod,
%   jacrange, intmod, biasres, numprec, interp, splineorder, intorder,
%   imprefval, impinval, lambda, biaslambda, imprefm, impinm, ssqlambda,
%   refderiv, estint (see fnirt documentation for further info)
%
% Main output, cout: warp field coefficient image.
%
% Note: scalar numeric arguments should be specified as numeric. All other
% variable arguments should be specified as strings (e.g., '5,5,5,5')

function fnirt(inputPath,referencePath,varargin)

% Define variable defaults
jsonOpts.indent = '\t';     % Use tab indentation for JSON outputs

% Edit variable arguments.  Note: optInputs checks for proper input.
varArgList = {'aff','inwarp','intin','cout','iout','fout','jout','logout','config',...
    'refmask','inmask','minmet','miter','subsamp','warpres','infwhm','reffwhm',...
    'regmod','jacrange','intmod','biasres','numprec','interp','splineorder','intorder',...
    'imprefval','impinval','lambda','biaslambda','imprefm','impinm','ssqlambda',...
    'refderiv','estint'};
for i=1:length(varArgList)
    argVal = fpp.util.optInputs(varargin,varArgList{i});
    if ~isempty(argVal)
        if ischar(argVal)
            eval([varArgList{i} ' = argVal;']);
        else
            eval([varArgList{i} ' = num2str(argVal);']);
        end
    else
        eval([varArgList{i} ' = [];']);
    end
end

% Define fnirt command
cmd = ['fnirt --in=' inputPath ' --ref=' referencePath];

% Add additional variables
for i=1:length(varArgList)
    if ~isempty(eval(varArgList{i}))
        eval(['cmd = [cmd '' --' varArgList{i} '=' eval(varArgList{i}) '''];']);
    end
end

% Run fnirt command
fpp.util.system(cmd);

% Write json output files
% Warp coefficients (main xfm file, generated by default)
if isempty(cout)
    cout = strrep(inputPath,inputExt,'_warpcoef.nii.gz');
end
jsonData.Type = 'Nonlinear';
jsonData.Software = 'FNIRT';
jsonData.Invertible = true;
jsonData.FromFile = inputPath;
jsonData.ToFile = referencePath;
jsonData.CommandLine = cmd;
jsonData.Description = 'Warp coefficient file generated by FNIRT.';
bids.util.jsonencode(fpp.bids.jsonPath(cout),jsonData,jsonOpts);
% Warp field
if ~isempty(fout)
    bids.util.jsonencode(fpp.bids.jsonPath(fout),jsonData,jsonOpts);
end
% Warp jacobian
if ~isempty(jout)
    bids.util.jsonencode(fpp.bids.jsonPath(jout),jsonData,jsonOpts);
end
% Warped output image
if ~isempty(fpp.bids.getMetadata(inputPath)) && ~isempty(iout)
    outputPath = iout;
    fpp.bids.jsonReconstruct(inputPath,outputPath,'midprepfmri');
    refJsonData = fpp.bids.getMetadata(referencePath);
    if isfield(refJsonData,'SpatialRef')
        spatialRef = refJsonData.SpatialRef;
    else
        spatialRef = fpp.bids.removeBidsDir(referencePath);
    end
    fpp.bids.jsonChangeValue(outputPath,'SpatialRef',spatialRef);
    if ~strcmp(inputPath,outputPath)
        fpp.bids.jsonChangeValue(outputPath,'Sources',fpp.bids.removeBidsDir(inputPath));
    end
end

end