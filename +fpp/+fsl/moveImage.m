
% Wrapper for FSL's flirt and applywarp, used to transform an MRI image
% from one space to another. Writes a json file for the output.
%
% If variable argument warp is specified, uses applywarp. Otherwise, uses
% flirt -applyxfm.
%
% fpp.fsl.moveImage(inputPath,referencePath,outputPath,preMat,varargin)
%
% Arguments:
% - inputPath (string): path to input image
% - referencPath (string): path to reference image (registration target)
% - outputPath (string): path to output image
%
% Optional arguments:
% - preMat (string): linear transformation to apply (necessary for linear
%   transformation, optional for nonlinear)
% 
% Variable arguments:
% - warp (string): path to warp coefficient image (cout from fnirt)
% - postMat (string): path to linear transform to apply after warp
% - isLabel (boolean): whether to treat image like workbench label file,
%       and copy label table
% - interp (string): interpolation type (nn/nearestneighbour, trilinear, 
%       sinc, spline)
% - datatype, mask, superlevel, paddingsize, abs, rel ,superlevel, 
%       usesqform - see flirt/applywarp documentation for info on these

function moveImage(inputPath,referencePath,outputPath,preMat,varargin)

% Define variable defaults. These two must be initialized to prevent
% confusion with MATLAB functions.
warp = [];
interp = [];

% Edit variable arguments.  Note: optInputs checks for proper input.
varArgList = {'isLabel','warp','datatype','postmat','mask','interp','superlevel','paddingsize',...
                'abs','rel','super','usesqform'};
indArgIsBoolean = 9:12;
for i=1:length(varArgList)
    argVal = fpp.util.optInputs(varargin,varArgList{i});
    if ~isempty(argVal)
        if ischar(argVal) || strcmp(varArgList{i},'isLabel')
            eval([varArgList{i} ' = argVal;']);
        else
            eval([varArgList{i} ' = num2str(argVal);']);
        end
    else
        eval([varArgList{i} ' = [];']);
    end
end

if isempty(isLabel), isLabel = 0; end

if ~isempty(warp)
    % Build applywarp command
    
    % Ensure interp is formatted correctly
    if ~isempty(interp) && (strcmp(interp,'nearestneighbour') || strcmp(interp,'nearestneighbor'))
        interp = 'nn';
    end
    
    % Define applywarp command
    cmd = ['applywarp --in=' inputPath ' --ref=' referencePath ' --out=' outputPath];
    
    if ~isempty(preMat)
        cmd = [cmd ' --premat=' preMat];
    end

    % Add non-boolean additional variables
    for i=setdiff(2:length(varArgList),indArgIsBoolean)
        if ~isempty(eval(varArgList{i}))
            eval(['cmd = [cmd '' --' varArgList{i} '=' eval(varArgList{i}) '''];']);
        end
    end
    % Add boolean additional variables
    for i=indArgIsBoolean
        if ~isempty(eval(varArgList{i}))
            eval(['cmd = [cmd '' --' varArgList{i} '''];']);
        end
    end
else
    % Build flirt command
    
    % Ensure interp is formatted correctly
    if ~isempty(interp) && (strcmp(interp,'nn') || strcmp(interp,'nearestneighbor'))
        interp = 'nearestneighbour';
    end
    
    % Define flirt command
    cmd = ['flirt -in ' inputPath ' -ref ' referencePath ' -out ' outputPath...
        ' -applyxfm -init ' preMat];
    
    % Add non-boolean additional variables
    for i=setdiff(2:length(varArgList),indArgIsBoolean)
        if ~isempty(eval(varArgList{i}))
            eval(['cmd = [cmd '' -' varArgList{i} ' ' eval(varArgList{i}) '''];']);
        end
    end
    % Add boolean additional variables
    for i=indArgIsBoolean
        if ~isempty(eval(varArgList{i}))
            eval(['cmd = [cmd '' -' varArgList{i} '''];']);
        end
    end
end

% Run flirt or applywarp command
fpp.util.system(cmd);

% Copy label table, for Workbench label files
if isLabel
    tmpLUT = fpp.bids.changeName(outputPath,'desc','tmpMoveImage209851343149','lut','.txt');
    fpp.wb.command('volume-label-export-table',inputPath,'1',tmpLUT);
    fpp.wb.command('volume-label-import',outputPath,tmpLUT,outputPath);
    fpp.util.system(['rm -rf ' tmpLUT]);
end

% Write json output files
if ~isempty(fpp.bids.getMetadata(inputPath)) && exist('outputPath','var')
    fpp.bids.jsonReconstruct(inputPath,outputPath);
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