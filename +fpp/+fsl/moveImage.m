
% Wrapper for FSL's flirt and applywarp, used to transform an MRI image
% from one space to another. Writes a json file for the output.
%
% If variable argument warp is specified, uses applywarp. Otherwise, uses
% flirt -applyxfm.

function moveImage(inputPath,referencePath,outputPath,preMat,varargin)

% Define variable defaults
jsonOpts.indent = '\t';     % Use tab indentation for JSON outputs
warp = [];                  % Prevent confusion with MATLAB warp fcn

% Edit variable arguments.  Note: optInputs checks for proper input.
varArgList = {'warp','datatype','postmat','mask','interp','superlevel','paddingsize',...
                'abs','rel','super','usesqform'};
indArgIsBoolean = 8:11;
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

if ~isempty(warp)
    % Build applywarp command
    
    % Ensure interp is formatted correctly
    if ~isempty(interp) && (strcmp(interp,'nearestneighbour') || strcmp(interp,'nearestneighbor'))
        interp = 'nn';
    end
    
    % Define applywarp command
    cmd = ['applywarp --in=' inputPath ' --ref=' referencePath ' --out=' outputPath];
    
    if ~isempty('preMat')
        cmd = [cmd ' --premat=' preMat];
    end

    % Add non-boolean additional variables
    for i=setdiff(1:length(varArgList),indArgIsBoolean)
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
    for i=setdiff(1:length(varArgList),indArgIsBoolean)
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

% Run fnirt command
fpp.util.system(cmd);

% Write json output files
if exist(fpp.bids.jsonPath(inputPath),'file') && exist('outputPath','var')
    fpp.bids.jsonReconstruct(inputPath,outputPath);
    fpp.bids.jsonChangeValue(outputPath,'SpatialRef',...
        strrep(referencePath,fpp.bids.checkBidsDir(referencePath),''));
end


end