
% Wrapper for FSL's invwarp, to invert a nonlinear warp image. Generates
% .json metadata for output warp coefficient image.
%
% fpp.fsl.invertWarp(inputWarp,outputWarp,referencePath,varargin)
%
% Arguments:
% - inputWarp (string): path to input warp coefficient file
% - outputWarp (string): path to output warp coefficient file
% - referencePath (string): path to the new reference image, i.e. the input
%       of the inputWarp transform
%
% Optional arguments:
% - jmin, jmax, rel, abs, noconstraint (see invwarp documentation for
%       further info)

function invertWarp(inputWarp,outputWarp,referencePath,varargin)

% Edit variable arguments.  Note: optInputs checks for proper input.
varArgList = {'jmin','jmax','rel','abs','noconstraint'};
indArgIsBoolean = 3:5;
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

cmd = ['invwarp --warp=' inputWarp ' --out=' outputWarp ' --ref=' referencePath];

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

% Run convert_xfm command
fpp.util.system(cmd);

% Write json output files
if ~isempty(fpp.bids.getMetadata(inputWarp))
    fpp.bids.jsonReconstruct(inputWarp,outputWarp,'xfm');
    jsonData = fpp.bids.getMetadata(fpp.bids.jsonPath(inputWarp));
    fpp.bids.jsonChangeValue(outputWarp,{'FromFile','ToFile','CommandLine'},...
        {jsonData.ToFile,jsonData.FromFile,cmd});
end

end