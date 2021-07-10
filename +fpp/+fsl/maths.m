
% Wrapper for FSL's fslmaths tool, generates JSON metadata for output
%
% fpp.fsl.maths(inputPath,cmd,outputPath,newDesc,outputDescription,appendDescription)
%
% Arguments:
% - inputPath: path to input image
% - cmd: fslmaths command text between input/output paths
% - outputPath: path to output image
%
% Optional arguments:
% - outputDescription: new value for JSON Description field
% - appendDescription (boolean): whether to append outputDescription
%   to existing description
% - odt (string): data type (default: float)

function maths(inputPath,cmd,outputPath,outputDescription,appendDescription,odt)

if ~exist('outputDescription','var') || isempty(outputDescription)
    outputDescription = '';
end
if ~exist('appendDescription','var') || isempty(appendDescription)
    appendDescription = 0;
end

cmdFull = ['fslmaths ' inputPath ' ' cmd ' ' outputPath];
if exist('odt','var') && ~isempty(odt)
    cmdFull = [cmdFull ' -odt ' odt];
end
fpp.util.system(cmdFull);

if ~isempty(fpp.bids.getMetadata(inputPath))
    if ~strcmp(inputPath,outputPath)
        fpp.bids.jsonReconstruct(inputPath,outputPath);
    end
    if ~isempty(outputDescription)
        fpp.bids.jsonChangeValue(outputPath,'Description',outputDescription,appendDescription);
    end
    if ~strcmp(inputPath,outputPath)
        fpp.bids.jsonChangeValue(outputPath,'Sources',fpp.bids.removeBidsDir(inputPath));
    end
end

end