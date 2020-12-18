
% Wrapper for FSL's fslmaths tool, generates JSON metadata for output
%
% fpp.fsl.fslMaths(inputPath,cmd,outputPath,newDesc,outputDescription,appendDescription)
%
% Arguments:
% - inputPath: path to input image
% - cmd: fslmaths command text between input/output paths
% - outputPath: path to output image
% - outputDescription (optional): new value for JSON Description field
% - appendDescription (boolean, optional): whether to append
%       outputDescription to existing description
% - odt (string): data type (default: float)

function fslMaths(inputPath,cmd,outputPath,outputDescription,appendDescription,odt)

if ~exist('outputDescription','var') || isempty(outputDescription)
    outputDescription = '';
end
if ~exist('appendDescription','var') || isempty(appendDescription)
    appendDescription = 0;
end

[~,~,inputExt] = fpp.util.fileParts(inputPath);
[~,~,outputExt] = fpp.util.fileParts(outputPath);

cmdFull = ['fslmaths ' inputPath ' ' cmd ' ' outputPath];
if exist('odt','var') && ~isempty(odt)
    cmdFull = [cmdFull ' -odt ' odt];
end
fpp.util.system(cmdFull);

inputJsonPath = strrep(inputPath,inputExt,'.json');
if exist(inputJsonPath,'file')
    outputJsonPath = strrep(outputPath,outputExt,'.json');
    fpp.util.system(['cp ' inputJsonPath ' ' outputJsonPath]);
    if ~isempty(outputDescription)
        fpp.bids.jsonChangeValue(outputJsonPath,'Description',outputDescription,appendDescription);
    end
end

end