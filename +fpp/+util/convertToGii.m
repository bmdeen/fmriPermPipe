
% 
% fpp.util.convertToGii(inputPath,outputPath,structure[,xfm,fieldsToChange,newValues])
%
% Script to generate GIFTI file (e.g., from Freesurfer surface file),
% including file conversion, structure specification, and optional
% application of an affine transform, and .json metadata file generation
% 
% Arguments:
% - inputPath (string): path to input surface file
%     OR (cell array, optional): inputPath{1} = shape or annot file,
%     inputPath{2} = associated surface file
% - outputPath (string): path to output GIFTI (.gii) file
% - structure (string): structure to specificy (e.g. CORTEX_LEFT)
%     OR (cell array, optional): structure{1} = structure, structure{2} =
%     surface type, structure{3} = surface secondary type
% - xfm (string, optional): path to affine xfm (wb_command "world" format) to apply
%     OR (cell array, optional): xfm{1} = xfm path, xfm{2} = flirt input
%     path, xfm{3} = flirt output path for FLIRT-based xfm
% - fieldsToChange (cell array of strings): JSON fields to edit
% - newValues (cell array): new values for those fields

function convertToGii(inputPath,outputPath,structure,xfm,fieldsToChange,newValues)

% Settings
jsonOpts.indent = '\t';     % Use tab indentation for JSON outputs
inputIsSurface = 1;         % Whether input is surface or shape/annot

% Check output extension
[~,~,outputExt] = fileparts(outputPath);
if ~strcmpi(outputExt,'.gii'), error('Output extention must be .gii'); end

% Convert to gifti
if iscell(inputPath)
    inputIsSurface = 0;
    [~,~,inputExt] = fpp.util.fileParts(inputPath{1});
    if strcmpi(inputExt,'.annot')
        flagStr = ['--annot ' inputPath{1}];
    elseif sum(regexpi(inputPath{1},'(sulc|curv|thickness)'))>0
        flagStr = ['-c ' inputPath{1}];
    else
        error('When specifying a cell array for inputPath, first argument must be .annot or sulc, curv, or thickness file.');
    end
    fpp.fs.mrisConvert(inputPath{2},outputPath,flagStr);
    inputPath = inputPath{1};
elseif ~strcmp(inputPath,outputPath)
    fpp.fs.mrisConvert(inputPath,outputPath);
end

% Add structure metadata
if iscell(structure)
    flagStr = '';
    if length(structure)>=2 && ~isempty(structure{2})
        flagStr = [flagStr ' -surface-type ' structure{2}];
    end
    if length(structure)>=3 && ~isempty(structure{3})
        flagStr = [flagStr ' -surface-type-secondary ' structure{3}];
    end
    fpp.wb.command('set-structure',outputPath,structure{1},[],flagStr);
else
    fpp.wb.command('set-structure',outputPath,structure);
end

% Apply affine xfm, if specified
if exist('xfm','var') && ~isempty(xfm) && inputIsSurface
    flirtStr = '';
    if iscell(xfm) && length(xfm)==3    % Include flirt source/target volume arguments
        flirtStr = [' -flirt ' xfm{2} ' ' xfm{3}];
    else
        xfm = {xfm};
    end
    fpp.wb.command('surface-apply-affine',outputPath,xfm{1},outputPath,flirtStr);
end

% Generate .json metadata, if specified
if exist('fieldsToChange','var') && ~isempty(fieldsToChange) && ...
        exist('newValues','var') && ~isempty(newValues)
    if ~isempty(fieldnames(fpp.bids.getMetadata(inputPath)))
        fpp.bids.jsonReconstruct(inputPath,outputPath);
        fpp.bids.jsonChangeValue(outputPath,fieldsToChange,newValues);
    else
        jsonData = struct();
        outputJsonPath = fpp.bids.jsonPath(outputPath);
        if ~iscell(fieldsToChange), fieldsToChange = {fieldsToChange}; end
        if ~iscell(newValues), newValues = {newValues}; end
        for f=1:length(fieldsToChange)
            jsonData.(fieldsToChange{f}) = newValues{f};
        end
        bids.util.jsonencode(outputJsonPath,jsonData,jsonOpts);
    end
end

end