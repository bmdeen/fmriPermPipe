
% Function to change values in specific json metadata fields for an input 
% data file. Can also remove fields, by setting the corresponding element 
% of newValues to an empty vector.
%
% fpp.bids.jsonChangeValue(inputPath,fieldsToChange,newValues,appendValue)
%
% Arguments:
% - inputPath (string): path to file to modify (data file, not .json file)
% - fieldsToChange (cell array of strings): fields to edit
% - newValues (cell array): new values for those fields
% - appendValue (optional, boolean): whether to append new value to
%       existing field value, via horizontal vector concatenation
%
% Dependencies: bids-matlab (required), bids-matlab-tools (recommended for
% JSONio)

function jsonChangeValue(inputPath,fieldsToChange,newValues,appendValue)

if ~iscell(fieldsToChange)
    fieldsToChange = {fieldsToChange};
    if iscell(newValues) && length(newValues)>1     % Cell array of values for a single field
        newValues = {newValues};
    end
end
if ~iscell(newValues), newValues = {newValues}; end
if ~exist('appendValue','var') || isempty(appendValue)
    appendValue = 0;    % Whether to append new value to existing value
end
jsonOpts.indent = '\t';     % Use tab indentation for JSON outputs

[~,~,inputExt] = fpp.util.fileParts(inputPath);
if strcmpi(inputExt,'.json')
    error('fpp.bids.jsonChangeValue must be run on data file, not json file.');
end

inputJsonPath = fpp.bids.jsonPath(inputPath);
jsonData = fpp.bids.getMetadata(inputPath);

for f=1:length(fieldsToChange)
    if isempty(newValues{f}) && isfield(jsonData,fieldsToChange{f})
        jsonData = rmfield(jsonData,fieldsToChange{f});
    elseif appendValue && isfield(jsonData,fieldsToChange{f})
        eval(['jsonData.' fieldsToChange{f} ' = [jsonData.' fieldsToChange{f} ' newValues{f}];']);
    else
        eval(['jsonData.' fieldsToChange{f} ' = newValues{f};']);
    end
end

bids.util.jsonencode(inputJsonPath,jsonData,jsonOpts);

end