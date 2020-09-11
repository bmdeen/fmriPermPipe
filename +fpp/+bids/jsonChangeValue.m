
% Function to change values in specific fields of a JSON file. Can also
% remove fields, by setting the corresponding element of newValues to an
% empty vector.
%
% Arguments:
% - inputJsonPath (string): path to JSON file, or corresponding data
% - fieldsToChange (cell array of strings): fields to edit
% - newValues (cell array): new values for those fields
% - appendValue (optional, boolean): whether to append new value to
%       existing field value, via horizontal vector concatenation

function jsonChangeValue(inputJsonPath,fieldsToChange,newValues,appendValue)

if ~iscell(fieldsToChange), fieldsToChange = {fieldsToChange}; end
if ~iscell(newValues), newValues = {newValues}; end
if ~exist('appendValue','var') || isempty(appendValue)
    appendValue = 0;    % Whether to append new value to existing value
end
jsonOpts.indent = '\t';     % Use tab indentation for JSON outputs

if length(inputJsonPath)>=5 && ~strcmpi(inputJsonPath(end-4:end),'.json')
    inputJsonPath = fpp.bids.jsonPath(inputJsonPath);
end
if ~exist(inputJsonPath,'file'), return; end

jsonData = bids.util.jsondecode(inputJsonPath);

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