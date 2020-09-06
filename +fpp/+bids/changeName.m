
% Function to change BIDS file name: key/value pairs, suffix, and
% extension.
% 
% - Can remove entities by setting new value to an empty string, and can
%   build a BIDS filename from empty string input.
% - If input filename has at least one BIDS key (sub, ses, etc), it will be
%   considered BIDS-formatted. Any unrecognized entities will be left at
%   of the filename (this includes the file suffix, and other unrecognized
%   key-value pairs, e.g. used by derivatives).
% - Deals with non-BIDS-formatted inputs by adding entities and a default 
%   suffix to the end of the file name. Any input with no recognized BIDS 
%   keys is considered non-BIDS-formatted (note: this includes names like 
%   "bold.nii.gz," because the script has no knowledge about which strings
%   are valid BIDS suffices.)
%
% Arguments:
%   inputPath (string) - filename to modify
%   keysToChange (cell array or string) - keys to modify
%   newValues (cell array or string) - new values for each specified key
%   newSuffix (string, optional) - n
%   newExtension (string, optional) - new file extension

function outputPath = changeName(inputPath,keysToChange,newValues,newSuffix,newExtension)

defaultSuffix = 'bold';
changeSuffix = 0;
if ~exist('keysToChange','var'), keysToChange = {};
elseif ~iscell(keysToChange), keysToChange = {keysToChange}; end
if ~exist('newValues','var'), newValues = {};
elseif ~iscell(newValues), newValues = {newValues}; end
if exist('newSuffix','var') && ~isempty(newSuffix) && ischar(newSuffix)
    defaultSuffix = newSuffix;
    changeSuffix = 1;
end
if ~isempty(inputPath)
    [inputDir,inputName,inputExtension] = fpp.util.fileParts(inputPath);
else
    inputDir = '';
    inputName = '';
    inputExtension = '.nii.gz';
end
if ~exist('newExtension','var') || ~ischar(newExtension) || ~strcmp(newExtension(1),'.')
    newExtension = inputExtension;
end

% Potential key names, in order of precedence in filename. "Key index"
% refers to position in this list.
keys = {'sub','ses','tpl','task','acq','ce','rec','dir','run','mod','echo','recording',...
    'proc','hemi','space','volspace','atlas','split','res','den','label','model',...
    'parameter','desc','subset','from','to','mode'};

% Check indices of keys to change, remove improper keys
keysToChangeInd = [];
for k=1:length(keysToChange)
    tmp = find(strcmp(keys,keysToChange{k}));
    if isempty(tmp), tmp = 0; end
    keysToChangeInd(k) = tmp;
end
keysToChange(keysToChangeInd==0) = [];
newValues(keysToChangeInd==0) = [];
keysToChangeInd(keysToChangeInd==0) = [];

% Decompose name into entities (and suffix), by separating substrings
% between underscores
entities = split(inputName,'_');

% Check which # entity each key is in
keyIndByEntity = inf(length(entities),1);   % Put substrings with no key (suffix or unrecognized) at the end of filename
for e=1:length(entities)
    for k=1:length(keys)
        if strfind(entities{e},[keys{k} '-'])
            keyIndByEntity(e) = k;
        end
    end
end

% If no key-value pairs are found, assume input filename is not 
% BIDS-formatted. In this case, add new entities and suffix to the end of
% the filename.
if sum(keyIndByEntity<Inf)==0
    [~,sortInd] = sort(keysToChangeInd);
    newValues = newValues(sortInd);
    entities = {};
    for i=1:length(newValues)
        if ~isempty(newValues{i})
            entities{end+1} = [keysToChange{i} '-' newValues{i}];
        end
    end
    entityString = join(entities,'_');
    entityString = entityString{1};
    
    outputPath = [];
    if ~isempty(inputDir), outputPath = [inputDir '/']; end
    if ~isempty(inputName)
        outputPath = [outputPath inputName];
    end
    if ~isempty(entityString)
        if ~isempty(inputName)
            outputPath = [outputPath '_' entityString];
        else
            outputPath = [outputPath entityString];
        end
    end
    if ~(isempty(inputName) && isempty(entityString))
        outputPath = [outputPath '_'];
    end
    outputPath = [outputPath defaultSuffix newExtension];
    return;
end

% Change values
indEntitiesToRemove = [];
for k=1:length(keysToChange)
    if ismember(keysToChangeInd(k),keyIndByEntity)  % key is already in filename
        if ~isempty(newValues{k})
            entities{keyIndByEntity==keysToChangeInd(k)} = [keysToChange{k} '-' newValues{k}];
        else
            indEntitiesToRemove = [indEntitiesToRemove find(keyIndByEntity==keysToChangeInd(k))];
        end
    elseif ~isempty(newValues{k})                   % add key to end of filename
        entities{end+1} = [keysToChange{k} '-' newValues{k}];
        keyIndByEntity(end+1) = keysToChangeInd(k);
    end
end
entities(indEntitiesToRemove) = [];
keyIndByEntity(indEntitiesToRemove) = [];

% Resort entities, so that anything new is put into proper order
[~,sortInd] = sort(keyIndByEntity);
entities = entities(sortInd);

% Change suffix, if new suffix was specified
if changeSuffix
    entities{end} = defaultSuffix;
end

% Recombine entities into file name
outputName = join(entities,'_');
outputName = outputName{1};

% Define full output path
if ~isempty(inputDir), inputDir = [inputDir '/']; end
outputPath = [inputDir outputName newExtension];
    
end
