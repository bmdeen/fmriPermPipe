
% Function to change or add BIDS key/value pairs in filename.

function outputPath = bidsChangeEntity(inputPath,keysToChange,newValues)

% Get filename subentities separated by "_". If no sub entities have "key-" at
% the start, just add relevant keys to the end of the filename. Otherwise,
% alter/add keys, maintaining order.

if ~iscell(keysToChange), keysToChange = {keysToChange}; end
if ~iscell(newValues), newValues = {newValues}; end
[inputDir,inputName,inputExt] = filepartsGZ(inputPath);

% Potential key names, in order
keys = {'sub','ses','task','acq','ce','rec','dir','run','mod','echo','recording',...
    'proc','space','split','desc'};

% Check index of keys to change, remove improper keys
for k=1:length(keysToChange)
    tmp = find(contains(keys,keysToChange{k}));
    if ~isempty(tmp)
        keysToChangeInd(k) = tmp;
    else
        keysToChangeInd(k) = 0;
    end
end
keysToChange(keysToChangeInd==0) = [];
newValues(keysToChangeInd==0) = [];
keysToChangeInd(keysToChangeInd==0) = [];

% Decompose name into entities and suffix, separated by underscore
entities = {};
underscoreInd = strfind(inputName,'_');
if ~isempty(underscoreInd)
    for e=1:length(underscoreInd)+1

        if e==1
            entities{e} = inputName(1:underscoreInd(e)-1);
        elseif e==length(underscoreInd)+1
            entities{e} = inputName(underscoreInd(e-1)+1:end);
        else
            entities{e} = inputName(underscoreInd(e-1)+1:underscoreInd(e)-1);
        end
    end
end

% Check which # entity each key is in
keyIndByEntity = inf(length(entities),1);   % Strings with no key (suffix or unrecognized) at the end
for e=1:length(entities)
    for k=1:length(keys)
        if strfind(entities{e},[keys{k} '-'])
            keyIndByEntity(e) = k;
        end
    end
end
disp(keyIndByEntity)

% If no key-value pairs are found, input filename is not BIDS-formated. In
% this case, just add entities to the end of the filename, assuming _bold
% suffix.
if sum(keyIndByEntity<Inf)==0
    entityString = '';
    for k=1:length(keysToChange)
        entityString = [entityString '_' keysToChange{k} '-' newValues{k}];
    end
    outputPath = [inputDir '/' inputName entityString '_bold' inputExt];
    return;
end

% Change values
for k=1:length(keysToChange)
    if ismember(keysToChangeInd(k),keyIndByEntity)
        entities{keyIndByEntity==keysToChangeInd(k)} = [keysToChange{k} '-' newValues{k}];
    else
        entities{end+1} = [keysToChange{k} '-' newValues{k}];
        keyIndByEntity(end+1) = keysToChangeInd(k);
    end
end

% Resort entities, so that anything new is put into proper order
[~,sortInd] = sort(keyIndByEntity);
entities = entities(sortInd);

outputName = entities{1};
for e=2:length(entities)
    outputName = [outputName '_' entities{e}];
end

if ~isempty(inputDir), inputDir = [inputDir '/']; end
outputPath = [inputDir outputName inputExt];
    
end
