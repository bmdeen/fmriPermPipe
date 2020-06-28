
% Function to convert Freesurfer-style para files to behav files used by
% fMRIPermPipe, with fsl_inputs variable specifying timing info.
%
% Arguments:
% inputPath: path to .para file to convert
% outputPath: path to .mat output file
% con_info: struct array with field "vals" and "name" that specifyings
%   default contrasts to use for analysis. Vals: vector of contrast
%   coefficients.
% nColumns (optional): number of columns in .para file. Default (4): 
%   columns for start time, condition, duration, condition label. Other 
%   option (5): columns for start time, condition, duration, regressor 
%   height, condition label.
% condsToRemove (optional): vector of conditions to remove, e.g. rest
%   conditions that shouldn't have a regressor. Script removes condition
%   with label 0 by default.
% newCondNames (optional): cell array of new names for each condition
% convertToBlock (optional): whether to combine adjacent rows of the para
%   file with the same condition label into a single block. To be used for
%   para files with a row for each TR rather than block.
% overwrite (optional): boolean, whether to overwrite existingg files

function para2Behav(inputPath,outputPath,con_info,nColumns,condsToRemove,newCondNames,convertToBlock,overwrite)

if ~(exist('nColumns','var') && isscalar(nColumns) && ismember(nColumns,[4 5]))
    nColumns = 4;
end
if ~(exist('condsToRemove','var') && isnumeric(condsToRemove) && isvector(condsToRemove))
    condsToRemove = [];
end
if ~(exist('convertToBlock','var') && (isnumeric(convertToBlock) || islogical(convertToBlock)) && ...
        isscalar(convertToBlock) && ismember(convertToBlock,[0 1]))
    convertToBlock = 0;
end

if ~(exist('overwrite','var') && (isnumeric(overwrite) || islogical(overwrite)) && ...
        isscalar(overwrite) && ismember(overwrite,[0 1]))
    overwrite = 0;
end

if convertToBlock==1
    fprintf('%s\n',['Converting single-TR rows to blocks. NOTE: this will fail '...
        'if the design contains adjacent blocks of the same condition.']);
end

if exist(outputPath,'file') && ~overwrite
    fprintf('%s\n',['Output path ' outputPath ' already exists.']);
    return;
end

% Need to modify this: include option for 4-column 

if nColumns==4
    [para(:,1), para(:,2), para(:,3), condList] = textread(inputPath,'%f%d%f%s');
    para(:,4) = ones(size(para,1),1);
else    % Five-column para file
    [para(:,1), para(:,2), para(:,3), para(:,4), condList] = textread(inputPath,'%f%d%f%f%s');
end
para(ismember(para(:,2),condsToRemove),2) = 0;      % Remove conditions listed in condsToRemove
%condList = condList(para(:,2)~=0);
%para = para(para(:,2)~=0,:);

for c=1:length(setdiff(unique(para(:,2)),0))
    tmpind = find(para(:,2)==c);
    condNames{c} = condList{tmpind(1)};
end
if exist('newCondNames','var') && iscell(newCondNames) && ...
        length(newCondNames)==length(condNames)
    condNames = newCondNames;
end

% Concatenate successive lines in the para file with the same condition,
% for para files that separate cond label for each TR.
if convertToBlock
    newPara = [];
    blockDuration = 0;
    for i=1:size(para,1)
        if i==1
            newPara(end+1,:) = para(i,:);
            lastCond = para(i,2);
        elseif para(i,2)~=lastCond % New condition!
            newPara(end,3) = blockDuration;
            blockDuration = 0;
            lastCond = para(i,2);
            newPara(end+1,:) = para(i,:);
        end
        blockDuration = blockDuration + para(i,3);
    end
    para = newPara;
end

for c = 1:length(condNames)
    cpara = para(para(:,2)==c,:);

    fsl_inputs(c).name = condNames{c};
    fsl_inputs(c).regressor = cpara(:,[1 3 4]);
end

save(outputPath,'fsl_inputs','con_info','condNames');

end