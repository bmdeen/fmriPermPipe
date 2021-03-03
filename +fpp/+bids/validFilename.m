
% Function to check if a filename is BIDS-formatted
%
% isValid = fpp.bids.validFilename(inputPath)

function isValid = validFilename(inputPath)

isValid = true;
[~,inputName,~] = fpp.util.fileParts(inputPath);
entities = split(inputName,'_');
for e=1:length(entities)-1
    if sum(regexp(entities{e},'[0-9a-zA-Z]+-[0-9a-zA-Z]+'))==0
        isValid = false; return;
    end
end
if sum(regexp(entities{end},'[0-9a-zA-Z]+'))==0
    isValid = false; return;
end

end