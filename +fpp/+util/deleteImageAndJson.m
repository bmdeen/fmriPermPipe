
% Simple utility to delete an image and its associated json metadata file.
%
% fpp.util.deleteImageAndJson(inputPath)

function deleteImageAndJson(inputPath)

fpp.util.system(['rm -rf ' inputPath]);
inputJsonPath = fpp.bids.jsonPath(inputPath);
if exist(inputJsonPath,'file')
    fpp.util.system(['rm -rf ' inputJsonpath]);
end

end