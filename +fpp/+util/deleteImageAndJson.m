
% Simple utility to delete an image and its associated json metadata file.
%
% fpp.util.deleteImageAndJson(inputPath)
%
% inputPath (string or cell array of strings): path(s) to image file(s) to
%   delete

function deleteImageAndJson(inputPath)

if iscell(inputPath)
    for i=1:length(inputPath)
        fpp.util.deleteImageAndJson(inputPath{i});
    end
else
    fpp.util.system(['rm -rf ' inputPath]);
    inputJsonPath = fpp.bids.jsonPath(inputPath);
    if exist(inputJsonPath,'file')
        fpp.util.system(['rm -rf ' inputJsonPath]);
    end
end

end