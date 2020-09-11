
% Function to change file path to .json file.

function outputPath = jsonPath(inputPath)

    [inputDir,inputName,~] = fpp.util.fileParts(inputPath);
    if ~isempty(inputDir)
        outputPath = [inputDir '/' inputName '.json'];
    else
        outputPath = [inputName '.json'];
    end

end