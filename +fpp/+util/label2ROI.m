
% Function to convert a NIFTI/CIFTI "label" or discrete segmentation image
% to an ROI, bincluding multiple subregions. Wrapper for wb_command
% -volume-label-to-roi and -cifti-label-to-roi. Outputs mask.nii.gz or
% mask.dscalar.nii file.
%
% fpp.util.label2ROI(parcPath,parcInds,outputPath)

function label2ROI(parcPath,parcInds,outputPath)

[outputDir,outputName,~] = fpp.util.fileParts(outputPath);
if isempty(outputDir)
    outputDir = pwd;
end
[~,~,inputExt] = fpp.util.fileParts(parcPath);
if ismember(lower(inputExt),{'.nii.gz','.nii'})
    inputType = 'volume';
elseif strcmpi(inputExt,'.dlabel.nii')
    inputType = 'cifti';
else
    error('parcPath input must be a NIFTI or dlabel.nii CIFTI file.');
end

sumEquation = '';
flagText = '';

% Generate ROIs for each specified label.
for i=1:length(parcInds)
    tmpPath{i} = [outputDir '/' outputName '_tmpLabel2ROI' int2str(i) '_20398511243890' inputExt];
    if strcmp(inputType,'volume')
        % Use fslmaths in case input isn't formatted as a Connectome Workbench label volume
        fpp.fsl.maths(parcPath,['-uthr ' int2str(parcInds(i)) ' -thr ' int2str(parcInds(i)) ' -bin'],tmpPath{i});
    else
        fpp.wb.command([inputType '-label-to-roi'],parcPath,[],tmpPath{i},['-key ' int2str(parcInds(i))]);
    end
    if i==1
        sumEquation = [sumEquation 'r1'];
    else
        sumEquation = [sumEquation '+r' int2str(i)];
    end
    flagText = [flagText ' -var r' int2str(i) ' ' tmpPath{i}];
end

% Add individual ROIs
fpp.wb.command([inputType '-math'],[],sumEquation,outputPath,flagText);

% Modify json metadata
if exist(fpp.bids.jsonPath(outputPath),'file')
    fpp.bids.jsonChangeValue(outputPath,'Sources',fpp.bids.removeBidsDir(parcPath));
end

% Remove temporary files
for i=1:length(parcInds)
    fpp.util.system(['rm -rf ' tmpPath{i}]);
    if exist(fpp.bids.jsonPath(tmpPath{i}),'file')
        fpp.util.system(['rm -rf ' fpp.bids.jsonPath(tmpPath{i})]);
    end
end

end