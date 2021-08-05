
% [psc,condNames] = fpp.func.roiExtractWorker(modelDir,roiPath)
%
% Basic function that extracts percent signal change responses across
% conditions from one run of one task, within a region-of-interest
% specified by a binary mask file.
%
% Arguments:
% - modelDir (string): modelperm directory to extract responses from
% - roipath (string): path to region-of-interest mask image

function [psc,condNames] = roiExtractWorker(modelDir,roiPath)

psc = [];

% Check if ROI input is NIFTI or CIFTI
[~,roiName,roiExt] = fpp.util.fileParts(roiPath);
if ~ismember(lower(roiExt),{'.nii.gz','.nii','.dscalar.nii'})
    error('roiPath must be a NIFTI or CIFTI dscalar file.');
end
imageType = 'volume';
isCifti = 0;
if strcmpi(roiExt,'.dscalar.nii')
    imageType = 'cifti';
    isCifti = 1;
end

% Determine model suffix and condition names
[~,modelName,~] = fpp.util.fileParts(modelDir);
suffix = fpp.bids.checkNameValue(modelName,'desc');
roiDesc = fpp.bids.checkNameValue(roiName,'desc');
condPath = [modelDir '/' fpp.bids.changeName(modelName,[],[],'conditions','.tsv')];
condTSV = bids.util.tsvread(condPath);
condNames = condTSV.cond_names;

% For CIFTI files, if ROI lacks volume and statistical outputs have it,
% generate a temporary ROI file for mean extraction.
tmpROI = 0;
if isCifti
    pscPath = [modelDir '/' fpp.bids.changeName(modelName,'desc',[suffix condNames{1}],'psc',roiExt)];
    pscHasVol = fpp.util.checkMRIProperty('HasVol',pscPath);
    roiHasVol = fpp.util.checkMRIProperty('HasVol',roiPath);
    if fpp.util.checkMRIProperty('HasVol',pscPath) && ~fpp.util.checkMRIProperty('HasVol',roiPath)
        tmpVolPath = [modelDir '/' fpp.bids.changeName(roiName,{'den','desc'},...
            {'',[roiDesc 'tmpRoiExtract230951835']},[],'.nii.gz')];
        tmpROIPath = [modelDir '/' fpp.bids.changeName(roiName,{'den','desc'},...
            {'',[roiDesc 'tmpRoiExtract230951835']})];
        % Extract volume from pscPath, multiply by zero
        fpp.wb.command('cifti-separate',pscPath,'COLUMN',[],['-volume-all ' tmpVolPath]);
        fpp.fsl.maths(tmpVolPath,'-mul 0',tmpVolPath);
        % Construct ROI with zero volume
        fpp.wb.command('cifti-create-dense-from-template',pscPath,[],tmpROIPath,...
            ['-cifti ' roiPath ' -volume-all ' tmpVolPath]);
        roiPath = tmpROIPath;   % set to tmp ROI path
        tmpROI = 1;
    end
end

% Extract PSC values within ROI
for c=1:length(condNames)
    pscPath = [modelDir '/' fpp.bids.changeName(modelName,'desc',[suffix condNames{c}],'psc',roiExt)];
    pscVal = fpp.wb.command([imageType '-stats'],pscPath,[],[],['-reduce MEAN -roi ' roiPath]);
    psc(c) = str2num(strtrim(pscVal));
end

if tmpROI
    fpp.util.deleteImageAndJson({tmpROIPath,tmpVolPath});
end

end