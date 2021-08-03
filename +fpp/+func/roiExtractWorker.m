
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
[~,~,roiExt] = fpp.util.fileParts(roiPath);
if ~ismember(lower(roiExt),{'.nii.gz','.nii','.dscalar.nii'})
    error('roiPath must be a NIFTI or CIFTI dscalar file.');
end
imageType = 'volume';
if strcmpi(roiExt,'.dscalar.nii')
    imageType = 'cifti';
end

% Determine model suffix and condition names
[~,modelName,~] = fpp.util.fileParts(modelDir);
suffix = fpp.bids.checkNameValue(modelName,'desc');
condPath = [modelDir '/' fpp.bids.changeName(modelName,[],[],'conditions','.tsv')];
condTSV = bids.util.tsvread(condPath);
condNames = condTSV.cond_names;

% Extract PSC values within ROI
for c=1:length(condNames)
    pscPath = [modelDir '/' fpp.bids.changeName(modelName,'desc',[suffix condNames{c}],'psc',roiExt)];
    pscVal = fpp.wb.command([imageType '-stats'],pscPath,[],[],['-reduce MEAN -roi ' roiPath]);
    psc(c) = str2num(strtrim(pscVal));
end

end