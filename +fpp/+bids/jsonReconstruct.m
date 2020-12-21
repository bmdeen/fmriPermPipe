
% Function to copy a json metadata file, keeping only specified fields.
% Derives metadata from all relevant JSON files in hierarchy, and collapses
% specified fields into single output file. Based on BIDS 1.4.1 spec, with
% some fields from BIDS derivatives BEP014 included as options.
%
% fpp.bids.jsonReconstruct(inputJsonPath,outputJsonPath,fieldsToKeep)
%
% Arguments:
% - inputJsonPath (string): path to json file to copy
% - outputJsonPath (string): path to output json field
% - fieldsToKeep (optional, cell array of strings): JSON fields to copy
%    OR (string): label for image type, determining which fields to keep.
%    Options: midprepFMRI, fMRI, Mask, Seg, MRI, Surf, Cifti, Xfm, Fmap,
%    keepAll, keepAllPreproc
%
% Dependencies: bids-matlab (required), bids-matlab-tools (recommended for
% JSONio)

function jsonReconstruct(inputJsonPath,outputJsonPath,fieldsToKeep)

jsonOpts.indent = '\t';     % Use tab indentation for JSON outputs
keepAllFields = 0;          % Whether to keep all fields in output JSON file.

if ~exist('fieldsToKeep','var') || isempty(fieldsToKeep)
    % Default JSON fields to maintain:
    fieldsToKeep = {'Manufacturer','ManufacturersModelName','DeviceSerialNumber',...
        'StationName','SoftwareVersions','HardcopyDeviceSoftwareVersion','MagneticFieldStrength',...
        'ReceiveCoilName','ReceiveCoilActiveElements','GradientSetType','MRTransmitCoilSequence',...
        'MatrixCoilMode','CoilCombinationMethod','PulseSequenceType','ScanningSequence',...
        'SequenceVariant','ScanOptions','SequenceName','PulseSequenceDetails',...
        'NonlinearGradientCorrection','NumberShots','ParallelReductionFactorInPlane',...
        'ParallelAcquisitionTechnique','PartialFourier','PartialFourierDirection',...
        'PhaseEncodingDirection','EffectiveEchoSpacing','TotalReadoutTime','EchoTime',...
        'InversionTime','SliceTiming','SliceEncodingDirection','DwellTime','FlipAngle',...
        'MultibandAccelerationFactor','NegativeContrast','MultibandAccelerationFactor',...
        'AnatomicalLandmarkCoordinates','InstitutionName','InstitutionAddress',...
        'InstitutionalDepartmentName','ContrastBolusIngredient','RepetitionTime',...
        'VolumeTiming','TaskName','NumberOfVolumesDiscardedByScanner','NumberOfVolumesDiscardedByUser',...
        'DelayTime','AcquisitionDuration','DelayAfterTrigger','Instructions','TaskDescription',...
        'CogAtlasID','CogPOID','IntendedFor','Description','Sources','RawSources',...
        'SpatialReference','SkullStripped','Resolution','Density','Type','Atlas',...
        'Manual','LabelMap','Method','Multiplexed','Software','SoftwareVersion',...
        'Invertible','FromFile','ToFile','FromFileSHA','ToFileSHA','CommandLine',...
        'CoordinateSystem','ReferenceMap','NonstandardReference','ReferenceIndex',...
        'TransformFile','PipelineDescription','SourceDatasets','BIDSVersion','Name',...
        'DatasetType','Authors','Space','BandpassFilter','Neighborhood','Threshold',...
        'Transformations','ROI'};
elseif ischar(fieldsToKeep)
    fieldsToKeepLabel = fieldsToKeep;
    fieldsToKeep = {'Description','Sources','RawSources'};   % Keep for all derivatives
    switch lower(fieldsToKeepLabel)
        case 'midprepfmri'  % fMRI data in middle stages of preprocessing (volume, surface, or cifti)
            fieldsToKeep = [fieldsToKeep {'SpatialRef','SkullStripped',...
                'Resolution','Density','TaskName','RepetitionTime','DelayAfterTrigger',...
                'NumberOfVolumesDiscardedByScanner','NumberOfVolumesDiscardedByUser',...
                'EchoTime','EchoNumber','SliceTiming','PhaseEncodingDirection',...
                'SliceEncodingDirection','EffectiveEchoSpacing','TotalReadoutTime','DwellTime'}];
        case 'fmri'     % Preprocessed fMRI data (volume, surface, or cifti)
            fieldsToKeep = [fieldsToKeep {'SpatialRef','SkullStripped',...
                'Resolution','Density','TaskName','RepetitionTime','DelayAfterTrigger',...
                'NumberOfVolumesDiscardedByScanner','NumberOfVolumesDiscardedByUser'}];
        case 'mask'     % Mask image (volume, surface, or cifti)
            fieldsToKeep = [fieldsToKeep {'SpatialRef','Resolution','Density','Type','Atlas'}];
        case 'seg'      % Segmentation image (volume, surface, or cifti)
            fieldsToKeep = [fieldsToKeep {'SpatialRef','Resolution','Density','Manual','Atlas'}];
        case 'mri'      % Preprocessed 3D MRI data
            fieldsToKeep = [fieldsToKeep {'SpatialRef','SkullStripped','Resolution'}];
        case 'surf'     % Surface file (surf or metric - func/shape)
            fieldsToKeep = [fieldsToKeep {'SpatialRef','Density'}];
        case 'citfi'    % Cifti file (dseries, dscalar, etc)
            fieldsToKeep = [fieldsToKeep {'SpatialRef','SkullStripped','Resolution','Density'}];
        case 'xfm'      % Transformation file
            fieldsToKeep = [fieldsToKeep {'Software','SoftwareVersion','Invertible',...
                'Multiplexed','FromFile','ToFile','FromFileSHA','ToFileSHA','CommandLine'}];
        case 'fmap'     % Field map or spin-echo EPI
            fieldsToKeep = [fieldsToKeep {'SpatialRef','PhaseEncodingDirection',...
                'EffectiveEchoSpacing','TotalReadoutTime','IntendedFor','Units'}];
        case 'keepall'	% Keep all fields
            keepAllFields = 1;
        case 'keepallpreproc'   % Keep all fields involved in preprocessing
            fieldsToKeep = [fieldsToKeep {'SpatialRef','SkullStripped',...
                'Resolution','Density','TaskName','RepetitionTime','DelayAfterTrigger',...
                'NumberOfVolumesDiscardedByScanner','NumberOfVolumesDiscardedByUser',...
                'EchoTime','EchoNumber','SliceTiming','PhaseEncodingDirection',...
                'SliceEncodingDirection','EffectiveEchoSpacing','TotalReadoutTime','DwellTime',...
                'Type','Manual','Atlas','Software','SoftwareVersion','Invertible',...
                'Multiplexed','FromFile','ToFile','FromFileSHA','ToFileSHA','CommandLine',...
                'IntendedFor','Units'}];
        otherwise
            error('If fieldsToKeep is specified as a string, it must correspond to one of the allowed keywords.');
    end
end

if length(inputJsonPath)>=5 && ~strcmpi(inputJsonPath(end-4:end),'.json')
    inputJsonPath = fpp.bids.jsonPath(inputJsonPath);
end
if length(outputJsonPath)>=5 && ~strcmpi(outputJsonPath(end-4:end),'.json')
    outputJsonPath = fpp.bids.jsonPath(outputJsonPath);
end
if ~exist(inputJsonPath,'file')
    return;
end
inputJsonData = fpp.bids.getMetadata(inputJsonPath);
outputJsonData = struct();

if keepAllFields
    outputJsonData = inputJsonData;
else
    for f=1:length(fieldsToKeep)
        thisField = fieldsToKeep{f};
        if isfield(inputJsonData,thisField)
            eval(['outputJsonData.' thisField ' = inputJsonData.' thisField ';']);
        end
    end
end

bids.util.jsonencode(outputJsonPath,outputJsonData,jsonOpts);

end