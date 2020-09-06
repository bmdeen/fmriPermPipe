
% Function to copy a json metadata file, keeping only specified fields.
% Derives metadata from all relevant JSON files in hierarchy, and collapses
% specified fields into single output file. Based on BIDS 1.4.0 spec.
%
% Arguments:
% - inputJsonPath (string): path to json file to copy
% - outputJsonPath (string): path to output json field
% - fieldsToKeep (cell array of strings): JSON fields to copy
%
% Dependencies: bids-matlab (required), bids-matlab-tools (recommended for
% JSONio)

function jsonReconstruct(inputJsonPath,outputJsonPath,fieldsToKeep)

jsonOpts.indent = '\t';     % Use tab indentation for JSON outputs

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
        'Manual','LabelMap','Method'};
end

if ~exist(inputJsonPath,'file')
    return;
end
inputJsonData = fpp.bids.getMetadata(inputJsonPath);
outputJsonData = struct();

for f=1:length(fieldsToKeep)
    thisField = fieldsToKeep{f};
    if isfield(inputJsonData,thisField)
        eval(['outputJsonData.' thisField ' = inputJsonData.' thisField ';']);
    end
end

bids.util.jsonencode(outputJsonPath,outputJsonData,jsonOpts);

end