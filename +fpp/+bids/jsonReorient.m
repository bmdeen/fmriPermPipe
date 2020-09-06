
% Script to alter fields in a .json metadata file pertinent to image
% orientation.
%
% Arguments:
%   inputJsonPath (string): path to json file to edit
%   inputOrientation (string): orientation string for input (e.g. LAS, RPI)
%   outputOrientation (string): orientation string for output
%
% Modifies fields: PhaseEncodingDirection, SliceEncodingDirection
% Removes fields: InPlanePhaseEncodingDirectionDICOM, ImageOrientationPatientDICOM
%
% Script conventions: x,y,z are dimension in physical space, i,j,k are
%   dimensions in image coordinate space.

function jsonReorient(inputJsonPath,inputOrientation,outputOrientation)

inputOrientation = upper(inputOrientation);
outputOrientation = upper(outputOrientation);
orientationOptions = {'RAS','LAS','RPS','LPS','RAI','LAI','RPI','LPI',...
                      'RSA','LSA','RSP','LSP','RIS','LIA','RIP','LIP',...
                      'ARS','ALS','PRS','PLS','ARI','ALI','PRI','PLI',...
                      'ASR','ASL','PSR','PSL','AIR','AIL','PIR','PIL',...
                      'SRA','SLA','SRP','SLP','IRA','ILA','IRP','ILP',...
                      'SAR','SAL','SPR','SPL','IAR','IAL','IPR','IPL'};
ijk = 'ijk';

% Check if input/output orientations are valid
if ~ismember(inputOrientation,orientationOptions) || ~ismember(outputOrientation,orientationOptions)
    error('ERROR: Orientation inputs are not valid three-letter orientation strings, e.g. ''LAS''');
end

% Load JSON file
json = bids.util.jsondecode(inputJsonPath);

% Remove fields that will no longer be accurate
if isfield(json,'ImageOrientationPatientDICOM')
    json = rmfield(json,'ImageOrientationPatientDICOM');
end
if isfield(json,'InPlanePhaseEncodingDirectionDICOM')
    json = rmfield(json,'InPlanePhaseEncodingDirectionDICOM');
end

% Compute which of x,y,z physical dimensions correspond to i (1), j (2), and k (3) image dimensions
dimsInputXYZ = [findstr(inputOrientation,'L') findstr(inputOrientation,'R') findstr(inputOrientation,'A') ...
    findstr(inputOrientation,'P') findstr(inputOrientation,'S') findstr(inputOrientation,'I')];
dimsOutputXYZ = [findstr(outputOrientation,'L') findstr(outputOrientation,'R') findstr(outputOrientation,'A') ...
    findstr(outputOrientation,'P') findstr(outputOrientation,'S') findstr(outputOrientation,'I')];
% Convert to which of i,j,k image dimensions correspond to x (1), y (2), and z (3) physical dimensions
for d=1:3
    dimsInputIJK(d) = find(dimsInputXYZ==d);
    dimsOutputIJK(d) = find(dimsOutputXYZ==d);
end

% Check which image space dimensions (i, j, k) have switched direction
for d=1:3
    dimSwappedIJK(d) = ~strcmp(inputOrientation(d),outputOrientation(dimsOutputIJK==dimsInputIJK(d)));
end

if isfield(json,'PhaseEncodingDirection')
    newPhaseString = ijk(dimsOutputIJK==dimsInputIJK(findstr(json.PhaseEncodingDirection(1),ijk)));
%     if (dimSwappedIJK(findstr(json.PhaseEncodingDirection(1),ijk)) && isempty(findstr(json.PhaseEncodingDirection,'-')))...
%             || (~dimSwappedIJK(findstr(json.PhaseEncodingDirection(1),ijk)) && ~isempty(findstr(json.PhaseEncodingDirection,'-')))
    if xor(dimSwappedIJK(findstr(json.PhaseEncodingDirection(1),ijk)),~isempty(findstr(json.PhaseEncodingDirection,'-')))
        newPhaseString(end+1) = '-';
    end
    json.PhaseEncodingDirection = newPhaseString;
end

if isfield(json,'SliceEncodingDirection')
    newSliceString = ijk(dimsOutputIJK==dimsInputIJK(findstr(json.SliceEncodingDirection(1),ijk)));
    if ~xor(dimSwappedIJK(findstr(json.SliceEncodingDirection(1),ijk)),isempty(findstr(json.SliceEncodingDirection(1),'-')))
        newSliceString(end+1) = '-';
    end
    json.SliceEncodingDirection = newSliceString;
end

bids.util.jsonencode(inputJsonPath,json);

end