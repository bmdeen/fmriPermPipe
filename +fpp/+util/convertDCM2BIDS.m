
% Function to convert DCMs from a single scan session to BIDS-formatted
% NIFTI and JSON files (wrapper for dcm2niix).
%
% Arguments:
%   subject (string): subject ID
%   dicomDir (string): input directory containing DCM files to convert
%   outputDir (string): directory in which to output NII files (BIDS
%       subject level)
%   scanlogPath (string): path to scanlog TSV file, specifying what each
%       dicom series corresponds to (see below)
%   overwrite (boolean, optional): whether to overwrite existing NIFTI
%       files
%
% Dependencies:
%   - dcm2niix
%   - bids-matlab
%   - bids-matlab-tools (optional, recommended for JSONio)
%
% Usage Notes:
%
%   scanlog.tsv has required fields series and type, and optional fields
%   ses, task, acq, dir, and run. Series corresponds to DICOM series
%   number, and type corresponds to BIDS image suffix. EG:
% 
%   [sub-01_scan-01_scanlog.tsv]
%   series  type	task    dir     run
%	11  T1w	n/a	n/a 1
%	12	T2w	n/a	n/a	1
%	13	epi	n/a	AP	1
%	14	epi	n/a	PA	1
%	20	sbref	facelocalizer	n/a	1
%	21	bold	facelocalizer	n/a	1
%	22	sbref	facelocalizer	n/a	2
%	23	bold	facelocalizer	n/a	2
%	24	sbref	facelocalizer	n/a	3
%	25	bold	facelocalizer	n/a	3
%	30	dwi	n/a	n/a	1
%	38	dwi	n/a	n/a	2
%
%   Unspecified values can be defined as "n/a" or "-"
%
%   NOTE: BIDS key "ses" is distinct from the label "scan," which corresponds
%   corresponds to a specific scan sesssion. This way, multiple scans can 
%   be combined into a single session, or one scan can be split into
%   multiple sessions if desired. Each scan must have its own DICOM folder
%   and scanlog.tsv file, while ses is defined as a field within the
%   scanlog file.
%

function convertDCM2BIDS(subject,dicomDir,outputDir,scanlogPath,overwrite)

if ~exist('overwrite','var'), overwrite = 0; end
jsonOpts.indent = '\t';     % Use tab indentation for JSON outputs

% Load and validate scanlog TSV
scanlog = bids.util.tsvread(scanlogPath);
if ~isfield(scanlog,'series') || ~isfield(scanlog,'type')
    error('Scanlog file must have series and type fields.');
end
if ~isnumeric(scanlog.series) || any(mod(scanlog.run,1)~=0) || any(scanlog.run<1)
    error('Series must be defined as positive integers for all acquisitions.');
end
if isfield(scanlog,'ses') && (any(strcmp(scanlog.ses,'-')) || any(strcmpi(scanlog.ses,'n/a')))
    error('If sessions are specified, they must be defined for all series.');
end
if ~exist('overwrite','var')
    overwrite = 0;
end

% Convert all fields of scanlog, other than series, to cell arrays of
% strings.
scanlogFields = setdiff(fieldnames(scanlog),'series');
for f=1:length(scanlogFields)
    if ~iscell(eval(['scanlog.' scanlogFields{f}]))
        eval(['scanlog.' scanlogFields{f} '=num2cell(scanlog.' scanlogFields{f} ');']);
    end
    for i=1:length(scanlog.series)
        value = eval(['scanlog.' scanlogFields{f} '{' int2str(i) '}']);
        if isnumeric(value)
            if mod(value,1)==0  % value is an integer
                if value<100
                    value = fpp.util.numPad(value,2);
                else
                    value = int2str(value);
                end
            else
                value = num2str(value);
            end
            eval(['scanlog.' scanlogFields{f} '{' int2str(i) '}=value;']);
        end
    end
end

% Run dcm2niix to convert DCM to NII
niixDir = [outputDir '/dcm2niix_tmp'];
mkdir(niixDir);
fpp.util.system(['dcm2niix -f %s -z y -b y -ba y -o ' niixDir ' ' dicomDir]);

% Convert to BIDS naming and filefpp.util.system structure
indSpinEcho = []; outputJsonPathsSpinEcho = {};
for i=1:length(scanlog.series)
  
    % Determine image types
    imageType = scanlog.type{i};
    if contains(imageType,'T1') || contains(imageType,'T2') ...
            || contains(imageType,'PD') || contains(imageType,'FLAIR') ...
            || contains(imageType,'FLASH') || contains(imageType,'angio')
        imageTypeID(i) = 1;
    elseif contains(imageType,'bold')
        imageTypeID(i) = 2;
    elseif contains(imageType,'dwi')
        imageTypeID(i) = 3;
    elseif contains(imageType,'fmap') || contains(imageType,'epi') ...
            || contains(imageType,'phasediff') || contains(imageType,'magnitude') ...
            || contains(imageType,'fieldmap')
        imageTypeID(i) = 4;
    elseif contains(imageType,'sbref')
        imageTypeID(i) = 5;
    else
        continue;   % Skip unrecognized image types
    end
    
    % Determine key/value pairs specified by scanlog.
    keys = {'sub'};
    values = {subject};
    % Parse additional keys in scanlog file
    if isfield(scanlog,'ses')
        outputSubDir = [outputDir '/ses-' scanlog.ses{i}];
        keys{end+1} = 'ses';
        values{end+1} = scanlog.ses{i};
    else
        outputSubDir = outputDir;
    end
    if isfield(scanlog,'task') && ~(strcmp(scanlog.task{i},'-') || strcmpi(scanlog.task{i},'n/a'))
        keys{end+1} = 'task';
        values{end+1} = scanlog.task{i};
    end
    if isfield(scanlog,'acq') && ~(strcmp(scanlog.acq{i},'-') || strcmpi(scanlog.acq{i},'n/a'))
        keys{end+1} = 'acq';
        values{end+1} = scanlog.acq{i};
    end
    if isfield(scanlog,'dir') && ~(strcmp(scanlog.dir{i},'-') || strcmpi(scanlog.dir{i},'n/a'))
        keys{end+1} = 'dir';
        values{end+1} = scanlog.dir{i};
    end
    if isfield(scanlog,'run') && ~(strcmp(scanlog.run{i},'-') || strcmpi(scanlog.run{i},'n/a'))
        keys{end+1} = 'run';
        values{end+1} = scanlog.run{i};
    end
    
    switch imageTypeID(i)
        case 1
            outputSubDir = [outputSubDir '/anat'];
        case 2
            outputSubDir = [outputSubDir '/func'];
        case 3
            outputSubDir = [outputSubDir '/dwi'];
        case 4
            outputSubDir = [outputSubDir '/fmap'];
        case 5
            if strcmp(scanlog.type{i+1},'bold')
                outputSubDir = [outputSubDir '/func'];
                sbrefForBold = 1;
            elseif strcmp(scanlog.type{i+1},'dwi')
                outputSubDir = [outputSubDir '/dwi'];
                sbrefForBold = 0;
            else
                continue;
            end
    end
    
    if ~exist(outputSubDir,'dir'), mkdir(outputSubDir); end
    
    % If multiple echoes exist, process each separately
    didNotWrite(i) = 0;
    meFiles = fpp.util.regExpDir([niixDir '/' int2str(scanlog.series(i)) '_e*.nii.gz'],[int2str(scanlog.series(i)) '_e[1-9].nii.gz']);
    if length(meFiles)>1
        keys{end+1} = 'echo';
        values{end+1} = 1;
        for e=1:length(meFiles)
            values{end} = int2str(e);
            inputNiftiPath = [niixDir '/' meFiles(e).name];
            inputJsonPath = strrep(inputNiftiPath,'.nii.gz','.json');
            outputNiftiPath = [outputSubDir '/' fpp.bids.changeName('',keys,values,imageType,'.nii.gz')];
            outputNiftiPathsRelative{i,e} = strrep(outputNiftiPath,[outputDir '/'],'');
            outputJsonPath = strrep(outputNiftiPath,'.nii.gz','.json');
            if ~exist(outputNiftiPath,'file') || overwrite
                fpp.util.system(['cp ' inputNiftiPath ' ' outputNiftiPath]);
                fpp.util.system(['cp ' inputJsonPath ' ' outputJsonPath]);
                if imageTypeID(i)==2
                    json = bids.util.jsondecode(outputJsonPath);
                    json.TaskName = scanlog.task{i};
                    bids.util.jsonencode(outputJsonPath,json,jsonOpts);
                elseif imageTypeID(i)==5 && sbrefForBold
                    json = bids.util.jsondecode(outputJsonPath);
                    json.TaskName = scanlog.task{i+1};
                    bids.util.jsonencode(outputJsonPath,json,jsonOpts);
                end
            else
                didNotWrite(i) = 1;
                break;
            end
        end
    else
        inputNiftiPath = [niixDir '/' int2str(scanlog.series(i)) '.nii.gz'];
        inputJsonPath = strrep(inputNiftiPath,'.nii.gz','.json');
        outputNiftiPath = [outputSubDir '/' fpp.bids.changeName('',keys,values,imageType,'.nii.gz')];
        outputNiftiPathsRelative{i,1} = strrep(outputNiftiPath,[outputDir '/'],'');
        outputJsonPath = strrep(outputNiftiPath,'.nii.gz','.json');
        if ~exist(outputNiftiPath,'file') || overwrite
            fpp.util.system(['cp ' inputNiftiPath ' ' outputNiftiPath]);
            fpp.util.system(['cp ' inputJsonPath ' ' outputJsonPath]);
            if imageTypeID(i)==2
                json = bids.util.jsondecode(outputJsonPath);
                json.TaskName = scanlog.task{i};
                bids.util.jsonencode(outputJsonPath,json,jsonOpts);
            elseif imageTypeID(i)==5 && sbrefForBold
                json = bids.util.jsondecode(outputJsonPath);
                json.TaskName = scanlog.task{i+1};
                bids.util.jsonencode(outputJsonPath,json,jsonOpts);
            end
        else
            didNotWrite(i) = 1;
        end
    end
    
    % For DTI data, copy bvec and bval files
    if imageTypeID(i)==3
        inputBvalPath = strrep(inputNiftiPath,'.nii.gz','.bval');
        inputBvecPath = strrep(inputNiftiPath,'.nii.gz','.bvec');
        outputBvalPath = strrep(outputNiftiPath,'.nii.gz','.bval');
        outputBvecPath = strrep(outputNiftiPath,'.nii.gz','.bvec');
        if exist(inputBvalPath,'file') && (~exist(outputBvalPath,'file') || overwrite)
            fpp.util.system(['cp ' inputBvalPath ' ' outputBvalPath]);
            fpp.util.system(['cp ' inputBvecPath ' ' outputBvecPath]);
        end
    end
    
    % For spin echo "field map" data, keep track of indices and json file
    % names, to subsequently edit IntendedFor field.
    if contains(imageType,'epi')
        indSpinEcho(end+1) = i;
        outputJsonPathsSpinEcho{end+1} = outputJsonPath;
    end
end

% Define IntendedFor field for spin echo "field map" acquisitions, as all
% func, dwi, and sbref images following the current set of consecutive spin
% echo series, and before the next set. This way, if multiple field maps
% were acquired during a scan session, the most recent will be used for
% distortion correction. Assumes that maps with opposing phase-encode
% directions were acquired consecutively (or at least, without any
% intervening BOLD/DWI acquisitions). NOTE: This is based on order of
% acquisitions listed in the scanlog file (not series number).
indFuncAndDwi = find(ismember(imageTypeID,[2 3 5]));
for i=1:length(indSpinEcho)
    if didNotWrite(i), continue; end
    ind = indSpinEcho(i);
    indIntendedFor = indFuncAndDwi(indFuncAndDwi>ind);
    if ~isempty(indIntendedFor) && i<length(indSpinEcho)
        tmp = indSpinEcho(i+1:end);
        indFirstSpinEchoAfterFunc = min(tmp(tmp>indIntendedFor(1)));
        if ~isempty(indFirstSpinEchoAfterFunc)
            indIntendedFor = indIntendedFor(indIntendedFor<indFirstSpinEchoAfterFunc);
        end
    end
    if ~isempty(indIntendedFor)
        intendedFor = outputNiftiPathsRelative(indIntendedFor,:)';
        intendedFor = intendedFor(:);
        intendedFor = intendedFor(~cellfun('isempty',intendedFor));
        json = bids.util.jsondecode(outputJsonPathsSpinEcho{i});
        json.IntendedFor = intendedFor;
        if length(json.IntendedFor)==1, json.IntendedFor = json.IntendedFor{1}; end
        bids.util.jsonencode(outputJsonPathsSpinEcho{i},json,jsonOpts);
    end
end

fpp.util.system(['rm -rf ' niixDir]);

end
