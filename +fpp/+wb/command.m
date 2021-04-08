
% cmdOut = fpp.wb.command(cmdType,inputPath,argText,outputPath,flagText[,outputDescription,appendDescription])
%
% Wrapper for HCP's wb_command function, which generates JSON metadata for
% the output data if it exists for the input.
%
% This script will call:
%   wb_command -cmdType inputPath argText outputPath flagText
%
% Optional arguments:
% - outputDescription (string): new Description in output json metadata
% - appendDescription (boolean): whether outputDescription should be
%       appended to existing description
%
% inputPath and outputPath should be paths to individual files, while
% argText/flagText can contain multiple flags, terms and/or paths.
%
% Run wb_command -help for more information on functionality.
%
% NOTE: Function doesn't currently generate .json output for commands:
%   cifti-create-scalar-series, volume-create
%

function cmdOut = command(cmdType,inputPath,argText,outputPath,flagText,...
    outputDescription,appendDescription)

jsonOpts.indent = '\t';     % Use tab indentation for JSON outputs
removeBidsDir = @(x) fpp.bids.removeBidsDir(x); % Define wrapper function for fpp.bids.removeBidsDir (for cellfun functionality)

% Process arguments
if ~exist('argText','var') || isempty(argText)
    argText = '';
end
if ~exist('outputPath','var') || isempty(outputPath)
    outputPath = '';
end
if ~exist('flagText','var') || isempty(flagText)
    flagText = '';
end
if ~exist('outputDescription','var') || isempty(outputDescription)
    outputDescription = '';
end
if ~exist('appendDescription','var') || isempty(appendDescription)
    appendDescription = 0;
end

% Load information about wb_command subcommand
[fppWbDir,~,~]		= fileparts(mfilename('fullpath'));			% path to the directory containing this script
cmdInfo = bids.util.tsvread([fppWbDir '/wb_command_info.tsv']);
if strcmp(cmdType(1),'-'), cmdType = cmdType(2:end); end
cmdInd = find(strcmpi(cmdType,cmdInfo.Command));
if isempty(cmdInd)
    error('Invalid cmdType inpute. See wb_command -list-commands for valid options.');
end
% Define input/output type for this command.
inputType = cmdInfo.InputType{cmdInd};
outputType = cmdInfo.OutputType{cmdInd};



%%% PART 1: RUN COMMAND

% Add quotes to math expression if needed
if sum(regexp(cmdType,'-math'))>0
    if ~isempty(inputPath) && isempty(argText)
        argText = inputPath;
        inputPath = [];
    end
    if ~strcmp(argText(1),'''')
        argText = ['''' argText];
    end
    if ~strcmp(argText(end),'''')
        argText = [argText ''''];
    end
end

cmd = ['wb_command -' cmdType ' ' inputPath ' ' argText ' ' outputPath ' ' flagText];
[~,cmdOut] = fpp.util.system(cmd);

% For commands that output text and no file, return output value.
if strcmp(outputType,'Disp')
    return;
end



%%% PART 2: PROCESS JSON METADATA

% Convert inputs to cell arrays
if ~isempty(inputPath), inputPath = {inputPath}; end
if ~isempty(outputPath), outputPath = {outputPath}; end
allArgs = {};
if ~isempty(argText), allArgs = [allArgs; split(argText,' ')]; end
if ~isempty(flagText), allArgs = [allArgs; split(flagText,' ')]; end

% Adjust input list relevant to json output
if cmdInfo.NonstandardInput(cmdInd) && cmdInfo.NonstandardOutput(cmdInd)
    % For commands listed below, need to redefine inputs and outputs:
    % convert-affine, convert-warpfield, file-convert, metric-convert
    switch lower(cmdType)
        case {'convert-affine','convert-warpfield'}
            inputType = 'Xfm';
            outputType = 'Xfm';
            inputOpts = {'-from-world','-from-itk','-from-flirt','-from-fnirt'};
            outputOpts = {'-to-world','-to-itk','-to-flirt','-to-fnirt'};
            for i=1:length(allArgs)
                if strcmp(allArgs{i},inputOpts)
                    inputInd = i+1;
                end
                if strcmp(allArgs{i},outputOpts)
                    outputInd = i;
                end
            end
        case 'file-convert'
            cmdOpts = {'-border-version-convert','-nifti-version-convert',...
                '-cifti-version-convert'};
            types = {'Border','Volume','Cifti'};
            for i=1:length(allArgs)
                if strcmp(allArgs{i},cmdOpts)
                    inputInd = i+1;
                    outputInd = i+3;
                    inputType = types(strcmp(allArgs{i},cmdOpts));
                    outputType = inputType;
                    break;
                end
            end
        case 'metric-convert'
            for i=1:length(allArgs)
                if strcmp(allArgs{i},'-to-nifti')
                    inputInd = i+1;
                    outputInd = i+2;
                    inputType = 'Metric';
                    outputType = 'Volume';
                    break;
                elseif strcmp(allArgs{i},'-from-nifti')
                    inputInd = i+1;
                    outputInd = i+3;
                    inputType = 'Volume';
                    outputType = 'Metric';
                    break;
                end
            end
    end
    inputPath = allArgs(inputInd);
    outputPath = allArgs(outputInd);
    
elseif cmdInfo.NonstandardInput(cmdInd)
    % Find gifti, nifti, cifti, .border, .scene, .spec, .wbsparse files in 
    % args/flags, add these to inputs.
    [files,types] = findData(allArgs);
    inputPath = [inputPath files];
    if strcmp(inputType,'n/a') && ~isempty(types)
        inputType = types{1};
    end
    
    % Shift the new input to be indexed first, so json comes from this:
    if ismember(cmdType,{'create-signed-distance-volume','foci-get-projection-vertex'})
        inputPath{end+1} = inputPath{1};
        inputPath(1) = [];
    end
    % Remove the first input, which shouldn't contribute to json output
    if ismember(cmdType,{'gifti-convert','label-to-border','metric-extrema'',''metric-false-correlation',...
            'metric-fill-holes','metric-find-clusters','metric-gradient','metric-remove-islands',...
            'metric-rois-from-extrema','metric-rois-to-border','metric-smoothing','metric-tfce'})
        inputPath(1) = [];
    end
end

% Return if no input file or no input json metadata
if isempty(inputPath) && ~ismember(cmdType,{'surface-create-sphere',...
        'cifti-create-scalar-series','volume-create'})
    return;
elseif ~isempty(inputPath) && isempty(fpp.bids.getMetadata(inputPath{1}))
    return;
end

% Modify list of outputs that need json files, based on additional args
if cmdInfo.NonstandardOutput(cmdInd) && ~cmdInfo.NonstandardInput(cmdInd)
    switch lower(cmdType)
        case {'set-structure','set-map-names'}
            % Set output to input
            outputPath = inputPath;
        otherwise
            % Find gifti, nifti, cifti, .border, .scene, .spec, .wbsparse files in 
            % args/flags, add these to outputs.
            [files,types] = findData(allArgs);
            outputPath = [outputPath files];
            if strcmp(outputType,'n/a') && ~isempty(types)
                outputType = types{1};
            end
    end
end

% If output doesn't need JSON file, return
if sum(strcmp(outputType,{'Zip','Text','Image','n/a'}))>0, return; end

% Determine fields to keep in output json, based on output type
% Note: most data types use midprepfmri option, which uses a general set of
% preprocessing-specific json fields
switch outputType
    case 'Volume'
        fieldsToKeep = 'midprepfmri';
    case 'Surface'
        fieldsToKeep = 'midprepfmri';
    case 'Label'
        fieldsToKeep = 'seg';
    case 'Metric'
        fieldsToKeep = 'midprepfmri';
    case 'Gifti'
        fieldsToKeep = 'midprepfmri';
    case 'Cifti'
        fieldsToKeep = 'midprepfmri';
    case 'Xfm'
        fieldsToKeep = 'xfm';
    case 'Scene'
        fieldsToKeep = 'keepall';   % No BIDS standard for these files, so leave all fields
    case 'Spec'
        fieldsToKeep = 'keepall';
    case 'Annot'
        fieldsToKeep = 'keepall';
    case 'Border'
        fieldsToKeep = 'seg';
    case 'Foci'
        fieldsToKeep = 'keepall';
    case 'Fiber'
        fieldsToKeep = 'midprepfmri';
    case 'Wbsparse'
        fieldsToKeep = 'keepall';
    case 'File'
        fieldsToKeep = 'keepall';
end
if ~isempty(inputPath)
    for i=1:length(outputPath)
        
        if strcmp(inputPath{1},outputPath{i}), continue; end  % Don't reconstruct json if output file has the same name as input
        
        fpp.bids.jsonReconstruct(inputPath{1},outputPath{i},fieldsToKeep);
        
        if length(inputPath)>1
            % For multiple inputs, define cell array of Sources / Raw Sources
            fpp.bids.jsonChangeValue(outputPath{i},'Sources',...
                cellfun(removeBidsDir,inputPath,'UniformOutput',false));
            rawSources = {};
            for j=1:length(inputPath)
                jsonData = fpp.bids.getMetadata(inputPath{j});
                if isfield(jsonData,'RawSources')
                    if iscell(jsonData.RawSources)
                        rawSources = [rawSources jsonData.RawSources];
                    else
                        rawSources = [rawSources {jsonData.RawSources}];
                    end
                end
            end
            if ~isempty(rawSources)
                fpp.bids.jsonChangeValue(outputPath{i},'RawSources',rawSources);
            end
        else
            fpp.bids.jsonChangeValue(outputPath{i},'Sources',removeBidsDir(inputPath{1}));
        end
    end
end

% Modify output json metadata where necessary, for first output file
if isempty(outputPath), return; end
if cmdInfo.NonstandardJsonDef(cmdInd)
    switch lower(cmdType)
        % Resampling: update SpatialRef, Resolution, and Density based on target
        case {'border-resample','label-resample','surface-cut-resample',...
                'metric-resample','surface-resample','volume-to-surface-mapping',...
                'volume-label-to-surface-mapping','volume-parcel-resampling',...
                'volume-parcel-resampling-generic','label-to-volume-mapping',...
                'metric-to-volume-mapping','volume-affine-resample',...
                'volume-warpfield-resample','cifti-resample','cifti-resample-dconn-memory'}
            if ismember(cmdType,{'volume-to-surface-mapping','volume-label-to-surface-mapping'})
                newTemplate = allArgs{1};
            else
                newTemplate = allArgs{2};
            end
            jsonDataTemplate = fpp.bids.getMetadata(newTemplate);
            % Compute resolution for volume/cifti outputs
            if ismember(cmdType,{'volume-parcel-resampling','volume-parcel-resampling-generic',...
                'label-to-volume-mapping','metric-to-volume-mapping','volume-affine-resample',...
             	'volume-warpfield-resample','cifti-resample','cifti-resample-dconn-memory'})
                useRes = 1;
            else
                useRes = 0;
            end
            % Compute density for surface/cifti outputs
            if ismember(cmdType,{'border-resample','label-resample','surface-cut-resample',...
                'metric-resample','surface-resample','volume-to-surface-mapping',...
                'volume-label-to-surface-mapping','cifti-resample','cifti-resample-dconn-memory'})
                useDen = 1;
            else
                useDen = 0;
            end
            if useRes && isfield(jsonDataTemplate,'Resolution')
                fpp.bids.jsonChangeValue(outputPath{1},'Resolution',jsonDataTemplate.Resolution);
            end
            if useDen && isfield(jsonDataTemplate,'Density')
                fpp.bids.jsonChangeValue(outputPath{1},'Density',jsonDataTemplate.Density);
            end
            if isfield(jsonDataTemplate,'SpatialRef')
                fpp.bids.jsonChangeValue(outputPath{1},'SpatialRef',jsonDataTemplate.SpatialRef);
            else
                fpp.bids.jsonChangeValue(outputPath{1},'SpatialRef',fpp.bids.removeBidsDir(newTemplate));
            end
            
            
        % Surface modification: set SpatialRef to output file
        case {'surface-flip-lr','surface-flip-normals'}
            fpp.bids.jsonChangeValue(outputPath{1},'SpatialRef',fpp.bids.removeBidsDir(outputPath{1}));
            
        % NOTE: Not currently changing SpatialRef for the following
        % commands, but this might depend on what exactly is meant by the
        % SpatialRef field - just surface coords, or their embedding in a
        % in a specific spatial coordinate frame).
        case {'surface-modify-sphere','surface-set-coordinates',...
             'surface-smoothing','surface-sphere-project-unproject',...
             'surface-apply-affine','surface-apply-warpfield','surface-match'}
            % Do nothing
            
            
        % Surface generation: define Density, set SpatialRef to output file
        case 'surface-create-sphere'
            jsonData.SpatialRef = fpp.bids.removeBidsDir(outputPath{1});
            jsonData.Density = [inputPath argText];
            bids.util.jsonencode(fpp.bids.jsonPath(outputPath{1}),jsonData,jsonOpts);
            
            
        % Cifti generation: combine SpatialRef info from volume and surfaces
        case {'cifti-create-dense-scalar','cifti-create-dense-timeseries','cifti-create-label'}
            spatialRefData = [];
            volInd = find(strcmp('-volume',allArgs))+1;
            if ~isempty(volInd)
                volJsonData = fpp.bids.getMetadata(allArgs{volInd});
                if isfield(volJsonData,'SpatialRef')
                    spatialRefData.VolumeReference = volJsonData.SpatialRef;
                end
            end
            lhInd = [find(strcmp('-left-metric',allArgs)) find(strcmp('-left-label',allArgs))] + 1;
            if ~isempty(lhInd)
                lhJsonData = fpp.bids.getMetadata(allArgs{lhInd});
                if isfield(lhJsonData,'SpatialRef')
                    spatialRefData.CIFTI_STRUCTURE_CORTEX_LEFT = lhJsonData.SpatialRef;
                end
            end
            rhInd = [find(strcmp('-right-metric',allArgs)) find(strcmp('-right-label',allArgs))] + 1;
            if ~isempty(rhInd)
                rhJsonData = fpp.bids.getMetadata(allArgs{rhInd});
                if isfield(rhJsonData,'SpatialRef')
                    spatialRefData.CIFTI_STRUCTURE_CORTEX_RIGHT = rhJsonData.SpatialRef;
                end
            end
            if ~isempty(spatialRefData)
                fpp.bids.jsonChangeValue(outputPath{1},'SpatialRef',spatialRefData);
            end
            
            
        % Cifti separation: modify SpatialRef of each individual output.
        case 'cifti-separate'
            % Define indices of -volume-all, -label, -metric, and -volume
            % flags in argument list.
            [volAllInd,labelInd,metricInd,volumeInd] = deal([],[],[],[]);
            [volFiles,lhFiles,rhFiles] = deal({},{},{});
            volAllInd = find(strcmp('-volume-all',allArgs));
            labelInd = find(strcmp('-label',allArgs));
            for i=1:length(labelInd)    % Remove -label following -volume-all
                if length(allArgs)<labelInd(i)+2 || strcmp(allArgs{labelInd(i)+2}(1),'-')
                    labelInd(i) = [];
                    break;
                end
            end
            metricInd = find(strcmp('-metric',allArgs));
            volumeInd = find(strcmp('-volume',allArgs));
            flagTypes = [ones(1,length(volAllInd)) 2*ones(1,length(labelInd))...
                3*ones(1,length(metricInd)) 4*ones(1,length(volumeInd))];
            allInd = [volAllInd; labelInd; metricInd; volumeInd];
            [allInd,sortInd] = sort(allInd);
            flagTypes = flagTypes(sortInd);
            allInd = [allInd; length(allArgs)];  % Add index of final argument
            % Loop through flags, add images to volfiles, lhFiles, rhFiles cell arrays
            for i=1:length(allInd)-1
                ind = allInd(i);
                switch flagTypes(i)
                    case 1  % -volume-all
                        volFiles = [volFiles allArgs(ind+2)];
                        if ind+4 <= allInd(i+1)
                            ind2 = find('-roi',allArgs(ind+3:allInd(i+1)-1))+ind+2;
                            if ~isempty(ind2), volFiles = [volFiles allArgs(ind2+1)]; end
                            ind2 = find('-label',allArgs(ind+3:allInd(i+1)-1))+ind+2;
                            if ~isempty(ind2), volFiles = [volFiles allArgs(ind2+1)]; end
                        end
                    case {2,3}  % -label, -metric
                        if strcmp(allArgs{ind+1},'CORTEX_LEFT')
                            lhFiles = [lhFiles allArgs(ind+2)];
                            if ind+4 <= allInd(i+1)
                                ind2 = find('-roi',allArgs(ind+3:allInd(i+1)-1))+ind+2;
                                if ~isempty(ind2), lhFiles = [lhFiles allArgs(ind2+1)]; end
                            end
                        elseif strcmp(allArgs{ind+1},'CORTEX_RIGHT')
                            rhFiles = [rhFiles allArgs(ind+2)];
                            if ind+4 <= allInd(i+1)
                                ind2 = find('-roi',allArgs(ind+3:allInd(i+1)-1))+ind+2;
                                if ~isempty(ind2), rhFiles = [rhFiles allArgs(ind2+1)]; end
                            end
                        end
                    case 4  % -volume
                        volFiles = [volFiles allArgs(ind+2)];
                        if ind+4 <= allInd(i+1)
                            ind2 = find('-roi',allArgs(ind+3:allInd(i+1)-1))+ind+2;
                            if ~isempty(ind2), volFiles = [volFiles allArgs(ind2+1)]; end
                        end
                end
            end
            % Load Cifti json data, check for relevant fields
            jsonData = fpp.bids.getMetadata(inputPath{1});
            if isfield(jsonData,'SpatialRef') && isfield(jsonData.SpatialRef,'VolumeReference')
                volRef = jsonData.SpatialRef.VolumeReference;
                for f=1:length(volFiles)
                    fpp.bids.jsonChangeValue(volFiles{f},'SpatialRef',volRef);
                end
            end
            if isfield(jsonData,'SpatialRef') && isfield(jsonData.SpatialRef,'CIFTI_STRUCTURE_CORTEX_LEFT')
                lhRef = jsonData.SpatialRef.CIFTI_STRUCTURE_CORTEX_LEFT;
                for f=1:length(lhFiles)
                    fpp.bids.jsonChangeValue(lhFiles{f},'SpatialRef',lhRef);
                end
            end
            if isfield(jsonData,'SpatialRef') && isfield(jsonData.SpatialRef,'CIFTI_STRUCTURE_CORTEX_RIGHT')
                rhRef = jsonData.SpatialRef.CIFTI_STRUCTURE_CORTEX_RIGHT;
                for f=1:length(rhFiles)
                    fpp.bids.jsonChangeValue(rhFiles{f},'SpatialRef',rhRef);
                end
            end
            
            
        % Affine transform modification and generation
        case {'convert-affine','convert-warpfield'}
            jsonData = fpp.bids.getMetadata(outputPath{1});
            if isfield(jsonData,'CommandLine')
                fpp.bids.jsonChangeValue(outputPath{1},'CommandLine',['; ' cmd],1);
            end
        case 'volume-warpfield-affine-regression'
            jsonData = fpp.bids.getMetadata(outputPath{1});
            if isfield(jsonData,'CommandLine')
                fpp.bids.jsonChangeValue(outputPath{1},'CommandLine',['; ' cmd],1);
            end
            if isfield(jsonData,'Software')
                fpp.bids.jsonChangeValue(outputPath{1},'Software','; wb_command -volume-warpfield-affine-regression',1);
            else
                fpp.bids.jsonChangeValue(outputPath{1},'Software','wb_command -volume-warpfield-affine-regression',1);
            end
            fpp.bids.jsonChangeValue(outputPath{1},{'Type','Invertible'},{'12-dof affine',true});
        case 'surface-affine-regression'
            jsonData.Type = '12 dof affine';
            jsonData.Software = 'wb_command -surface-affine-regression';
            jsonData.Invertible = true;
            jsonData.FromFile = fpp.bids.removeBidsDir(inputPath);
            jsonData.ToFile = fpp.bids.removeBidsDir(allArgs{1});
            jsonData.CommandLine = cmd;
            jsonData.Description = 'Affine transformation file generated by wb_command -surface-affine-regression.';
            bids.util.jsonencode(fpp.bids.jsonPath(outputPath{1}),jsonData,jsonOpts);
            
            
        % Volumetric reorientation: reorient json, define SpatialRef as itself
        case 'volume-reorient'
            inOrientation = fpp.util.getImageOrientation(inputPath{1});
            outOrientation = argsAll{1};
            fpp.bids.jsonReorient(outputPath{1},inOrientation,outOrientation);
            fpp.bids.jsonChangeValue(outputPath{1},'SpatialRef',fpp.bids.removeBidsDir(outputPath{1}));
        case 'volume-set-space'
            inOrientation = fpp.util.getImageOrientation(inputPath{1});
            outOrientation = fpp.util.getImageOrientation(outputPath{1});
            fpp.bids.jsonReorient(outputPath{1},inOrientation,outOrientation);
            fpp.bids.jsonChangeValue(outputPath{1},'SpatialRef',fpp.bids.removeBidsDir(outputPath{1}));
            
            
        % Not currently implemented
        case 'cifti-create-scalar-series'
            % Generates CIFTI from text.
            warning('Note: -cifti-create-scalar-series option does not currently generate JSON metadata.');
        case 'volume-create'
            % Generates volume 
            warning('Note: -volume-create option does not currently generate JSON metadata.');
            
    end
end

% Modify description, if specified
if ~isempty(outputDescription)
    for i=1:length(outputPath)
        fpp.bids.jsonChangeValue(outputPath{i},'Description',outputDescription,appendDescription);
    end
end

end

% Auxiliary function to find arguments that are data file paths
function [files,types] = findData(allArgs)

extensionsToFind = {'nii','gii','gz','border','scene','spec','wbsparse'};
typesToFind = {'Cifti','Gifti','Volume','Border','Scene','Spec','Wbsparse'};
files = {};
types = {};   % Infered output image type. Doesn't need to be exact, just used to determine json metadata to copy

% Loop through each argument, find files
for i=1:length(allArgs)
    [~,~,ext] = fileparts(allArgs{i});
    extInd = find(strcmpi(ext(2:end),extensionsToFind));    % Check if argument has neuroimaging file extension
    if ~isempty(extInd)
        files = [files allArgs(i)];
        types = [types typesToFind(extInd)];
    end
end

end


