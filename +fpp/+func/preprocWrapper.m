
% preprocWrapper(studyDir,subjects,varargin)
%
%%% CURRENTLY UNDERGOING MAJOR OVERHAUL
%%% CONVERTING TO BIDS NAMING CONVENTION
%
%
%
%
% 
% Wrapper for fpp.func.preprocMRI
%
% Example usage: fpp.func.preprocWrapper('/pathto/studyDirectory','SUB*','fwhm',4,'plotResults',0)
% 
% Required arguments:
% - studyDir (string): path to study directory
% - subjects (string or cell array of strings): list of subjects to
%       analyze. Can use asterisk-based regular expression. Examples:
%       {'SUB01','SUB02'} or 'SUB*'.
% 
% Variable arguments: (can also use any variable arg for preprocMRI)
% - expt (string): name of experiment to process.
% - runList (vector of integers in (0,Inf)): list of runs to process.
% - outputSuffix (string): suffix for .prep directory.  Can be used to run
%       the script multiple times with different options.

function preprocWrapper(studyDir,subjects,varargin)

% Load/check config variables.
configError = fpp.util.checkConfig;
if ~isempty(configError)
    fprintf('%s\n',configError);
    return;
end
if ~exist(studyDir,'dir')
    fprintf('%s\n',['ERROR: Study directory ' studyDir ' does not exist.']);
    return;
end

% Basic parameters
expt = [];                      % Experiment to analyze.
runList = [];                   % Runs to analyze.
outputSuffix = '';              % New suffix for preproc output

% Edit variable arguments.  Note: optInputs checks for proper input.
varArgList = {'expt','runList','outputSuffix'};
for i=1:length(varArgList)
    argVal = fpp.util.optInputs(varargin,varArgList{i});
    if ~isempty(argVal)
        eval([varArgList{i} ' = argVal;']);
    end
end

% Convert list of subjects from "dir" inputs to actual names
subjectsNew = {};
if ~iscell(subjects), subjects = {subjects}; end
for s=1:length(subjects)
    if ~ischar(subjects{s})
        fprintf('%s\n','ERROR: subject argument must be a string or cell array of strings.');
        return;
    end
    if strfind(subjects{s},'*')
        subjectsTmp = dir([studyDir '/' subjects{s}]);
        for j=1:length(subjectsTmp)
            subjectsNew{end+1} = subjectsTmp(j).name;
        end
    else
        subjectsNew{end+1} = subjects{s};
    end
end
subjects = subjectsNew;

for s = 1:length(subjects)
    
    subject = subjects{s};
    subjDir = [studyDir '/' subject];
    funcDir = [subjDir '/func'];
    analDir = [subjDir '/analysis'];
    targetInit = [funcDir '/' targetName '.nii.gz'];
    
    scanlogFiles = fpp.util.regExpDir(subjDir,'.*scanlo.*([^~]$)');
    
    if exist(subjDir,'dir')==0
        fprintf('%s\n\n',['Subject directory does not exist for ' subject '. Skipping this subject.']);
        continue;
    elseif isempty(scanlogFiles)
        fprintf('%s\n\n',['Missing scanlog file for subject ' subject '. Skipping this subject.']);
        continue;
    elseif exist(funcDir,'dir')==0
        fprintf('%s\n\n',['Functional directory does not exist for ' subject '. Skipping this subject.']);
        continue;
    elseif ~genTarget && exist(targetInit,'file')==0
        fprintf('%s\n\n',['Registration target does not exist for ' subject '. Skipping this subject.']);
        continue;
    end
    
    if exist(analDir,'dir')==0, mkdir(analDir); end
    
    for scan = 1:length(scanlogFiles)
        
        fileID = fopen([subjDir '/' scanlogFiles(scan).name]);
        fileData = textscan(fileID,'%d%s%d\n');
        fclose(fileID);
        dataLengths = cellfun(@length,fileData);
        if length(unique(dataLengths))>1 || sum(dataLengths==0)>1
            fprintf('%s\n',['Scanlog is invalid for ' subject ', scan ' int2str(scan) '. Skipping this scan.']);
            continue;
        end
        [acqs, expts, runs] = deal(fileData{1},fileData{2},fileData{3});
        
        for i = 1:length(acqs)
            
            % Skip run if it's anatomical or DTI.
            if sum(regexpi(expts{i},'anat'))>0 || sum(regexpi(expts{i},'dti'))>0
                continue;
            end
            
            % Skip run if it's not the specified experiment/run number.
            if ~isempty(expt) && ~strcmpi(expts{i},expt)
                continue;
            end
            if ~isempty(runList) && ~ismember(runs(i),runList)
                continue;
            end
            
            runSuffix = ['-' int2str(runs(i))];
            rfdInit = [funcDir '/' expts{i} runSuffix '.nii.gz'];
            if ~exist(rfdInit,'file')
                fprintf('%s\n',['Raw func data does not exist for ' subject ', expt ' expts{i} '-' int2str(runs(i)) '.']);
                continue;
            end
            
            outputDir = [analDir '/' expts{i} runSuffix outputSuffix '.prep'];
            
            fpp.func.preperoc(rfdInit,outputDir,varargin{:});
            
        end
    end
end

end