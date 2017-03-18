
% convertDCM(study,subjects,varargin)
% 
% Converts DCM files to NIFTI images, labeled by experiment and run number,
% and realign to LAS orientation.  Loops through subjects.
% 
% Example usage: convertDCM('studyName','SUB*','overwrite',true)
% 
% SETUP: Prior to using this function with a given scanning system, the 
%        line beginning "dcm = dir" must be edited. See inline comments.
% 
% Required arguments:
% - study (string): name of the study to analyze, used to define the study 
%       directory as [getenv('FMRI_BASE_DIR') '/' study].
% - subjects (string or cell array of strings): list of subjects to
%       analyze. Can use asterisk-based regular expression. Examples:
%       {'SUB01','SUB02'} or 'SUB*'.
% 
% Variable arguments:
% - overwrite (boolean; default=0): whether to overwrite files that have
%       already been written by this function.

function convertDCM(study,subjects,varargin)

addpath([strrep(mfilename('fullpath'),mfilename,'') '/utils']);

% Load/check config variables.
[configError, studyDir, fslPrefix] = checkConfig(study);
if ~isempty(configError)
    fprintf('%s\n',configError);
    return;
end;

overwrite = false;

% Convert list of subjects from "dir" inputs to actual names
subjectsNew = {};
if ~iscell(subjects), subjects = {subjects}; end;
for s=1:length(subjects)
    if ~ischar(subjects{s})
        fprintf('%s\n','ERROR: subject argument must be a string or cell array of strings.');
        return;
    end;
    if strfind(subjects{s},'*')
        subjectsTmp = dir([studyDir '/' subjects{s}]);
        for j=1:length(subjectsTmp)
            subjectsNew{end+1} = subjectsTmp(j).name;
        end;
    else
        subjectsNew{end+1} = subjects{s};
    end;
end;
subjects = subjectsNew;

% Edit variable arguments.  Note: optInputs checks for proper input.
varArgList = {'overwrite'};
for i=1:length(varArgList)
    argVal = optInputs(varargin,varArgList{i});
    if ~isempty(argVal)
        eval([varArgList{i} ' = argVal;']);
    end;
end;

for s = 1:length(subjects)
    
    subject = subjects{s};
    subjDir = [studyDir '/' subject];
    
    % Note: script assumes that dicom directories and scanlogs are
    % consistently ordered by "dir" function.  Can be accomplished with a
    % consistent naming convention: scanlog, scanlog2, etc and dicom,
    % dicom2, etc.
    scanlogFiles = regExpDir(subjDir,'.*scanlo.*([^~]$)');
    dcmDirs = dir([subjDir '/dicom*']);
    
    if exist(subjDir,'dir')==0
        fprintf('%s\n',['Subject directory does not exist for ' subject '. '...
            'Skipping this subject.']);
        continue;
    elseif isempty(scanlogFiles)
        fprintf('%s\n\n',['Missing scanlog file for subject ' subject '. Skipping this subject.']);
        continue;
    elseif length(dcmDirs)~=length(scanlogFiles)
        fprintf('%s\n',['Missing scanlog file or dicom directory for subject ' ...
            subject '. Skipping this subject.']);
        continue;
    end;
    
    anatDir = [subjDir '/anat'];
    funcDir = [subjDir '/func'];
    if exist(anatDir,'dir')==0, mkdir(anatDir); end;
    if exist(funcDir,'dir')==0, mkdir(funcDir); end;
    
    for scan = 1:length(scanlogFiles)
        
        [acqs, expts, runs] = textread([subjDir '/' scanlogFiles(scan).name],'%d%s%d\n');
        
        dcmDir = [subjDir '/' dcmDirs(scan).name];
        
        for i = 1:length(acqs)
            
            acq = acqs(i);
            runSuffix = ['-' int2str(runs(i))];
            
            if sum(regexpi(expts{i},'dti'))>0
                outputDir = [subjDir '/dti' runSuffix];
                if ~exist(outputDir,'dir'), mkdir(outputDir); end;
                outputPath = [outputDir '/' expts{i} runSuffix '.nii.gz'];
            elseif sum(regexpi(expts{i},'anat'))>0
                outputPath = [anatDir '/' expts{i} runSuffix '.nii.gz'];
            else
                outputPath = [funcDir '/' expts{i} runSuffix '.nii.gz'];
            end;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% NOTE: the below line must be changed to fit your system's
            %%% dicom naming convention. This dir command should return
            %%% only the first dicom image for series number acq.
            %%% Asterisk-based regular expressions can be used, and the 
            %%% utility numPad may also be useful: numPad(input,len) 
            %%% returns a string consisting of the integer input, padded
            %%% with zeros on the left up to length len.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            dcm = dir([dcmDir '/*_S*' numPad(acq,3) 'I*001.DCM']);
            
            if length(dcm)>1
                fprintf('%s\n\n',['Multiple DCMs fit string template for subject ' ...
                    subject ', acquisition #' int2str(acq) '. Skipping this run. Check '
                    '"dcm = dir..." code in convertDCM.']);
                continue;
            elseif isempty(dcm)
                fprintf('%s\n\n',['No DCMs found for subject ' subjects{s} ', '...
                    'acquisition #' int2str(acq) '. Skipping this run. Check '...
                    '"dcm = dir..." code in convertDCM.']);
                continue;
            end;
            
            dcmPath = [dcmDir '/' dcm.name];
            
            if ~exist(outputPath,'file') || overwrite
                
                % Convert DCM to NIFTI
                system(['mri_convert ' dcmPath ' ' outputPath]);
                
                % Change orientation to LAS
                system(['analyze2bxh ' outputPath ' ' funcDir '/tmp_img_realign.bxh']);
                system(['bxhreorient --orientation=LAS ' funcDir '/tmp_img_realign.bxh ' funcDir '/tmp_img_realign-1.bxh']);
                system(['bxh2analyze --niigz ' funcDir '/tmp_img_realign-1.bxh ' funcDir '/tmp_img_realign-2']);
                system(['rm -rf ' outputPath '; mv ' funcDir '/tmp_img_realign-2.nii.gz ' outputPath]);
                system(['rm -rf ' funcDir '/tmp_img_realign*']);
                
                fprintf('%s\n',['Converted data for subject ' subject ', expt: ' expts{i} runSuffix '.']);
        
            end;
            
        end;
    
    end;
end;