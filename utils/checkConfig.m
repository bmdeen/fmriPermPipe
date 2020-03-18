
function [configError, studyDir, fslPrefix] = checkConfig(study)

configError = '';
studyDir = '';
fslPrefix = '';

baseDir = getenv('FMRI_BASE_DIR');
if isempty(baseDir)
    configError = 'ERROR: FMRI_BASE_DIR not defined in bash environment.';
    return;
elseif ~exist(baseDir,'dir')
    configError = 'ERROR: fMRI base directory does not exist.';
    return;
end

studyDir = [baseDir '/' study];
if ~exist(studyDir,'dir')
    configError = ['ERROR: Study directory ' studyDir ' does not exist.'];
    return;
end

fslDir = getenv('FSL_DIR');
if isempty(fslDir)
    configError = 'ERROR: FSL_DIR is not defined in bash environment.';
    return;
end

% On some systems, it may be necessary to source FSL before running FSL
% commands via matlab, by uncommenting the below code.  Note that this will
% increase iteration time for model2ndPerm, so it is better to avoid if 
% possible.
% 
% fslPrefix = ['LD_LIBRARY_PATH=$LD_LIBRARY_PATH:' fslDir '; PATH=' fslDir '/bin:${PATH}; ' ...
%              '. ${FSLDIR}/etc/fslconf/fsl.sh; export PATH LD_LIBRARY_PATH; '];

fslCommands = {'fslmaths','fslmerge','fslsplit','fslroi','flirt','fnirt','applywarp',...
    'invwarp','convert_xfm','mcflirt','fslstats','fslval','fsl_tsplot','fslmeants',...
    'bet2','susan','cluster'};
cmdError = 0;
for c=1:length(fslCommands)
    [status,~] = system([fslPrefix fslCommands{c}]);
    if status==127
        cmdError=1;
        break;
    end
end
if cmdError==1
    configError = ['ERROR: FSL scripts not functioning properly from within '...
        'MATLAB. Check that FSL is properly installed and sourced.  If it is, '...
        'try uncommenting the fslPrefix code in checkConfig.'];
    return;
end

fsCommands = {'mri_convert','bbregister','tkregister2'};
cmdError = 0;
for c=1:length(fsCommands)
    [status,~] = system([fsCommands{c}]);
    if status==127
        cmdError=1;
        break;
    end
end
if cmdError==1 || ~exist('MRIread','file') || ~exist('MRIwrite','file')
    configError = ['ERROR: Freesurfer scripts not functioning properly from '...
        'within MATLAB. Check that Freesurfer is properly installed and sourced.'];
    return;
end

bxhCommands = {'analyze2bxh','bxh2analyze','bxhreorient'};
cmdError = 0;
for c=1:length(bxhCommands)
    [status,~] = system(['which ' bxhCommands{c}]);
    if status
        cmdError=1;
        break;
    end
end
if cmdError==1
    configError = 'ERROR: BXH/XCEDE scripts are not properly installed or sourced.';
    return;
end

end