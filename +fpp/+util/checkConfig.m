
function checkConfig

configError = '';
cmdPrefix = 'export LD_LIBRARY_PATH=``; ';  % Prevent matlab from sourcing its own libraries

fslDir = getenv('FSL_DIR');
if isempty(fslDir)
    configError = 'ERROR: FSL_DIR is not defined in bash environment.';
end

fslCommands = {'fslmaths','fslmerge','fslsplit','fslroi','flirt','fnirt','applywarp',...
    'invwarp','convert_xfm','mcflirt','fslstats','fslval','fsl_tsplot','fslmeants',...
    'bet2','susan','cluster'};
cmdError = 0;
for c=1:length(fslCommands)
    [status,~] = system([cmdPrefix fslCommands{c}]);
    if status==127
        cmdError=1;
        break;
    end
end
if cmdError==1
    configError = ['ERROR: FSL scripts not functioning properly from within '...
        'MATLAB. Check that FSL is properly installed and sourced.  You may '...
        'need to run MATLAB from a terminal on some systems.'];
end

fsCommands = {'mri_convert','mris_convert','bbregister','tkregister2'};
cmdError = 0;
for c=1:length(fsCommands)
    [status,~] = system([cmdPrefix fsCommands{c}]);
    if status==127
        cmdError=1;
        break;
    end
end
if cmdError==1 || ~exist('MRIread','file') || ~exist('MRIwrite','file')
    configError = ['ERROR: Freesurfer scripts not functioning properly from '...
        'within MATLAB. Check that Freesurfer is properly installed and sourced.'];
end

wbCommands = {'wb_command'};
cmdError = 0;
for c=1:length(wbCommands)
    [status,~] = system([cmdPrefix wbCommands{c}]);
    if status==127
        cmdError=1;
        break;
    end
end
if cmdError==1
    configError = ['ERROR: Connectome Workbench commands are not functioning properly from within '...
        'MATLAB. Check that Connectome Workbencth is properly installed and sourced.  You may '...
        'need to run MATLAB from a terminal on some systems.'];
end

if ~isempty(configError)
    error(configError);
end

end