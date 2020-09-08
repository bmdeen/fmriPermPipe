
function configError = checkConfig

configError = '';

fslDir = getenv('FSL_DIR');
if isempty(fslDir)
    configError = 'ERROR: FSL_DIR is not defined in bash environment.';
    return;
end

fslCommands = {'fslmaths','fslmerge','fslsplit','fslroi','flirt','fnirt','applywarp',...
    'invwarp','convert_xfm','mcflirt','fslstats','fslval','fsl_tsplot','fslmeants',...
    'bet2','susan','cluster'};
cmdError = 0;
for c=1:length(fslCommands)
    [status,~] = system([fslCommands{c}]);
    if status==127
        cmdError=1;
        break;
    end
end
if cmdError==1
    configError = ['ERROR: FSL scripts not functioning properly from within '...
        'MATLAB. Check that FSL is properly installed and sourced.  You may '...
        'need to run MATLAB from a terminal on some systems.'];
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

end