
% preprocMRI(study,subjects,varargin)
% 
% Preprocesses fMRI data, including including artifact detection/scrubbing,
% motion correction, skull-stripping, spatial smoothing, intensity
% normalization, high-pass temporal filtering (optional), and registration 
% to a functional template image.  Loops through subjects, experiments, and
% runs by default.  Should be run after convertDCM.
% 
% Example usage: preprocMRI('studyName','SUB*','fwhm',4,'plotResults',0)
%
% Critical output files (in *.prep directory):
% - filtered_func_data.nii.gz: preprocessed data
% - prefiltered_func_data_*.nii.gz: partially preprocessed data
% - example_func.nii.gz: middle functional image (used for registration)
% - mean_func.nii.gz: time-averaged functional image
% - mask.nii.gz: brain mask image
% - tsnr_*.nii.gz: temporal SNR maps
% - reg: registration files between target_func and example_func
% - reg_target: various images transformed to target_func space
% - art/badvols: indices of artifact time points removed from data
% - art/art_results.png: image depicting artifact detection results
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
% - expt (string): name of experiment to process.
% - runList (vector of integers in (0,Inf)): list of runs to process.
% - outputSuffix (string): suffix for .prep directory.  Can be used to run
%       the script multiple times with different options.
% - fwhm (scalar in [0,Inf); default=5): FWHM of spatial smoothing kernel 
%       (mm).
% - smThresh (scalar in (0,Inf)): brightness threshold for SUSAN-based
%       smoothing.
% - faValue (scalar in (0,1); default=.3): fractional intensity threshold
%       for BET2-based brain extraction. Smaller value -> larger mask.
% - targetName (string; default='target_func'): name of image in func 
%       directory to register example_func to, without file extension.
% - genTarget (boolean; default=1): whether to generate target_func image
%       from example_func of first run analyzed if it doesn't exist.
% - forceTR (boolean; default=0): whether to force TR to default value, 
%       instead of reading from header.  Should be used only if TR
%       information in image header is corrupted.
% - defaultTR (scalar in (0,Inf); default=2): default TR (s), used if no TR
%       info is found in functional data header, or if forcetr==1.
% - disdaqs = 3 (integer in [0,Inf); default=3): number of disdaq volumes
%       (at beginning of run) to remove.
% - transCutoff (scalar in [0,Inf]; default=.5): artifact detection cutoff
%       for total translation (mm).
% - rotCutoff (scalar in [0,Inf]; default=.5): artifact detection cutoff
%   for total rotation (degrees).
% - transSingleAxisCutoff (scalar in [0,Inf]; default=Inf): artifact 
%       detection cutoff for translation along individual axes (mm).
% - rotSingleAxisCutoff (scalar in [0,Inf]; default=Inf): artifact 
%       detection cutoff for rotation along individual axes (degrees).
% - tptsAfter (integer in [0,Inf); default=0): remove this many time points
%       after pairs of volumes with motion.
% - stdCutoff (scalar in (0,Inf); default=3.5): artifact detection cutoff
%       for mean signal z-score
% - multiEcho (boolean; default=0): whether data is multi-echo
% - teVals (vector of values in (0,Inf)): TE values (ms) for multi-echo 
%       data.
% - echoesToUse (integer in (0,Inf]; default=inf): number of TEs to use for
%       T2* estimation (inf = use all). Better to use fewer than N(TEs) if
%       signal floors at higher TEs.
% - tempFilt (boolean; default=0): whether to perform high-pass temporal
%       filtering.
% - filtCutoff (scalar in (0,Inf); default=.01: highpass filter cutoff (Hz)
% - filtOrder (integer in (0,Inf); default=60): Hanning-window FIR filter
%       order.
% - genTSNR (boolean; default=0): whether to generate TSNR maps.
% - plotResults (boolean; default=1): whether to display result plots.

function preprocMRI(study,subjects,varargin)

addpath([strrep(mfilename('fullpath'),mfilename,'') '/utils']);

% Load/check config variables.
[configError, studyDir, fslPrefix] = checkConfig(study);
if ~isempty(configError)
    fprintf('%s\n',configError);
    return;
end;

% Basic parameters
overwrite = 0;                  % Whether to overwrite output
expt = [];                      % Experiment to analyze.
runList = [];                   % Runs to analyze.
outputSuffix = '';              % New suffix for preproc output

% General preproc parameters
fwhm = 5;                       % FWHM of spatial smoothing kernel (mm)
smThresh = [];                  % SUSAN brightness threshold
faValue = .3;                   % BET2 fractional intensity threshold. Smaller value -> larger mask.
targetName = 'target_func';     % Image in func directory to register to
genTarget = 1;                  % Whether to generate target_func image from example_func if it doesn't exist
forceTR = 0;                    % Whether to force TR to default value, instead of reading from header
defaultTR = 2;                  % Default TR (s), used if no TR info in header, or if forcetr==1
newMedian = 10000;              % Target value for intensity normalization (median of data is set to this)

% Artifact detection parameters
disdaqs = 3;                    % Number of disdaq volumes (at beginning of run) to remove
transCutoff = .5;               % Artifact detection cutoff: total translation (mm)
rotCutoff = .5;                 % Artifact detection cutoff: total rotation (degrees)
transSingleAxisCutoff = inf;    % Artifact detection cutoff: translation along individual axes (mm)
rotSingleAxisCutoff = inf;      % Artifact detection cutoff: rotation along individual axes (degrees)
tptsAfter = 0;                  % Remove this many time points after pairs of volumes with motion
stdCutoff = 3.5;                % Artifact detection cutoff: mean signal z-score

% Multi-echo parameters
multiEcho = 0;                  % Whether data is multi-echo
teVals = [];                    % TE values (ms)
echoesToUse = inf;              % Number of TEs to use for T2* estimation (inf = use all)

% Temporal filtering parameters
tempFilt = 0;                   % Whether to perform temporal filtering
filtCutoff = .01;               % Highpass filter cutoff (Hz)
filtOrder = 60;                 % Hanning-window FIR filter order

% Data generation parameters
genTSNR = 1;                    % Whether to generate TSNR maps
plotResults = 1;                % Whether to display result plots

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
varArgList = {'overwrite','expt','runList','outputSuffix','targetName','fwhm',...
    'tempFilt','genTarget','multiEcho','teVals','echoesToUse','transCutoff',...
    'rotCutoff','transSingleAxisCutoff','rotSingleAxisCutoff','tptsAfter',...
    'stdCutoff','disdaqs','defaultTR','forceTR','genTSNR','plotResults',...
    'faValue','smThresh'};
for i=1:length(varArgList)
    argVal = optInputs(varargin,varArgList{i});
    if ~isempty(argVal)
        eval([varArgList{i} ' = argVal;']);
    end;
end;

if multiEcho
    if isempty(teVals)
        fprintf('%s\n','ERROR: For multi-echo mode, must specify TE values.');
        return;
    end;
    
    TEs=length(teVals);
    for i=1:TEs
        TEsuffices{i} = ['_TE' int2str(i)];
    end;
    
    if size(teVals,1)<size(teVals,2), teVals = teVals'; end;
    if echoesToUse>TEs, echoesToUse = TEs; end;
else
    TEs = 1;
    TEsuffices{1} = '';
end;

for s = 1:length(subjects)
    
    subject = subjects{s};
    subjDir = [studyDir '/' subject];
    funcDir = [subjDir '/func'];
    analDir = [subjDir '/analysis'];
    targetInit = [funcDir '/' targetName '.nii.gz'];
    
    scanlogFiles = regExpDir(subjDir,'.*scanlo.*([^~]$)');
    
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
    end;
    
    if exist(analDir,'dir')==0, mkdir(analDir); end;
    
    for scan = 1:length(scanlogFiles)
        
        fileID = fopen([subjDir '/' scanlogFiles(scan).name]);
        fileData = textscan(fileID,'%d%s%d\n');
        fclose(fileID);
        dataLengths = cellfun(@length,fileData);
        if length(unique(dataLengths))>1 || sum(dataLengths==0)>1
            fprintf('%s\n',['Scanlog is invalid for ' subject ', scan ' int2str(scan) '. Skipping this scan.']);
            continue;
        end;
        [acqs, expts, runs] = deal(fileData{1},fileData{2},fileData{3});
        
        for i = 1:length(acqs)
            
            % Skip run if it's anatomical or DTI.
            if sum(regexpi(expts{i},'anat'))>0 || sum(regexpi(expts{i},'dti'))>0
                continue;
            end;
            
            % Skip run if it's not the specified experiment/run number.
            if ~isempty(expt) && ~strcmpi(expts{i},expt)
                continue;
            end;
            if ~isempty(runList) && ~ismember(runs(i),runList)
                continue;
            end;
            
            runSuffix = ['-' int2str(runs(i))];
            rfdInit = [funcDir '/' expts{i} runSuffix '.nii.gz'];
            if ~exist(rfdInit,'file')
                fprintf('%s\n',['Raw func data does not exist for ' subject ', expt ' expts{i} '-' int2str(runs(i)) '.']);
                continue;
            end;
            
            outputDir = [analDir '/' expts{i} runSuffix outputSuffix '.prep'];
            
            if exist(outputDir,'dir')
                if overwrite, system(['rm -rf ' outputDir]);
                else continue; end;
            end;
            mkdir(outputDir);
            outputMat = [outputDir '/preproc_options.mat'];
            
            % Copy raw functional data
            rfd = [outputDir '/raw_func.nii.gz'];
            system(['cp ' rfdInit ' ' rfd]);
                
            % Check TR from header
            tr = determineTR(rfd,forceTR,defaultTR,fslPrefix);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% STEP 1: Artifact time point detection
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf('%s\n',[subject ' - ' expts{i} runSuffix ': Step 1, artifact identification']);
            
            % Define raw dataset for artifact detection.
            rfdForAD = [outputDir '/raw_func_forAD.nii.gz'];
            if multiEcho    % Use high-contrast middle-echo image for initial MC
                echoForAD = 1+floor((TEs-1)/2);
                [~,vols] = system([fslPrefix 'fslval ' rfd ' dim4']);
                vols = str2num(strtrim(vols));
                system([fslPrefix 'fslsplit ' rfd ' ' outputDir '/split -t']);
                mergeCmd{1} = [fslPrefix 'fslmerge -tr ' rfdForAD];
                for t=1:vols
                    if mod(t,TEs)==echoForAD
                        mergeCmd{1} = [mergeCmd{1} ' ' outputDir '/split' numPad(t-1,4) '.nii.gz'];
                    end;
                end;
                mergeCmd{1} = [mergeCmd{1} ' ' num2str(tr)];
                system(mergeCmd{1});
                system(['rm -rf ' outputDir '/split*']);
            else
                system(['cp ' rfd ' ' rfdForAD]);
            end;
            
            % Initial motion correction, for artifact time point identification
            mcBase = 'raw_func_mcf';
            system([fslPrefix 'mcflirt -in ' rfdForAD ' -out ' outputDir '/' mcBase ...
                ' -mats -plots -rmsrel -rmsabs']);
            mcDir = [outputDir '/art'];
            mkdir(mcDir);
            system(['mv -f ' outputDir '/' mcBase '.mat ' outputDir '/' mcBase '.par '...
                outputDir '/' mcBase '_abs.rms ' outputDir '/' mcBase '_abs_mean.rms '...
                outputDir '/' mcBase '_rel.rms ' outputDir '/' mcBase '_rel_mean.rms ' ...
                mcDir]);
            system(['rm -rf ' outputDir '/' mcBase ' ' outputDir '/' mcBase '.nii.gz']);
            % Write MC plots
            system([fslPrefix 'fsl_tsplot -i ' mcDir '/' mcBase '.par '...
                '-t ''MCFLIRT estimated rotations (radians)'' -u 1 --start=1 --finish=3 '...
                '-a x,y,z -w 640 -h 144 -o ' mcDir '/rot.png']);
            system([fslPrefix 'fsl_tsplot -i ' mcDir '/' mcBase '.par '...
                '-t ''MCFLIRT estimated translations (mm)'' -u 1 --start=4 --finish=6 '...
                '-a x,y,z -w 640 -h 144 -o ' mcDir '/trans.png']);
            
            % Artifact detection
            moPath = [mcDir '/' mcBase '.par'];
            bvPath = [mcDir '/badvols'];
            motion = load(moPath);
            badVols = [];
            moDiff = zeros(size(motion));
            moDiff(2:end,:) = diff(motion);
            trans = sqrt(sum(moDiff(:,4:6).^2,2));
            rot = acos((cos(moDiff(:,1)).*cos(moDiff(:,2)) + cos(moDiff(:,1)).*cos(moDiff(:,3)) + ...
                cos(moDiff(:,2)).*cos(moDiff(:,3)) + sin(moDiff(:,1)).*sin(moDiff(:,2)).*sin(moDiff(:,3)) - 1)/2)*180/pi;
            % Remove high-movement time points
            badVols = find(trans>transCutoff | rot>rotCutoff | moDiff(:,4)>transSingleAxisCutoff | moDiff(:,5)>transSingleAxisCutoff | ...
                moDiff(:,6)>transSingleAxisCutoff | moDiff(:,1)>rotSingleAxisCutoff | moDiff(:,2)>rotSingleAxisCutoff | moDiff(:,3)>rotSingleAxisCutoff);
            badVols = sort(union(badVols,badVols-1));   % Include time points before and after a movement
            % Remove time points after motion volumes, if desired
            for j = 1:tptsAfter
                badVols = setdiff(sort(union(badVols,badVols+1)),size(motion,1)+1);
            end;
            badVols = sort(union(badVols,1:disdaqs));
            % Remove high-stddev time points
            [~,ts] = system([fslPrefix 'fslmeants -i ' rfdForAD]);
            ts = str2num(strtrim(ts));
            ts = zscore(ts(disdaqs+1:end));
            badVols = sort(union(badVols,find(abs(ts)>stdCutoff)+disdaqs));
            % Write badVols file.
            fid = fopen(bvPath,'w+');
            fprintf(fid,'%d\n',badVols);
            fclose(fid);
            
            % Remove bad volumes, and separate data from different TEs
            [~,vols] = system([fslPrefix 'fslval ' rfd ' dim4']);
            vols = str2num(strtrim(vols));
            system([fslPrefix 'fslsplit ' rfd ' ' outputDir '/split -t']);
            for e=1:TEs
                pfdOrig{e} = [outputDir '/prefiltered_func_data' TEsuffices{e} '_orig.nii.gz'];
                mergeCmd{e} = [fslPrefix 'fslmerge -tr ' pfdOrig{e}];
            end;
            for t=1:vols
                tReal = ceil(t/TEs);
                if ~ismember(tReal,badVols)
                    echoNum = mod(t,TEs); if echoNum==0, echoNum=TEs; end;
                    mergeCmd{echoNum} = [mergeCmd{echoNum} ' ' outputDir '/split' numPad(t-1,4) '.nii.gz'];
                end;
            end;
            for e=1:TEs
                mergeCmd{e} = [mergeCmd{e} ' ' num2str(tr)];
                system(mergeCmd{e});
            end;
            system(['rm -rf ' outputDir '/split*']);
            system(['rm -rf ' rfdForAD]);
            
            % Plot/save artifact detection results
            if multiEcho, vols = vols/TEs; end;
            
            % Compute total translation/rotation relative to reference
            transRef = sqrt(sum(motion(disdaqs+1:end,4:6).^2,2));
            rotRef = acos((cos(motion(disdaqs+1:end,1)).*cos(motion(disdaqs+1:end,2)) + cos(motion(disdaqs+1:end,1)).*cos(motion(disdaqs+1:end,3)) + ...
                cos(motion(disdaqs+1:end,2)).*cos(motion(disdaqs+1:end,3)) + sin(motion(disdaqs+1:end,1)).*sin(motion(disdaqs+1:end,2)).*sin(motion(disdaqs+1:end,3)) - 1)/2)*180/pi;
            
            artifactFig = figure;
            if ~plotResults, set(artifactFig,'visible','off'); end;
            
            subplot(3,1,1);
            plot(disdaqs+1:vols,ts);
            ylimVal = ylim;
            if ~isempty(badVols(disdaqs+1:end))
                for t=badVols(disdaqs+1:end)'
                    line([t t],[ylimVal(1) ylimVal(2)],'Color','r');
                end;
            end;
            set(gca,'YLim',ylimVal);
            title([subject ', ' expts{i} runSuffix outputSuffix ': Mean time series (z-score)'],...
                'interpreter','none');
            
            subplot(3,1,2);
            plot(disdaqs+1:vols,transRef);
            ylimVal = ylim;
            if ~isempty(badVols(disdaqs+1:end))
                for t=badVols(disdaqs+1:end)'
                    line([t t],[ylimVal(1) ylimVal(2)],'Color','r');
                end;
            end;
            set(gca,'YLim',ylimVal);
            title('Total translation (mm)','interpreter','none');
            
            subplot(3,1,3);
            plot(disdaqs+1:vols,rotRef);
            ylimVal = ylim;
            if ~isempty(badVols(disdaqs+1:end))
                for t=badVols(disdaqs+1:end)'
                    line([t t],[ylimVal(1) ylimVal(2)],'Color','r');
                end;
            end;
            set(gca,'YLim',ylimVal);
            title('Total rotation (deg)','interpreter','none');
            
            saveas(artifactFig,[mcDir '/art_results.png']);
            if ~plotResults, close(artifactFig); end;
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% STEP 2: Motion correction
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf('%s\n',[subject ' - ' expts{i} runSuffix ': Step 2, motion correction']);
            
            % Compute # of volumes, exfunc #
            [~,volsNew] = system([fslPrefix 'fslval ' pfdOrig{1} ' dim4']);
            volsNew = str2num(strtrim(volsNew));
            exfuncVolNum = ceil(volsNew/2);   % Example func volume #, indexed by 0
            
            % Motion correction
            for e=1:TEs
                mcBase = ['prefiltered_func_data' TEsuffices{e} '_mcf'];
                pfdMCF{e} = [outputDir '/' mcBase '.nii.gz'];
                system([fslPrefix 'mcflirt -in ' pfdOrig{e} ' -out ' outputDir '/' mcBase ...
                    ' -mats -plots -refvol ' int2str(exfuncVolNum) ' -rmsrel -rmsabs']);
                if e==1
                    mcDir = [outputDir '/mc'];
                else
                    mcDir = [outputDir '/mc' int2str(e)];
                end;
                mkdir(mcDir);
                system(['mv -f ' outputDir '/' mcBase '.mat ' outputDir '/' mcBase '.par '...
                    outputDir '/' mcBase '_abs.rms ' outputDir '/' mcBase '_abs_mean.rms '...
                    outputDir '/' mcBase '_rel.rms ' outputDir '/' mcBase '_rel_mean.rms ' ...
                    mcDir]);
                % Write MC plots
                system([fslPrefix 'fsl_tsplot -i ' mcDir '/' mcBase '.par '...
                    '-t ''MCFLIRT estimated rotations (radians)'' -u 1 --start=1 --finish=3 '...
                    '-a x,y,z -w 640 -h 144 -o ' mcDir '/rot.png']);
                system([fslPrefix 'fsl_tsplot -i ' mcDir '/' mcBase '.par '...
                    '-t ''MCFLIRT estimated translations (mm)'' -u 1 --start=4 --finish=6 '...
                    '-a x,y,z -w 640 -h 144 -o ' mcDir '/trans.png']);
            end;
            
            % NOTE: Slice timing correction would go here in FSL, or before
            % motion correction in SPM or AFNI.  Having already removed
            % artifact volumes, it would not be valid to implement slice
            % timing correction with standard tools at this stage.
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% STEP 2.5: T2* calculation and weighting
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if multiEcho
                fprintf('%s\n',[subject ' - ' expts{i} runSuffix ': Step 2.5, T2*-based weighting']);
                
                for e=1:TEs
                    meanFuncSep{e} = [outputDir '/mean_func' TEsuffices{e} '.nii.gz'];
                    system([fslPrefix 'fslmaths ' pfdMCF{e} ' -Tmean ' meanFuncSep{e}]);
                end;
                brainMask = [outputDir '/mask.nii.gz'];     % Temporary mask
                system([fslPrefix 'bet2 ' meanFuncSep{1} ' ' outputDir '/mask -f ' num2str(faValue) ' -m -n']);
                system(['mv ' outputDir '/mask_mask.nii.gz ' brainMask]);
                meanFuncCombined = [outputDir '/mean_func_TEsep.nii.gz'];
                mergeCmd{1} = [fslPrefix 'fslmerge -tr ' meanFuncCombined];
                for e=1:TEs
                    mergeCmd{1} = [mergeCmd{1} ' ' meanFuncSep{e}];
                end;
                mergeCmd{1} = [mergeCmd{1} ' ' num2str(tr)];
                system(mergeCmd{1});
                for e=1:TEs
                    system(['rm -rf ' meanFuncSep{e}]);
                end;
                
                r2StarPath = [outputDir '/EstR2Star.nii.gz'];
                t2StarPath = [outputDir '/EstT2Star.nii.gz'];
                pfdWeighted = [outputDir '/prefiltered_func_data_mcf.nii.gz'];
                
                funcData = MRIread(meanFuncCombined);
                brainMaskData = MRIread(brainMask);
                dims = size(funcData.vol);
                
                logFuncMat = reshape(log(funcData.vol),[prod(dims(1:3)) TEs]);
                logFuncMat = logFuncMat(:,1:echoesToUse);
                X = [ones(echoesToUse,1) teVals(1:echoesToUse)];
                betas = inv(X'*X)*X'*logFuncMat';
                r2Vec = -betas(2,:)';      % Estimated T2* (ms)
                r2Vec(brainMaskData.vol==0) = 0;        % Apply brain mask
                r2Vol = reshape(r2Vec,dims(1:3));
                t2Vec = -betas(2,:)'.^-1;      % Estimated T2* (ms)
                t2Vec(brainMaskData.vol==0) = 0;        % Apply brain mask
                t2Vol = reshape(t2Vec,dims(1:3));
                
                newData = funcData;
                newData.vol = r2Vol;
                MRIwrite(newData,r2StarPath);
                newData.vol = t2Vol;
                MRIwrite(newData,t2StarPath);
                
                % Weight by estimated T2*
                % For cases where estimated T2*<=0, use equal weighting of all three TEs.
                for e=1:TEs
                    funcSepData{e} = MRIread(pfdMCF{e});
                end;
                weightVol = zeros([dims(1:3) TEs]);
                for e=1:TEs
                    weightVol(:,:,:,e) = teVals(e)*(t2Vol>0).*exp(-teVals(e)./t2Vol) + (t2Vol<=0);
                end;
                weightVol = weightVol./repmat(sum(weightVol,4),[1 1 1 TEs]);
                weightedFuncVol = zeros([dims(1:3) volsNew]);
                for e=1:TEs
                    weightedFuncVol = weightedFuncVol + funcSepData{e}.vol.*repmat(weightVol(:,:,:,e),[1 1 1 volsNew]);
                end;
                
                newData.vol = weightedFuncVol;
                MRIwrite(newData,pfdWeighted);
            else
                pfdWeighted = pfdMCF{1};
            end;
            
            % Extract example_func image
            exfunc = [outputDir '/example_func.nii.gz'];
            system([fslPrefix 'fslroi ' pfdWeighted ' ' exfunc ' ' int2str(exfuncVolNum) ' 1']);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% STEP 3: Skull-stripping
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf('%s\n',[subject ' - ' expts{i} runSuffix ': Step 3, skull-stripping']);
            
            % Compute mean functional image
            meanFunc = [outputDir '/mean_func.nii.gz'];
            system([fslPrefix 'fslmaths ' pfdWeighted ' -Tmean ' meanFunc]);
            
            % Skull-strip functionals
            brainMask = [outputDir '/mask.nii.gz'];
            pfdBET = [outputDir '/prefiltered_func_data_bet.nii.gz'];
            system([fslPrefix 'bet2 ' meanFunc ' ' outputDir '/mask -f ' num2str(faValue) ' -m -n']);
            system(['mv ' outputDir '/mask_mask.nii.gz ' brainMask]);
            system([fslPrefix 'fslmaths ' pfdWeighted ' -mul ' brainMask ' ' pfdBET]);
            
            % Tweak mask based on image values, based on FSL FEAT's method.
            [~,threshVals] = system([fslPrefix 'fslstats ' pfdBET ' -p 2 -p 98']);
            threshVals = str2double(strtrim(regexp(strtrim(threshVals),' ','split')));
            system([fslPrefix 'fslmaths ' pfdBET ' -thr ' num2str(threshVals(2)/10) ...
                ' -Tmin -bin ' brainMask ' -odt char']);
            
            % Compute data median for intensity normalization
            [~,funcMedian] = system([fslPrefix 'fslstats ' pfdWeighted ' -k ' brainMask ' -p 50']);
            funcMedian = str2num(strtrim(funcMedian));
            
            system([fslPrefix 'fslmaths ' brainMask ' -dilF ' brainMask]);    % dilate mask
            system([fslPrefix 'fslmaths ' pfdWeighted ' -mas ' brainMask ' ' pfdBET]);   % reapply mask
            system([fslPrefix 'fslmaths ' pfdBET ' -Tmean ' meanFunc]);      % recompute mean func
            
            % Mask example func image, for registration
            exfuncMasked = [outputDir '/example_func_masked.nii.gz'];
            system([fslPrefix 'fslmaths ' exfunc ' -mul ' brainMask  ' ' exfuncMasked]);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% STEP 4: Spatial smoothing
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf('%s\n',[subject ' - ' expts{i} runSuffix ': Step 4, spatial smoothing']);
            
            % Spatial smoothing, based on FSL FEAT's method.
            pfdSmooth = [outputDir '/prefiltered_func_data_smooth.nii.gz'];
            if isempty(smThresh)
                smThresh = .75*(funcMedian-threshVals(1)); % SUSAN brightness threshold
            end;
            smSigma = fwhm/2.355;
            if fwhm>0
                system([fslPrefix 'susan ' pfdBET ' ' num2str(smThresh) ' ' ...
                    num2str(smSigma) ' 3 1 1 ' meanFunc ' ' num2str(smThresh) ' ' ...
                    strrep(pfdSmooth,'.nii.gz','')]);
            else
                system(['cp ' pfdBET ' ' pfdSmooth]);
            end;
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% STEP 5: Intensity normalization
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf('%s\n',[subject ' - ' expts{i} runSuffix ': Step 5, intensity normalization']);
            
            % Intensity normalization (or "grand mean scaling" in SPM): normalize image median.
            ffd = [outputDir '/filtered_func_data.nii.gz'];
            system([fslPrefix 'fslmaths ' pfdSmooth ' -mul ' num2str(newMedian/funcMedian) ' ' ffd]);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% STEP 5.5: Temporal filtering
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if tempFilt
                fprintf('%s\n',[subject ' - ' expts{i} runSuffix ': Step 5.5, temporal filtering']);
                pfdIntNorm = [outputDir '/prefiltered_func_data_intnorm.nii.gz'];
                system(['mv ' ffd ' ' pfdIntNorm]);
                funcData = MRIread(pfdIntNorm);
                brainMaskData = MRIread(brainMask);
                dims = size(funcData.vol);
                brainMaskInd = find(brainMaskData.vol==1);
                
                % Load data, fill missing time points with zeros (similar to zero padding)
                funcMat = reshape(funcData.vol,[prod(dims(1:3)) dims(4)])';
                funcMat = funcMat(:,brainMaskInd);
                funcMatMean = mean(funcMat);
                funcMat = bsxfun(@minus,funcMat,mean(funcMat));
                goodVols = setdiff(1:vols,badVols);
                funcMatPadded = zeros(goodVols(end),size(funcMat,2));
                funcMatPadded(goodVols,:) = funcMat;
                
                % Generate Hanning-window based FIR filter
                [kernel,~] = fir1(filtOrder,2*tr*filtCutoff,'high');
                
                % Filter each voxel's time series
                for v=1:size(funcMat,2)
                    funcMatPadded(:,v) = conv(funcMatPadded(:,v),kernel,'same');
                end;
                
                funcMat = funcMatPadded(goodVols,:);
                funcMat = bsxfun(@plus,funcMat,funcMatMean);
                
                % Write data
                outData = funcData;
                outData.vol = zeros(size(outData.vol));
                for t=1:size(funcMat,1)
                    tmp = zeros(dims(1:3));
                    tmp(brainMaskInd) = funcMat(t,:);
                    outData.vol(:,:,:,t) = tmp;
                end;
                MRIwrite(outData,ffd);
            end;
            
            % Compute mean of filtered func data
            ffdMean = [outputDir '/mean_func_filtered.nii.gz'];
            system([fslPrefix 'fslmaths ' ffd ' -Tmean ' ffdMean]);
            
            % Compute TSNR map
            if genTSNR
                
                inputPaths = {rfd,pfdWeighted,ffd};
                outputPaths = {[outputDir '/tsnr_raw.nii.gz'],[outputDir ...
                    '/tsnr_mcf.nii.gz'],[outputDir '/tsnr_filtered.nii.gz']};
                
                % For ME data, not meaningful to compute TSNR of raw data
                if multiEcho
                    inputPaths = inputPaths(2:3);
                    outputPaths = outputPaths(2:3);
                end;
                
                for i2 = 1:length(inputPaths)
                    system([fslPrefix 'fslmaths ' inputPaths{i2} ' -Tmean ' outputDir '/tmp_tsnr_calc-1']);
                    system([fslPrefix 'fslmaths ' inputPaths{i2} ' -Tstd ' outputDir '/tmp_tsnr_calc-2']);
                    system([fslPrefix 'fslmaths ' outputDir '/tmp_tsnr_calc-1 -div ' outputDir ...
                        '/tmp_tsnr_calc-2 ' outputPaths{i2}]);
                    system(['rm -rf ' outputDir '/tmp_tsnr_calc-*']);
                end;
                
            end;
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% STEP 6: Registration
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf('%s\n',[subject ' - ' expts{i} runSuffix ': Step 6, registration']);
            prepRegDir = [outputDir '/reg'];
            prepRegTargDir = [outputDir '/reg_target'];
            target = [prepRegDir '/target_func.nii.gz'];
            exfuncRegCopy = [prepRegDir '/example_func.nii.gz'];
            brainMask2TargetImg = [prepRegTargDir '/mask.nii.gz'];
            exfunc2Target = [prepRegDir '/example_func2target_func.mat'];
            exfunc2TargetImg = [prepRegDir '/example_func2target_func.nii.gz'];
            target2Exfunc = [prepRegDir '/target_func2example_func.mat'];
            
            if genTarget && ~exist(targetInit,'file')
                system(['cp ' exfunc ' ' targetInit]);
            end;
            
            mkdir(prepRegDir); mkdir(prepRegTargDir);
            system(['cp ' exfunc ' ' exfuncRegCopy]);
            system(['cp ' targetInit ' ' target]);
            system([fslPrefix 'flirt -in ' exfuncRegCopy ' -out ' exfunc2TargetImg ...
                ' -ref ' target ' -omat ' exfunc2Target ' -cost corratio -dof 6' ...
                ' -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -interp trilinear']);
            system([fslPrefix 'convert_xfm -omat ' target2Exfunc ' -inverse ' ...
                exfunc2Target]);
            system([fslPrefix 'flirt -in ' brainMask ' -out ' brainMask2TargetImg ...
                ' -ref ' target ' -applyxfm -init ' exfunc2Target ' -interp nearestneighbour']);
            
            inputPaths = {ffd,ffdMean,pfdBET};
            outputPaths = {[prepRegTargDir '/filtered_func_data.nii.gz'],...
                [prepRegTargDir '/mean_func_filtered.nii.gz'],[prepRegTargDir '/prefiltered_func_data_bet.nii.gz']};
            for p=1:length(inputPaths)
                system([fslPrefix 'flirt -in ' inputPaths{p} ' -out ' outputPaths{p} ...
                    ' -ref ' target ' -applyxfm -init ' exfunc2Target]);
                system([fslPrefix 'fslmaths ' outputPaths{p} ' -mul ' ...
                    brainMask2TargetImg ' ' outputPaths{p}]);
            end;
            
            save(outputMat,'fwhm','smThresh','faValue','tr','disdaqs',...
                'transCutoff','rotCutoff','transSingleAxisCutoff',...
                'rotSingleAxisCutoff','tptsAfter','stdCutoff','multiEcho',...
                'teVals','echoesToUse','tempFilt','filtCutoff','filtOrder');
            
            fprintf('%s\n\n',['Finished preproc for subject ' subject ', expt ' expts{i} runSuffix outputSuffix '.']);
            
        end;
    end;
end;

end