
function argVal = optInputs(argList,argName)

%%% Find intended value for argument argName
argVal = [];
for i = 1:length(argList)
    if strcmp(argList{i},argName) && length(argList)>i
        argVal = argList{i+1};
        break;
    end
end

if isempty(argVal), return; end

%%% Check input compatibility

argGood = 0;

% String
if ismember(argName,{'expt','inputSuffix','outputSuffix','funcTemplatePath',...
        'highresName','maskName','t2StarPath','r2StarPath','outputPath',...
        'roiPathWM','roiPathCSF','motionParPath','badVolPath','disdaqPath',...
        'unsmoothedDataPath','outputMat','maskPath','anatRegDir','funcRegDir',...
        'fsDir','roiDir','funcRoiDir','anatRoiDir','targetName'})
    if ischar(argVal)
        argGood = 1;
    end
    
% 0 or 1 (integer or logical)
elseif ismember(argName,{'overwrite','clustCorrect','permuteRest',...
        'subtractHalfTR','useNuisLinear','useNuisMotion','useNuisWMPCA',...
        'useNuisCSFPCA','useNuisGMR','tempFilt','forceTR','genTarget',...
        'genTSNR','multiEcho','plotResults','writeResiduals','erodeROI',...
        'pcaUnsmoothed','removeBadVols','useSTC','useTedana','undistort'})
    if (isnumeric(argVal) || islogical(argVal)) && isscalar(argVal) && ismember(argVal,[0 1])
        argGood = 1;
    end
    
% Scalar integer in [1,2]
elseif ismember(argName,{'hrfType'})
    if isnumeric(argVal) && isscalar(argVal) && ismember(argVal,[1 2])
        argGood = 1;
    end
    
% Scalar integer in [0,1,2]
elseif ismember(argName,{'nuisanceType'})
    if isnumeric(argVal) && isscalar(argVal) && ismember(argVal,0:2)
        argGood = 1;
    end
    
% Scalar integer in [0,Inf)
elseif ismember(argName,{'tptsAfter','disdaqs','looRun','pcaOrder','echoForMoCorr'})
    if isnumeric(argVal) && isscalar(argVal) && mod(argVal,1)==0 && argVal>=0
        argGood = 1;
    end
    
% Scalar integer in (0,Inf)
elseif ismember(argName,{'filtOrder','permIters'})
    if isnumeric(argVal) && isscalar(argVal) && mod(argVal,1)==0 && argVal>0
        argGood = 1;
    end
    
% Scalar integer in (0,Inf]
elseif ismember(argName,{'echoesToUse'})
    if isnumeric(argVal) && isscalar(argVal) && (mod(argVal,1)==0||argVal==inf) && argVal>0
        argGood = 1;
    end
    
% Vector of integers in (0,Inf)
elseif ismember(argName,{'runList','conList'})
    if isnumeric(argVal) && isvector(argVal) && sum(mod(argVal,1)==0)==length(argVal) && sum(argVal==0)==0
        argGood = 1;
    end
    
% Scalar in (0,1)
elseif ismember(argName,{'voxThresh','clustThresh','faValue'})
    if isnumeric(argVal) && isscalar(argVal) && argVal>0 && argVal<1
        argGood = 1;
    end
    
% Scalar in (0,Inf)
elseif ismember(argName,{'defaultTR','upsampledTR','smThresh'})
    if isnumeric(argVal) && isscalar(argVal) && argVal>0 && argVal~=inf
        argGood = 1;
    end
    
% Scalar in (0,Inf]
elseif ismember(argName,{'transCutoff','rotCutoff','transSingleAxisCutoff',...
        'rotSingleAxisCutoff','stdCutoff'})
    if isnumeric(argVal) && isscalar(argVal) && argVal>0
        argGood = 1;
    end
    
% Scalar in [0,Inf)
elseif ismember(argName,{'fwhm'})
    if isnumeric(argVal) && isscalar(argVal) && argVal>=0 && argVal~=inf
        argGood = 1;
    end
    
% Vector of values in (0,Inf)
elseif ismember(argName,{'teVals'})
    if isnumeric(argVal) && isvector(argVal) && sum(argVal>0)==length(argVal) && sum(argVal==inf)==0
        argGood = 1;
    end
    
% 1- or 2- element vector of values in (0,Inf)
elseif ismember(argName,{'filtCutoff'})
    if isnumeric(argVal) && isvector(argVal) && sum(argVal>0)==length(argVal) && ...
            sum(argVal==inf)==0 && ismember(length(argVal),[1 2])
        argGood = 1;
    end

% Numeric vector or matrix (2D)
elseif ismember(argName,{'customNuisRegr'})
    if isnumeric(argVal) && ~isscalar(argVal) && numel(size(argVal))==2
        argGood = 1;
    end
end

% Cell array of strings
elseif ismember(argName,{'spinEchoPaths','spinEchoPhaseEncodeDirections'})
    if iscellstr(argVal)
        argGood = 1;
    end
end

if argGood==0, argVal = []; end

end