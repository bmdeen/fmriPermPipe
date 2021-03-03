
% FPP utility to check validity of optional inputs (varargin).
%
% fpp.util.optInputs(argList,argName)

function argVal = optInputs(argList,argName)

%%% Find intended value for argument argName
argVal = [];
for i = 1:length(argList)
    if strcmpi(argList{i},argName) && length(argList)>i
        argVal = argList{i+1};
        break;
    end
end

if isempty(argVal), return; end

%%% Check input compatibility
argGood = 0;

% String
if ismember(lower(argName),lower({'expt','inputSuffix','outputSuffix','funcTemplatePath',...
        'highresName','maskName','t2StarPath','r2StarPath','outputPath',...
        'roiPathWM','roiPathCSF','motionParPath','badVolPath','disdaqPath',...
        'unsmoothedDataPath','outputMat','maskPath','anatRegDir','funcRegDir',...
        'fsDir','roiDir','funcRoiDir','anatRoiDir','targetName','funcDataPhaseEncodeDirection',...
        'funcTemplateSpace','init','datatype','cost','searchcost','interp','sincwindow',...
        'aff','inwarp','intin','cout','iout','fout','jout','logout','config',...
        'refmask','inmask','minmet','miter','subsamp','warpres','infwhm','reffwhm',...
        'regmod','jacrange','intmod','biasres','numprec','interp','warp','postmat',...
        'premat','midmat','mask','warp1','warp2','jacobian','confoundPath','outlierPath',...
        'analysisDir','fieldMapParamPath','topupWarpPath','topupJacobianPath','spinEchoPath'}))
    if ischar(argVal)
        argGood = 1;
    else
        errorMsg = ['Argument ' argName ' must be a string.'];
    end
    
% 0 or 1 (integer or logical)
elseif ismember(lower(argName),lower({'overwrite','clustCorrect','permuteRest',...
        'subtractHalfTR','useNuisLinear','useNuisMotion','useNuisWMPCA',...
        'useNuisCSFPCA','useNuisGMR','tempFilt','genTarget',...
        'genTSNR','multiEcho','plotResults','writeResiduals','erodeROI',...
        'pcaUnsmoothed','removeBadVols','useSTC','useTedana','undistort',...
        'useTaskTemplate','applyxfm','imprefm','impinm','ssqlambda','refderiv',...
        'estint','verbose','abs','rel','absout','relout','super','usesqform',...
        'constrainj','noconstraint','useDespike','nosearch','deleteMidprep'}))
    if (isnumeric(argVal) || islogical(argVal)) && isscalar(argVal) && ismember(argVal,[0 1])
        argGood = 1;
    else
        errorMsg = ['Argument ' argName ' must be a boolean.'];
    end
    
% Scalar integer in [1,2]
elseif ismember(lower(argName),lower({'hrfType'}))
    if isnumeric(argVal) && isscalar(argVal) && ismember(argVal,[1 2])
        argGood = 1;
    else
        errorMsg = ['Argument ' argName ' must be 1 or 2.'];
    end
    
% Scalar integer in [0,1,2]
elseif ismember(lower(argName),lower({'nuisanceType'}))
    if isnumeric(argVal) && isscalar(argVal) && ismember(argVal,0:2)
        argGood = 1;
    else
        errorMsg = ['Argument ' argName ' must be an integer from 0 to 2.'];
    end
    
% Scalar integer in [0,Inf)
elseif ismember(lower(argName),lower({'tptsAfter','disdaqs','looRun','pcaOrder',...
        'superlevel','paddingsize','randSeed','permIters'}))
    if isnumeric(argVal) && isscalar(argVal) && mod(argVal,1)==0 && argVal>=0
        argGood = 1;
    else
        errorMsg = ['Argument ' argName ' must be an integer >= 0.'];
    end
    
% Scalar integer in (0,Inf)
elseif ismember(lower(argName),lower({'filtOrder','dof','splineorder',...
        'intorder','echoForMoCorr'}))
    if isnumeric(argVal) && isscalar(argVal) && mod(argVal,1)==0 && argVal>0
        argGood = 1;
    else
        errorMsg = ['Argument ' argName ' must be an integer > 0.'];
    end
    
% Scalar integer in (0,Inf]
elseif ismember(lower(argName),lower({'echoesToUse'}))
    if isnumeric(argVal) && isscalar(argVal) && (mod(argVal,1)==0||argVal==inf) && argVal>0
        argGood = 1;
    else
        errorMsg = ['Argument ' argName ' must be an integer > 0 or Inf.'];
    end
    
% Vector of integers in (0,Inf)
elseif ismember(lower(argName),lower({'runList','conList','outlierInd'}))
    if isnumeric(argVal) && isvector(argVal) && sum(mod(argVal,1)==0)==length(argVal) && sum(argVal==0)==0
        argGood = 1;
    else
        errorMsg = ['Argument ' argName ' must be a vector of integers > 0.'];
    end
    
% Scalar in (0,1)
elseif ismember(lower(argName),lower({'voxThresh','clustThresh','faValue'}))
    if isnumeric(argVal) && isscalar(argVal) && argVal>0 && argVal<1
        argGood = 1;
    else
        errorMsg = ['Argument ' argName ' must be a scalar in (0,1).'];
    end
    
% Scalar in (0,Inf)
elseif ismember(lower(argName),lower({'upsampledTR','smThresh'}))
    if isnumeric(argVal) && isscalar(argVal) && argVal>0 && argVal~=inf
        argGood = 1;
    else
        errorMsg = ['Argument ' argName ' must be a scalar > 0.'];
    end
    
% Scalar in (0,Inf]
elseif ismember(lower(argName),lower({'transCutoff','rotCutoff','transSingleAxisCutoff',...
        'rotSingleAxisCutoff','stdCutoff'}))
    if isnumeric(argVal) && isscalar(argVal) && argVal>0
        argGood = 1;
    else
        errorMsg = ['Argument ' argName ' must be a scalar > 0 or Inf.'];
    end
    
% Scalar in [0,Inf)
elseif ismember(lower(argName),lower({'fwhm','imprefval','impinval','lambda','biaslambda',...
        'jmin','jmax'}))
    if isnumeric(argVal) && isscalar(argVal) && argVal>=0 && argVal~=inf
        argGood = 1;
    else
        errorMsg = ['Argument ' argName ' must be a scalar >= 0.'];
    end
    
% Vector of values in (0,Inf)
elseif ismember(lower(argName),lower({'teVals','sliceTimes'}))
    if isnumeric(argVal) && isvector(argVal) && sum(argVal>0)==length(argVal) && sum(argVal==inf)==0
        argGood = 1;
    else
        errorMsg = ['Argument ' argName ' must be a vector of numbers > 0.'];
    end
    
% 1- or 2- element vector of values in (0,Inf)
elseif ismember(lower(argName),lower({'filtCutoff'}))
    if isnumeric(argVal) && isvector(argVal) && sum(argVal>0)==length(argVal) && ...
            sum(argVal==inf)==0 && ismember(length(argVal),[1 2])
        argGood = 1;
    else
        errorMsg = ['Argument ' argName ' must be a 1- or 2- element vector of numbers > 0.'];
    end
    
% 2-element vector of values in [0,180]
elseif ismember(lower(argName),lower({'searchrx','searchry','searchrz'}))
    if isnumeric(argVal) && isvector(argVal) && length(argVal)==2 && ...
            sum(argVal<=180)==length(argVal) && sum(argVal==inf)==0
        argGood = 1;
    else
        errorMsg = ['Argument ' argName ' must be a 2-element vector of numbers in [0,180].'];
    end

% Numeric vector or matrix (2D)
elseif ismember(lower(argName),lower({'customNuisRegr'}))
    if isnumeric(argVal) && ~isscalar(argVal) && numel(size(argVal))==2
        argGood = 1;
    else
        errorMsg = ['Argument ' argName ' must be a numeric vector or 2D matrix.'];
    end

% Cell array of strings
elseif ismember(lower(argName),lower({'spinEchoPaths','spinEchoPhaseEncodeDirections',...
        'condNames','confoundNames','contrastNames','phaseEncodeDirections'}))
    if iscellstr(argVal)
        argGood = 1;
    else
        errorMsg = ['Argument ' argName ' must be a cell array of strings.'];
    end
end

%%% Throw error if argVal is not compatible
if argGood==0
    argVal = []; 
    error(errorMsg);
end

end