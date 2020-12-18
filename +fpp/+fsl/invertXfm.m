
% Wrapper for FSL's convert_xfm -concat, to combine affine transforms.
% Generates .json metadata for output transform.
%
% fpp.fsl.invertXfm(inputXfm,outputXfm)

function invertXfm(inputXfm,outputXfm)

cmd = ['convert_xfm -omat ' outputXfm ' -inverse ' inputXfm];

% Run convert_xfm command
fpp.util.system(cmd);

% Write json output files
if exist(fpp.bids.jsonPath(inputXfm),'file')
    fpp.bids.jsonReconstruct(inputXfm,outputXfm);
    jsonData = fpp.bids.getMetadata(fpp.bids.jsonPath(inputXfm));
    fpp.bids.jsonChangeValue(outputXfm,{'FromFile','ToFile','CommandLine'},...
        {jsonData.ToFile,jsonData.FromFile,cmd});
end

end