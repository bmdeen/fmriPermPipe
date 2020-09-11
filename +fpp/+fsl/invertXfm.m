
% Wrapper for FSL's convert_xfm -concat, to combine affine transforms.
% Generates .json metadata for output transform.

function invertXfm(inputXfm,outputXfm)

cmd = ['convert_xfm -omat ' outputXfm ' -inverse ' inputXfm];

% Run convert_xfm command
fpp.util.system(cmd);

% Write json output files
if exist(fpp.bids.jsonPath(inputXfm),'file')
    fpp.bids.jsonReconstruct(inputXfm,outputXfm);
    fpp.bids.jsonChangeValue(outputXfm,{'FromFile','ToFile','CommandLine'},...
        {fpp.bids.getMetadata(inputXfm,'ToFile'),fpp.bids.getMetadata(inputXfm,'FromFile'),cmd});
end

end