
% Function to identify motion-based artifact time points, based on
% framewise translation and rotation.
%
% artifactTPs = fpp.func.preproc.defineMotionArtifactTimePoints(confoundTSV,...
%   fdCutoff,transCutoff,rotCutoff,tptsAfter)

function artifactTPs = defineMotionArtifactTimePoints(confoundTSV,fdCutoff,transCutoff,rotCutoff,tptsAfter)

% Artifact detection parameters
if ~exist('fdCutoff','var')
    fdCutoff = .5;              % Artifact detection cutoff: total translation (mm)
end
if ~exist('transCutoff','var')
    transCutoff = Inf;       	% Artifact detection cutoff: total translation (mm)
end
if ~exist('rotCutoff','var')
    rotCutoff = Inf;           	% Artifact detection cutoff: total rotation (degrees)
end
if ~exist('tptsCutoff','var')
    tptsAfter = 0;            	% Remove this many time points after pairs of volumes with motionParams
end

artifactTPs = [];

% Remove high-movement time points
artifactTPs = find(confoundTSV.framewise_translation>transCutoff | ...
    confoundTSV.framewise_rotation>rotCutoff | confoundTSV.framewise_displacement>fdCutoff);
artifactTPs = sort(union(artifactTPs,artifactTPs-1));   % Include time points before and after a movement
% Remove time points after motionParams volumes, if desired
for j = 1:tptsAfter
    artifactTPs = setdiff(sort(union(artifactTPs,artifactTPs+1)),length(confoundTSV.framewise_translation)+1);
end

end