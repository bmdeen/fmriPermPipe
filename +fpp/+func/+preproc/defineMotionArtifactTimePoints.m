
% Function to identify motionParams-based artifact time points, based on
% intervolume translation and rotation.

function artifactTPs = defineMotionArtifactTimePoints(motionParams,transCutoff,rotCutoff,tptsAfter)

% Artifact detection parameters
if ~exist('transCutoff','var')
    transCutoff = .5;               % Artifact detection cutoff: total translation (mm)
end
if ~exist('rotCutoff','var')
    rotCutoff = .5;                 % Artifact detection cutoff: total rotation (degrees)
end
if ~exist('tptsCutoff','var')
    tptsAfter = 0;                  % Remove this many time points after pairs of volumes with motionParams
end

artifactTPs = [];

moDiff = zeros(size(motionParams));
moDiff(2:end,:) = diff(motionParams);
trans = sqrt(sum(moDiff(:,4:6).^2,2));
rot = acos((cos(moDiff(:,1)).*cos(moDiff(:,2)) + cos(moDiff(:,1)).*cos(moDiff(:,3)) + ...
    cos(moDiff(:,2)).*cos(moDiff(:,3)) + sin(moDiff(:,1)).*sin(moDiff(:,2)).*sin(moDiff(:,3)) - 1)/2)*180/pi;
% Remove high-movement time points
artifactTPs = find(trans>transCutoff | rot>rotCutoff);
artifactTPs = sort(union(artifactTPs,artifactTPs-1));   % Include time points before and after a movement
% Remove time points after motionParams volumes, if desired
for j = 1:tptsAfter
    artifactTPs = setdiff(sort(union(artifactTPs,artifactTPs+1)),size(motionParams,1)+1);
end

end