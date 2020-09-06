
% Function to take 3-column regressor format and generate time series,
% corresponding to boxcar function convolved with HRF.

function regressor = constructRegressor(taskRegr,hrfType,numVols,tr,upsampledTR,subtractHalfTR)

exptDuration = tr*numVols;
restRegr = fpp.func.analysis.defineRestBlocks(taskRegr,exptDuration,0);

allRegr = [taskRegr; restRegr];
[~,indSort] = sort(allRegr(:,1));
allRegr = allRegr(indSort,:);
onsets = allRegr(:,1);
durs = allRegr(:,2);
seq = allRegr(:,3);

if subtractHalfTR
    onsets = onsets-tr/2;   % Note: this code makes first onset -tr/2
end

hrf = constructHRF(upsampledTR,hrfType);

% Upsample sequence vector
onsetBuffer = tr/upsampledTR;    % Extra time points at start of sequence, to avoid negative onsets from subtractHalfTR
totalTime = onsets(end)+durs(end)+10;
sr = 1/upsampledTR; % Sampling rate
nSmps = round(totalTime*sr + onsetBuffer + 1);
seqUpsampled = zeros(nSmps,1);
for j = 1:length(seq)
    if seq(j)~=0
        inds = round(onsets(j)*sr + (1:(durs(j)*sr)) + onsetBuffer);
        seqUpsampled(inds) = seq(j);
    end
end

seqUpsampled = seqUpsampled(onsetBuffer+1:end);
trRatio = tr/upsampledTR;
x = conv(seqUpsampled,hrf);
regressor = x(1:trRatio:numVols*trRatio)/sum(hrf);
if hrfType==3, regressor = -regressor; end  % Invert sign for MION imaging

end