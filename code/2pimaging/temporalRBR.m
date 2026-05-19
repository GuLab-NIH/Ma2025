function [rollRBRMean,rollRBR] = temporalRBR(corrInfo,rollRuns)
%% temporalRBR
% Calculates a temporal (by run) rolling avergae of run-by-run consistency
%
% Inputs:
%   corrInfo - RBR correlation structure
%   rollRuns - number of runs to use for rolling average
%


%% Extract data

% get run numbers
nRuns = size(corrInfo.noNaN,1);
nLoop = nRuns-rollRuns+1;

% define toOthers RBR
RBR = corrInfo.toOthers;
nPairs = size(RBR,1);
nCells = size(RBR,2);

% calculate run pairs
runPairs = zeros(nPairs,2);
idx = 0;
for ii = 1:nRuns-1
    for jj = ii+1:nRuns
        % skip matched runs
        if ii==jj; continue; end

        % store run pair
        idx = idx+1;
        runPairs(idx,:) = [ii,jj];
    end
end


%% Calculate rolling RBR correlation

% initialize output arrays
rollRBR = zeros(nLoop,nCells);

for ii = 1:nLoop
    % define current run range
    runRng = ii:ii+rollRuns-1;

    % define valid pairs
    pairsUse = all(ismember(runPairs,runRng),2);

    % calculate rolling RBR
    rollRBR(ii,:) = mean(RBR(pairsUse,:),1,'omitnan');
end

% calculate mean rolling RBR
rollRBRMean = mean(rollRBR,2,'omitnan');

end

