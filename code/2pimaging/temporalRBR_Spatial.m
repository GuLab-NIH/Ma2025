function [rollRBRMean,rollRBR] = temporalRBR_Spatial(corrInfoLocation,rollRuns)
%% temporalRBR
% Calculates a temporal (by run) rolling avergae of run-by-run consistency
% by location
%
% Inputs:
%   corrInfoLocation - RBR correlation structure
%   rollRuns - number of runs to use for rolling average
%


%% Extract data

% get run numbers
nBins = size(corrInfoLocation.toOthers{1},2);
nPairs = size(corrInfoLocation.toOthers{1},1);
nRuns = (1+sqrt(1+8*nPairs))/2;
nLoop = nRuns-rollRuns+1;

% check run number
if nPairs==0 || isempty(nPairs)
    error('Check code')
end

% define toOthers RBR
RBR = corrInfoLocation.toOthers;
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
rollRBR = zeros(nLoop,nBins,nCells);

for ii = 1:nLoop
    % define current run range
    runRng = ii:ii+rollRuns-1;

    % define valid pairs
    pairsUse = all(ismember(runPairs,runRng),2);

    for cc = 1:nCells
        % calculate rolling RBR
        rollRBR(ii,:,cc) = mean(RBR{cc}(pairsUse,:),1,'omitnan');
    end
end

% calculate mean rolling RBR
rollRBRMean = mean(rollRBR,3,'omitnan');

end

