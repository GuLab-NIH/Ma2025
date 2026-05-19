function [rollSlow,rollLick,rollLickB] = rollingBehavior(logPath,rollRuns)
%% binnedBehavior.m
% calculate velocity, acceleration, licking, and slowing as a funciton of
% track position
%
% Inputs:
%   logPath - behavior log
%   rollRuns - number of runs to use for rolling average
%


%% Process inputs

% load behavior data
lickThresh = 0.5;
logData = readLog(logPath);
params = logParams(logData,lickThresh);

% define slowing settings
binSlow = 70;
postRewSlow = 30;
nRuns = 0;
trStart = 0;
trEnd = round(max(params.y));

% define licking parameters
preRewLick=20;
postRewLick = 30;

% calculate predictive slowing
SpeedChangesPercentile = speedChange_Percentile_optoRevision(logPath,binSlow,postRewSlow,nRuns,trStart,trEnd);
predSlow = SpeedChangesPercentile.percentile;

% get run numbers
nRuns = length(predSlow);
nLoop = nRuns-rollRuns+1;


%% Calculate rolling behavior

% initialize output arrays
rollSlow = zeros(nLoop,1);
rollLick = zeros(nLoop,1);
rollLickB = zeros(nLoop,1);

for ii = 1:nLoop
    % define current run range
    runRng = ii:ii+rollRuns-1;

    % calculate rolling slowing
    rollSlow(ii) = mean(predSlow(runRng),'omitnan');

    % calculate predictive licking for current runs
    curLick = predLicking_allRuns_allLicking_setLaps(params.y,params.lickIdx,...
        params.rewardIdx,preRewLick,postRewLick,runRng,0);
    rollLick(ii) = curLick.perPredInAll;
    rollLickB(ii) = curLick.perPredLick;
end

end

