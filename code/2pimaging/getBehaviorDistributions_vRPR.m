%% getBehaviorDistributions
% Get behavior distributions for all mice (velocity, acceleration, licking)

clear; close all; clc

p1 = 'D:\4m_optoPaper';
cd(p1)

nLogsTrue = 10;
nBins = 80;
folds = {'ID20220628A','TM221012','TM221019','TM221026'};
nFolds = length(folds);

% initialize data arrays
dataBeh = struct();
dataBeh.Vel = zeros(nLogsTrue,nFolds,nBins);
dataBeh.Acc = zeros(nLogsTrue,nFolds,nBins);
dataBeh.Lick = zeros(nLogsTrue,nFolds,nBins);
dataBeh.laps = zeros(nLogsTrue,nFolds,nBins);
dataBeh.VelMean = zeros(nLogsTrue,nFolds);
dataBeh.VelMean = zeros(nLogsTrue,nFolds);

for ff = 1:nFolds
    % move to current folder
    cd([p1 '\' folds{ff} '\VirmenLogs' ])

    % find virmen logs
    logs = dir('*T*.txt');
    nLogs = length(logs);
    if nLogs~=nLogsTrue
        error('Invalid number of logs')
    end

    for ii = 1:nLogs
        % calculate distributions
        [distrVel,distrAcc,distrLick,distrRuns] = binnedBehavior_vRPR(logs(ii).name);

        % calculate mean behavior
        [meanVel,meanAcc] = meanBehavior(logs(ii).name);

        % store data
        dataBeh.Vel(ii,ff,:) = distrVel;
        dataBeh.Acc(ii,ff,:) = distrAcc;
        dataBeh.Lick(ii,ff,:) = distrLick;
        dataBeh.laps(ii,ff,:) = distrRuns;
        dataBeh.VelMean(ii,ff) = meanVel;
        dataBeh.AccMean(ii,ff) = meanAcc;
    end
end

cd(p1)
save('RevisionAnalysis\behDistribution_vRPR.mat','dataBeh')

