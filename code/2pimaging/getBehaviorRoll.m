%% getBehaviorDistributions
% Get behavior distributions for all mice (velocity, acceleration, licking)

clear; close all; clc

p1 = 'D:\4m_optoPaper';
cd(p1)

nLogsTrue = 10;
rollRuns = [3 4 5 6 7 8];
nR = length(rollRuns);
folds = {'ID20220628A','TM221012','TM221019','TM221026'};
nFolds = length(folds);

% load percentile slowing info
load('Z:\papers\2025RewardOptogenetics\4mData\VirmenLogs\allLickRBREachDay.mat','allLickRBREachDay')

for rr = 1:nR
    % initialize data arrays
    dataBeh = struct();
    dataBeh.Slow = cell(nLogsTrue,nFolds);
    dataBeh.Lick = cell(nLogsTrue,nFolds);
    dataBeh.LickB = cell(nLogsTrue,nFolds);

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
            [rollSlow,rollLick,rollLickB] = rollingBehavior(logs(ii).name,rollRuns(rr));

            % calculate perctile predictive licking
            rollLickC = [movmean(allLickRBREachDay{ff,ii},rollRuns(rr),'omitnan','Endpoints','discard') NaN]';

            % store data
            dataBeh.Slow{ii,ff} = rollSlow;
            dataBeh.Lick{ii,ff} = rollLick;
            dataBeh.LickB{ii,ff} = rollLickB;
            dataBeh.LickC{ii,ff} = rollLickC;
        end
    end

    cd(p1)
    save(['RevisionAnalysis\data\behRolling_' num2str(rollRuns(rr)) '.mat'],'dataBeh')
end
