function [meanVel,meanAcc] = meanBehavior(logPath)
%% binnedBehavior.m
% calculate mean velocity, acceleration

% load virmen log
logData = readLog(logPath);

% read behavior parameters
T = (logData(:,1)-logData(1,1))*24*60*60;
Y = logData(:,3);

% identify laps
[runIdx,nRuns] = identifyLaps(Y);


%% Calculate behavior distributions

matVel = nan(nRuns,1);
matAcc = nan(nRuns,1);
spThresh = 1;
accThresh = 500;

for rr = 1:nRuns
    curT = T(runIdx(rr,1):runIdx(rr,2));
    curY = smooth(Y(runIdx(rr,1):runIdx(rr,2)),5);

    % calculate velocity and acceleration
    curVel = [diff(curY)./diff(curT); NaN];
    curAcc = diff(smooth(curVel,5))./diff(curT);

    % calculate lap mean velocity and acceleration
    meanVelCur = mean(curVel(curVel>spThresh),'omitnan');
    meanAccCur = mean(curAcc(abs(curAcc)<accThresh),'omitnan');

    % store values
    matVel(rr) = meanVelCur;
    matAcc(rr) = meanAccCur;
end

% take across lap averages
meanVel = mean(matVel,'omitnan');
meanAcc = mean(matAcc,1,'omitnan');

end

