function [distrVel,distrAcc,distrLick,distrRuns] = binnedBehavior_vRPR(logPath)
%% binnedBehavior.m
% calculate velocity, acceleration, licking, and slowing as a funciton of
% track position

% set parameters
lickThresh = 0.5;
binSize = 5;
trackLength = 400;
binsX = 0:binSize:trackLength;
nBins = length(binsX)-1;

% load virmen log
logData = readLog(logPath);

% read behavior parameters
T = (logData(:,1)-logData(1,1))*24*60*60;
Y = logData(:,3);
L = logData(:,9)>=lickThresh;
R = logData(:,8);

% identify laps
[runIdx,nRuns] = identifyLaps(Y);


%% Calculate behavior distributions

matVel = nan(nRuns,nBins);
matTime = nan(nRuns,nBins);
matLick = nan(nRuns,nBins);

for rr = 1:nRuns
    % find reward index
    idxR = runIdx(rr,1)+find(R(runIdx(rr,1):runIdx(rr,2)))-1;
    if isempty(idxR)
        idxR = runIdx(rr,2);
    end

    curT = T(runIdx(rr,1):idxR-1);
    diffT = movmean([diff(curT); NaN],2);
    curY = smooth(Y(runIdx(rr,1):idxR-1),5);
    curL = L(runIdx(rr,1):idxR-1);

    disp(curY(end))

    % calculate velocity and acceleration
    curVel = [diff(curY)./diff(curT); NaN];

    for bb = 1:nBins
        % define current X
        xCur = curY>=binsX(bb) & curY<=binsX(bb+1);

        if sum(xCur)==0
            continue
        end

        % calculate mean velocity
        matVel(rr,bb) = mean(curVel(xCur),'omitnan');

        % calculate bin time
        matTime(rr,bb) = sum(diffT(xCur));

        % determine lick presence
        if any(curL(xCur)==1)
            matLick(rr,bb) = 1;
        else
            matLick(rr,bb) = 0;
        end
    end
end

% take average time (from center of one bin to the center of the next)
timeMean = (matTime(:,1:end-1)+matTime(:,2:end))/2;

% correct for NaN times
nanTimes = [isnan(+matTime(:,2:end)) nan(size(matTime,1),1)];
timeMean(nanTimes==1) = matTime(nanTimes==1);

% calculate acceleration
matAcc = [diff(matVel,1,2)./timeMean  nan(nRuns,1)];

% calculate distributions
distrVel = mean(matVel,1,'omitnan');
distrAcc = mean(matAcc,1,'omitnan');
distrLick = mean(matLick,1,'omitnan');
distrRuns = sum(~isnan(matLick),1);


end

