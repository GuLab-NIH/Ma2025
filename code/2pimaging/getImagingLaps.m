function imageLapIdx = getImagingLaps()
%% getImagingLaps
% The variables in this code are dependent on the virmen synchronization
% settings of a particular rig and may need to ba adjusted. This code may
% not interact well with active learning code, although extraction B is an
% attempt to resolve this in some situations.
%

% load voltage recording
load('m.mat','m')

%%

% define offset
offSetA = 0.3;
scaleA = 100;
offSetB = 0;
scaleB = 20;
lapRoll = 50;

% identify synchronized points
thresh = 4.8;
a = m(:,2)>thresh;
b = contiguous(a,1);
c = b{1,2};
middleFrameIdx=ceil(mean(c,2));

% split galvo trace
md = m(middleFrameIdx,:);

%%
% extract data from md
y = md(:,3)*100;
lap = md(:,4)*scaleA+offSetA;

% perform smoothing of lap
lapSmooth = movmean(lap,lapRoll);

% find lap start and stop
startLap = round(min(lapSmooth))+1;
stopLap = round(max(lapSmooth))+1;

% check lap range
if startLap==stopLap
    % repeat lap extraction on alternate channel
    lap = md(:,7)*scaleB+offSetB;

    % perform smoothing of lap
    lapSmooth = movmean(lap,lapRoll);

    % find lap start and stop
    startLap = round(min(lapSmooth))+1;
    stopLap = round(max(lapSmooth))+1;

    if startLap==stopLap
        error('Lap range is 0')
    end
end


%%

load('RunByRun_dfof\dfofMInterpM.mat','dfofMInterpM')
nLaps = size(dfofMInterpM{1},1);

if stopLap+1-startLap-nLaps==2
    imageLapIdx = startLap+1:stopLap-1;
else
    error('First or last lap was not removed')
end

end

