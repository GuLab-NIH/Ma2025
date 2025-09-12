%% this one: used 40 cm after stim zone (because there is already an expansion), and the third condition is combining stim zone and exp zone (rather than exp zone alone) 
% 
% 50 cm before after and expanded stim zone, aligning stim zone: more for RS: only take the stim region with at least 100 cm from other loc (the beginning region can have less)
load('allSRBRUseRunLowerLick.mat')
%stim zones of RS data
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\stimZoneIdxAllDataFix.mat');
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\positionSegAllData.mat');
%these data order: CHR2R2, GRPRS, CHR2CS, GRPCS
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\idxCHR2CSLowerLick.mat')
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\idxCHR2RS.mat')
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\idxGFPCS.mat')
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\idxGFPRS.mat')
%THESE idx oder: CHR2CS, CHR2RS, GFPCS, GFPRS
allIdxAll={};
allIdxAll{1}=idxCHR2CSLowerLick;
allIdxAll{2}=idxCHR2RS;
allIdxAll{3}=idxGFPCS;
allIdxAll{4}=idxGFPRS;
%for chr2RS
figure
meanSpeedAll={};

positionDisThresh=99;
NPoints=[];

NSession=[10 12 12 11]; %NUMBER OF SESSIONS IN EACH CONDITION: CHR2 CS, RS, GRP CS RS
for N=1;

for NRunStart=1:10; %run start: include this run
totalDataStart=(NRunStart-1)*(NSession(N))+1;
NRunEnd=NRunStart; %run end: include this run
totalDataEnd=NRunEnd*(NSession(N));

mouseIdx=allIdxAll{N};

position=positionSegAllData(mouseIdx);
stim=stimZoneIdxAllDataFix(mouseIdx);
positionAll={};

stimLoc=[];
stimLoc(1,:)=[63+3.5 63+3.5+15];
stimLoc(2,:)=[163+3.5 163+3.5+15];
stimLoc(3,:)=[279+3.5 279+3.5+15];

% stimLoc(1,:)=[63 63+3.5+15];
% stimLoc(2,:)=[163 163+3.5+15];
% stimLoc(3,:)=[279 279+3.5+15];

for n=1:length(position); %per session
    for m=1:length(position{n}) %per run
        positionAll{n,m}=stimLoc; %all use the same stim loc
    end
end



startBinStim={};%each row is a session/run: %run 1: SESSION 1-N,, then run 2 for session1:N... the three numbers are the start point of the three stim
 for run=1:size(positionAll,2);
for session=1:size(positionAll,1)
   
        if ~isempty(positionAll{session,run})

        startPoint=positionAll{session,run}(:,1);
        startBinStim{end+1}=ceil(startPoint);
        else
            startBinStim{end+1}=[];

    end
    end
end

%reorganize the RBR speed
RBRSpeed=[]; %run 1: SESSION 1-N,, then run 2 for session1:N..
S=allSPostRBRUseRun{N};
for run=1:length(S); %EACH RUN
    for session=1:length(mouseIdx); %EACH SESSION
        RBRSpeed(end+1,:)=S{run}(session,:);
    end
end

%use speed of all stim
%remove the first 15cm: mouse lick for reward: so the start should be
%minimal 25cm. remove close to reward region, so the start should be max
%310 (mouse start to slow down) -15 (this is after zone) - 15 (this is the
%reward zone) = 280;
%so only use the stim zone that started between 25 and 280cm
%
% stimStart=25;
% stimEnd=280;
exp=20;
beforeAfterZone=50; %first 15min: many sessions show a large deceleration
stimStart=150;
% stimEnd=310-beforeAfterZone-15-exp+1;
stimEnd=250;

speed=[];
% for n=1:length(startBinStim)
%first NRun


for n=totalDataStart:totalDataEnd;%
    for m=1:length(startBinStim{n})
        startThis=startBinStim{n}(m);
        if startThis-stimStart>0 & startThis-stimEnd<0;
            speed(end+1,:)=RBRSpeed(n,startThis-beforeAfterZone:startThis+15+beforeAfterZone-1);
        end
    end
end



zones={};
zones{1}=[1:1:(beforeAfterZone-exp)];
zones{2}=[beforeAfterZone-exp+1:1:beforeAfterZone];
zones{3}=[beforeAfterZone+1:1:beforeAfterZone+15];
zones{4}=[beforeAfterZone+15+1:1:beforeAfterZone+15+exp];
zones{5}=[beforeAfterZone+15+exp+1:1:size(speed,2)];
meanSpeed=[];
meanSpeed=[nanmean(speed(:,zones{1}),2) nanmean(speed(:,zones{2}),2) nanmean(speed(:,zones{3}),2) nanmean(speed(:,zones{4}),2) nanmean(speed(:,zones{5}),2)];
NPoints(N,1)=size(meanSpeed,1);
meanSpeedAll{N}{NRunStart}=meanSpeed;


%normalize
% speed=speed./meanSpeed(:,1);
subplot(10,2,2*(NRunStart-1)+1)
semshade(speed,0.3,'g',[1:1:size(speed,2)]);
M=nanmean(speed,1);
E=nansem(speed,1);
y2=max(M)+max(E);
y1=min(M)-min(E);
hold on
line([beforeAfterZone+1 beforeAfterZone+1],[y1 y2])
hold on
line([beforeAfterZone+15 beforeAfterZone+15],[y1 y2])
hold on
line([beforeAfterZone+15+exp beforeAfterZone+15+exp],[y1 y2])
hold on
line([beforeAfterZone-exp beforeAfterZone-exp],[y1 y2])

% if N==1 | N==2;
% ylim([25 44])
% else
%     ylim([28 44]);
% end



subplot(10,2,NRunStart*2);
bar([1:1:5],nanmean(meanSpeed,1),'FaceColor','w');
hold on
errorbar([1:1:5],nanmean(meanSpeed,1),nansem(meanSpeed,1),'k')
ylim([0 45])
[~,p1]=ttest(meanSpeed(:,1),meanSpeed(:,2));
[~,p2]=ttest(meanSpeed(:,1),meanSpeed(:,3));
[~,p3]=ttest(meanSpeed(:,1),meanSpeed(:,4));
[~,p4]=ttest(meanSpeed(:,1),meanSpeed(:,5));
title(num2str([p1 p2 p3 p4]))
end
end
after=meanSpeedAll{N};
tightfig;
saveas(gcf,'aligniningStimZoneLowerLickFarApart50cmExp_combineStimDelay5PostEarly_stim2.fig')
print -painters -depsc aligniningStimZoneLowerLickFarApart50cmExp_combineStimDelay5PostEarly_stim2.eps

%% last run of stim: this one: used 40 cm after stim zone (because there is already an expansion), and the third condition is combining stim zone and exp zone (rather than exp zone alone) 
% 
% 50 cm before after and expanded stim zone, aligning stim zone: more for RS: only take the stim region with at least 100 cm from other loc (the beginning region can have less)
load('allSRBRUseRunLowerLick.mat')
%stim zones of RS data
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\stimZoneIdxAllDataFix.mat');
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\positionSegAllData.mat');
%these data order: CHR2R2, GRPRS, CHR2CS, GRPCS
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\idxCHR2CSLowerLick.mat')
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\idxCHR2RS.mat')
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\idxGFPCS.mat')
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\idxGFPRS.mat')
%THESE idx oder: CHR2CS, CHR2RS, GFPCS, GFPRS
allIdxAll={};
allIdxAll{1}=idxCHR2CSLowerLick;
allIdxAll{2}=idxCHR2RS;
allIdxAll{3}=idxGFPCS;
allIdxAll{4}=idxGFPRS;
%for chr2RS
figure
meanSpeedAll={};

positionDisThresh=99;
NPoints=[];

NSession=[10 12 12 11];
for N=1;
    for NRunStart=1:10; %run start: include this run
totalDataStart=(NRunStart-1)*(NSession(N))+1;
NRunEnd=NRunStart; %run end: include this run
totalDataEnd=NRunEnd*(NSession(N));

mouseIdx=allIdxAll{N};

position=positionSegAllData(mouseIdx);
stim=stimZoneIdxAllDataFix(mouseIdx);
positionAll={};

stimLoc=[];
stimLoc(1,:)=[63+3.5 63+3.5+15];
stimLoc(2,:)=[163+3.5 163+3.5+15];
stimLoc(3,:)=[279+3.5 279+3.5+15];
% stimLoc(1,:)=[63 63+3.5+15];
% stimLoc(2,:)=[163 163+3.5+15];
% stimLoc(3,:)=[279 279+3.5+15];


for n=1:length(position); %per session
    for m=1:length(position{n}) %per run
        positionAll{n,m}=stimLoc; %all use the same stim loc
    end
end



startBinStim={};%each row is a session/run: %run 1: SESSION 1-N,, then run 2 for session1:N... the three numbers are the start point of the three stim
 for run=1:size(positionAll,2);
for session=1:size(positionAll,1)
   
        if ~isempty(positionAll{session,run})

        startPoint=positionAll{session,run}(:,1);
        startBinStim{end+1}=ceil(startPoint);
        else
            startBinStim{end+1}=[];

    end
    end
end

%reorganize the RBR speed
RBRSpeed=[]; %run 1: SESSION 1-N,, then run 2 for session1:N..
S=allSDurRBRUseRun{N};
for run=1:length(S); %EACH RUN
    for session=1:length(mouseIdx); %EACH SESSION
        RBRSpeed(end+1,:)=S{run}(session,:);
    end
end

%use speed of all stim
%remove the first 15cm: mouse lick for reward: so the start should be
%minimal 25cm. remove close to reward region, so the start should be max
%310 (mouse start to slow down) -15 (this is after zone) - 15 (this is the
%reward zone) = 280;
%so only use the stim zone that started between 25 and 280cm
%
% stimStart=25;
% stimEnd=280;
exp=20;;
beforeAfterZone=50; %first 15min: many sessions show a large deceleration
stimStart=150;
% stimEnd=310-beforeAfterZone-15-exp+1;
stimEnd=250;

speed=[];
% for n=1:length(startBinStim)


for n=totalDataStart:totalDataEnd;%
    for m=1:length(startBinStim{n})
        startThis=startBinStim{n}(m);
        if startThis-stimStart>0 & startThis-stimEnd<0;
            speed(end+1,:)=RBRSpeed(n,startThis-beforeAfterZone:startThis+15+beforeAfterZone-1);
        end
    end
end



zones={};
zones{1}=[1:1:(beforeAfterZone-exp)];
zones{2}=[beforeAfterZone-exp+1:1:beforeAfterZone];
zones{3}=[beforeAfterZone+1:1:beforeAfterZone+15];
zones{4}=[beforeAfterZone+15+1:1:beforeAfterZone+15+exp];
zones{5}=[beforeAfterZone+15+exp+1:1:size(speed,2)];
meanSpeed=[];
meanSpeed=[nanmean(speed(:,zones{1}),2) nanmean(speed(:,zones{2}),2) nanmean(speed(:,zones{3}),2) nanmean(speed(:,zones{4}),2) nanmean(speed(:,zones{5}),2)];
NPoints(N,1)=size(meanSpeed,1);
meanSpeedAll{N}{NRunStart}=meanSpeed;


%normalize
% speed=speed./meanSpeed(:,1);
subplot(10,2,2*(NRunStart-1)+1)
semshade(speed,0.3,'g',[1:1:size(speed,2)]);
M=nanmean(speed,1);
E=nansem(speed,1);
y2=max(M)+max(E);
y1=min(M)-min(E);
hold on
line([beforeAfterZone+1 beforeAfterZone+1],[y1 y2])
hold on
line([beforeAfterZone+15 beforeAfterZone+15],[y1 y2])
hold on
line([beforeAfterZone+15+exp beforeAfterZone+15+exp],[y1 y2])
hold on
line([beforeAfterZone-exp beforeAfterZone-exp],[y1 y2])

% if N==1 | N==2;
% ylim([25 44])
% else
%     ylim([28 44]);
% end



subplot(10,2,NRunStart*2);
bar([1:1:5],nanmean(meanSpeed,1),'FaceColor','w');
hold on
errorbar([1:1:5],nanmean(meanSpeed,1),nansem(meanSpeed,1),'k')
ylim([0 45])
[~,p1]=ttest(meanSpeed(:,1),meanSpeed(:,2));
[~,p2]=ttest(meanSpeed(:,1),meanSpeed(:,3));
[~,p3]=ttest(meanSpeed(:,1),meanSpeed(:,4));
[~,p4]=ttest(meanSpeed(:,1),meanSpeed(:,5));
title(num2str([p1 p2 p3 p4]))
    end
end
dur=meanSpeedAll{N};
tightfig;
saveas(gcf,'aligniningStimZoneLowerLickFarApart50cmExp_combineStimDelay5DurLast_stim2.fig')
print -painters -depsc aligniningStimZoneLowerLickFarApart50cmExp_combineStimDelay5DurLast_stim2.eps

%% this one: used 40 cm after stim zone (because there is already an expansion), and the third condition is combining stim zone and exp zone (rather than exp zone alone) 
% 
% 50 cm before after and expanded stim zone, aligning stim zone: more for RS: only take the stim region with at least 100 cm from other loc (the beginning region can have less)
load('allSRBRUseRunLowerLick.mat')
%stim zones of RS data
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\stimZoneIdxAllDataFix.mat');
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\positionSegAllData.mat');
%these data order: CHR2R2, GRPRS, CHR2CS, GRPCS
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\idxCHR2CSLowerLick.mat')
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\idxCHR2RS.mat')
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\idxGFPCS.mat')
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\idxGFPRS.mat')
%THESE idx oder: CHR2CS, CHR2RS, GFPCS, GFPRS
allIdxAll={};
allIdxAll{1}=idxCHR2CSLowerLick;
allIdxAll{2}=idxCHR2RS;
allIdxAll{3}=idxGFPCS;
allIdxAll{4}=idxGFPRS;
%for chr2RS
figure
meanSpeedAll={};

positionDisThresh=99;
NPoints=[];

NSession=[10 12 12 11]; %NUMBER OF SESSIONS IN EACH CONDITION: CHR2 CS, RS, GRP CS RS
for N=1;

for NRunStart=1:10; %run start: include this run
totalDataStart=(NRunStart-1)*(NSession(N))+1;
NRunEnd=NRunStart; %run end: include this run
totalDataEnd=NRunEnd*(NSession(N));

mouseIdx=allIdxAll{N};

position=positionSegAllData(mouseIdx);
stim=stimZoneIdxAllDataFix(mouseIdx);
positionAll={};

stimLoc=[];
stimLoc(1,:)=[63+3.5 63+3.5+15];
stimLoc(2,:)=[163+3.5 163+3.5+15];
stimLoc(3,:)=[279+3.5 279+3.5+15];

% stimLoc(1,:)=[63 63+3.5+15];
% stimLoc(2,:)=[163 163+3.5+15];
% stimLoc(3,:)=[279 279+3.5+15];


for n=1:length(position); %per session
    for m=1:length(position{n}) %per run
        positionAll{n,m}=stimLoc; %all use the same stim loc
    end
end



startBinStim={};%each row is a session/run: %run 1: SESSION 1-N,, then run 2 for session1:N... the three numbers are the start point of the three stim
 for run=1:size(positionAll,2);
for session=1:size(positionAll,1)
   
        if ~isempty(positionAll{session,run})

        startPoint=positionAll{session,run}(:,1);
        startBinStim{end+1}=ceil(startPoint);
        else
            startBinStim{end+1}=[];

    end
    end
end

%reorganize the RBR speed
RBRSpeed=[]; %run 1: SESSION 1-N,, then run 2 for session1:N..
S=allSPreRBRUseRun{N};
for run=1:length(S); %EACH RUN
    for session=1:length(mouseIdx); %EACH SESSION
        RBRSpeed(end+1,:)=S{run}(session,:);
    end
end

%use speed of all stim
%remove the first 15cm: mouse lick for reward: so the start should be
%minimal 25cm. remove close to reward region, so the start should be max
%310 (mouse start to slow down) -15 (this is after zone) - 15 (this is the
%reward zone) = 280;
%so only use the stim zone that started between 25 and 280cm
%
% stimStart=25;
% stimEnd=280;
exp=20;;
beforeAfterZone=50; %first 15min: many sessions show a large deceleration
stimStart=150;
% stimEnd=310-beforeAfterZone-15-exp+1;
stimEnd=250;

speed=[];
% for n=1:length(startBinStim)
%first NRun


for n=totalDataStart:totalDataEnd;%
    for m=1:length(startBinStim{n})
        startThis=startBinStim{n}(m);
        if startThis-stimStart>0 & startThis-stimEnd<0;
            speed(end+1,:)=RBRSpeed(n,startThis-beforeAfterZone:startThis+15+beforeAfterZone-1);
        end
    end
end



zones={};
zones{1}=[1:1:(beforeAfterZone-exp)];
zones{2}=[beforeAfterZone-exp+1:1:beforeAfterZone];
zones{3}=[beforeAfterZone+1:1:beforeAfterZone+15];
zones{4}=[beforeAfterZone+15+1:1:beforeAfterZone+15+exp];
zones{5}=[beforeAfterZone+15+exp+1:1:size(speed,2)];
meanSpeed=[];
meanSpeed=[nanmean(speed(:,zones{1}),2) nanmean(speed(:,zones{2}),2) nanmean(speed(:,zones{3}),2) nanmean(speed(:,zones{4}),2) nanmean(speed(:,zones{5}),2)];
NPoints(N,1)=size(meanSpeed,1);
meanSpeedAll{N}{NRunStart}=meanSpeed;


%normalize
% speed=speed./meanSpeed(:,1);
subplot(10,2,2*(NRunStart-1)+1)
semshade(speed,0.3,'g',[1:1:size(speed,2)]);
M=nanmean(speed,1);
E=nansem(speed,1);
y2=max(M)+max(E);
y1=min(M)-min(E);
hold on
line([beforeAfterZone+1 beforeAfterZone+1],[y1 y2])
hold on
line([beforeAfterZone+15 beforeAfterZone+15],[y1 y2])
hold on
line([beforeAfterZone+15+exp beforeAfterZone+15+exp],[y1 y2])
hold on
line([beforeAfterZone-exp beforeAfterZone-exp],[y1 y2])

% if N==1 | N==2;
% ylim([25 44])
% else
%     ylim([28 44]);
% end



subplot(10,2,NRunStart*2);
bar([1:1:5],nanmean(meanSpeed,1),'FaceColor','w');
hold on
errorbar([1:1:5],nanmean(meanSpeed,1),nansem(meanSpeed,1),'k')
% ylim([0 45])
[~,p1]=ttest(meanSpeed(:,1),meanSpeed(:,2));
[~,p2]=ttest(meanSpeed(:,1),meanSpeed(:,3));
[~,p3]=ttest(meanSpeed(:,1),meanSpeed(:,4));
[~,p4]=ttest(meanSpeed(:,1),meanSpeed(:,5));
title(num2str([p1 p2 p3 p4]))
end
end
before=meanSpeedAll{N};
tightfig;
saveas(gcf,'aligniningStimZoneLowerLickFarApart50cmExp_combineStimDelay5PreEarly_stim2.fig')
print -painters -depsc aligniningStimZoneLowerLickFarApart50cmExp_combineStimDelay5PreEarly_stim2.eps

%% predictive slowing

beforeSpeedDiff=[];
for n=1:length(before);
    % beforeSpeedDiff(:,n)=(before{n}(:,2)-before{n}(:,1))./(before{n}(:,1)+before{n}(:,2));
    beforeSpeedDiff(:,n)=(before{n}(:,2)-before{n}(:,1));
end

durSpeedDiff=[];
for n=1:length(dur);
    % durSpeedDiff(:,n)=(dur{n}(:,2)-dur{n}(:,1))./(dur{n}(:,1)+dur{n}(:,2));
    durSpeedDiff(:,n)=(dur{n}(:,2)-dur{n}(:,1));
end

afterSpeedDiff=[];
for n=1:length(after);
    % afterSpeedDiff(:,n)=(after{n}(:,2)-after{n}(:,1))./(after{n}(:,1)+after{n}(:,2));
    afterSpeedDiff(:,n)=(after{n}(:,2)-after{n}(:,1));
end

all=[beforeSpeedDiff durSpeedDiff afterSpeedDiff];


figure
subplot(311)
errorbar([1:1:size(all,2)],nanmean(all,1),nansem(all,1))

pAll=[];
for n=1:size(all,2)
    [~,pAll(n)]=ttest(all(:,n),0);
    if pAll(n)<=0.05;
        hold on
        text(n,11,'*','Color','r');
    end
end


subplot(312)

precedingbeforeSpeed=[];
for n=1:length(before);
    % beforeSpeedDiff(:,n)=(before{n}(:,2)-before{n}(:,1))./(before{n}(:,1)+before{n}(:,2));
    precedingbeforeSpeed(:,n)=(before{n}(:,2));
end

precedingdurSpeed=[];
for n=1:length(dur);
    % durSpeedDiff(:,n)=(dur{n}(:,2)-dur{n}(:,1))./(dur{n}(:,1)+dur{n}(:,2));
    precedingdurSpeed(:,n)=(dur{n}(:,2));
end

precedingafterSpeed=[];
for n=1:length(after);
    % afterSpeedDiff(:,n)=(after{n}(:,2)-after{n}(:,1))./(after{n}(:,1)+after{n}(:,2));
    precedingafterSpeed(:,n)=(after{n}(:,2));
end

allPreceding=[precedingbeforeSpeed precedingdurSpeed precedingafterSpeed];


beforebeforeSpeed=[];
for n=1:length(before);
    % beforeSpeedDiff(:,n)=(before{n}(:,2)-before{n}(:,1))./(before{n}(:,1)+before{n}(:,2));
    beforebeforeSpeed(:,n)=(before{n}(:,1));
end

beforedurSpeed=[];
for n=1:length(dur);
    % durSpeedDiff(:,n)=(dur{n}(:,2)-dur{n}(:,1))./(dur{n}(:,1)+dur{n}(:,2));
    beforedurSpeed(:,n)=(dur{n}(:,1));
end

beforeafterSpeed=[];
for n=1:length(after);
    % afterSpeedDiff(:,n)=(after{n}(:,2)-after{n}(:,1))./(after{n}(:,1)+after{n}(:,2));
    beforeafterSpeed(:,n)=(after{n}(:,1));
end

allBefore=[beforebeforeSpeed beforedurSpeed beforeafterSpeed];


errorbar([1:1:size(allBefore,2)],nanmean(allBefore,1),nansem(allBefore,1),'k')
hold on
errorbar([1:1:size(allPreceding,2)],nanmean(allPreceding,1),nansem(allPreceding,1),'g')

earlyruns=[1:1:5];
lateruns=[6:1:10];

beforeEarly=[];
for n=1:length(earlyruns);
    beforeEarly=[beforeEarly;before{earlyruns(n)}(:,1:2)];
end

beforeLate=[];
for n=1:length(lateruns);
    beforeLate=[beforeLate;before{lateruns(n)}(:,1:2)];
end

durEarly=[];
for n=1:length(earlyruns);
    durEarly=[durEarly;dur{earlyruns(n)}(:,1:2)];
end

durLate=[];
for n=1:length(lateruns);
    durLate=[durLate;dur{lateruns(n)}(:,1:2)];
end

afterEarly=[];
for n=1:length(earlyruns);
    afterEarly=[afterEarly;after{earlyruns(n)}(:,1:2)];
end

afterLate=[];
for n=1:length(lateruns);
    afterLate=[afterLate;after{lateruns(n)}(:,1:2)];
end

allRuns={};
allRuns{1}=beforeEarly;
allRuns{2}=beforeLate;
allRuns{3}=durEarly;
allRuns{4}=durLate;
allRuns{5}=afterEarly;
allRuns{6}=afterLate;
% allRuns=cell2mat(allRuns);

M=[];
E=[];
p=[];
for n=1:length(allRuns);
    M(end+1:end+2)=nanmean(allRuns{n},1);
    E(end+1:end+2)=nansem(allRuns{n},1);
    [~,p(n)]=ttest(allRuns{n}(:,1),allRuns{n}(:,2));
end

subplot(313)
bar([1 2 4 5 7 8 10 11 13 14 16 17],M)
hold on
errorbar([1 2 4 5 7 8 10 11 13 14 16 17],M,E,'k.')
title([num2str(p)])

saveas(gcf,'predictSlowMemory_stim2.fig')
print -painters -depsc predictSlowMemory_stim2.eps









