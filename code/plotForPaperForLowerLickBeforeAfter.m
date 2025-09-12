%this code only calculates the slowing but separately calculated before and
%after to see wehther the mice learned

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


for N=1:length(allIdxAll);
mouseIdx=allIdxAll{N};

position=positionSegAllData(mouseIdx);
stim=stimZoneIdxAllDataFix(mouseIdx);
positionAll={};


for n=1:length(position); %per session
    for m=1:length(position{n}) %per run
        positionAll{n,m}=position{n}{m}(stim{n}{m});
    end
end

%ONLY TAKE THE POSITION THAT ARE FAR Enough

for n=1:size(positionAll,1);
    for m=1:size(positionAll,2);
    positionKeep=[];
    loc=nanmean(positionAll{n,m},2);
    loc=[0;loc;400];
    dis=diff(loc);
    for stimN=1:3;

  % %calculated the oringial distances of CS, when N=1;
% disAll=[];
% for n=1:length(positionAll);
%     positionKeep=[];
%     loc=nanmean(positionAll{n,1},2);
%     loc=[0;loc;400];
%     disAll(n,:)=diff(loc);
% end
% 
% min(disAll,[],1)
% 72.2173497465605	99.3820130340334	115.315843591600	109.841665459812   

        if stimN==1;%the first stim: can be close to track start but need to be away from the next one
            if dis(stimN+1)>=positionDisThresh&dis(stimN)>=72;
                positionKeep(end+1,:)=positionAll{n,m}(stimN,:);
            end
        elseif stimN==2;%the second stim: needs to be away from the one before and after
            if dis(stimN)>=positionDisThresh&dis(stimN+1)>=positionDisThresh;
               positionKeep(end+1,:)=positionAll{n,m}(stimN,:);
            end
        else %the 3rd stim: needs to be away from the one before
            if dis(stimN)>=positionDisThresh&dis(stimN+1)>=positionDisThresh;
                positionKeep(end+1,:)=positionAll{n,m}(stimN,:);
            end

        end
    end
    positionAll{n,m}=positionKeep;
    end
end

startBinStim={};%each row is a session/run: %run 1-10 for session1, then run 1-10 for session2... the three numbers are the start point of the three stim
for session=1:size(positionAll,1)
    for run=1:size(positionAll,2);
        if ~isempty(positionAll{session,run})

        startPoint=positionAll{session,run}(:,1);
        startBinStim{end+1}=ceil(startPoint);
        else
            startBinStim{end+1}=[];

    end
    end
end

%reorganize the RBR speed
RBRSpeed=[]; %run 1-10 for session1, then run 1-10 for session2...
S=allSDurRBRUseRun{N};
for session=1:length(mouseIdx); %starting from each session
    for run=1:length(S); %then each run
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
stimStart=15+beforeAfterZone;
% stimEnd=310-beforeAfterZone-15-exp+1;
stimEnd=340;

speed=[];
for n=1:length(startBinStim)
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
meanSpeedAll{N}=meanSpeed;


%normalize
% speed=speed./meanSpeed(:,1);
subplot(4,2,2*(N-1)+1)
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

if N==1 | N==2;
ylim([22 35])
else
    ylim([28 44]);
end



subplot(4,2,N*2);
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
save('meanSpeedAll5.mat','meanSpeedAll')
saveas(gcf,'aligniningStimZoneLowerLickFarApart50cmExp_combineStimDelay5.fig')
print -painters -depsc aligniningStimZoneLowerLickFarApart50cmExp_combineStimDelay5.eps

%% dowmsample of above 
% Downsample all to CHR2RS: 50 cm before after and expanded stim zone, aligning stim zone: more for RS: only take the stim region with at least 100 cm from other loc (the beginning region can have less)
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

%below we done sample all to 35 data points because after limiting the stim
%for CHR2RS with enough distances from others, there are only 35 data
%points left, so we down sample all to 35

for N=1:length(allIdxAll);
mouseIdx=allIdxAll{N};

position=positionSegAllData(mouseIdx);
stim=stimZoneIdxAllDataFix(mouseIdx);
positionAll={};


for n=1:length(position); %per session
    for m=1:length(position{n}) %per run
        positionAll{n,m}=position{n}{m}(stim{n}{m});
    end
end

%ONLY TAKE THE POSITION THAT ARE FAR Enough

for n=1:size(positionAll,1);
    for m=1:size(positionAll,2);
    positionKeep=[];
    loc=nanmean(positionAll{n,m},2);
    loc=[0;loc;400];
    dis=diff(loc);
    for stimN=1:3;

  % %calculated the oringial distances of CS, when N=1;
% disAll=[];
% for n=1:length(positionAll);
%     positionKeep=[];
%     loc=nanmean(positionAll{n,1},2);
%     loc=[0;loc;400];
%     disAll(n,:)=diff(loc);
% end
% 
% min(disAll,[],1)
% 72.2173497465605	99.3820130340334	115.315843591600	109.841665459812   

        if stimN==1;%the first stim: can be close to track start but need to be away from the next one
            if dis(stimN+1)>=positionDisThresh&dis(stimN)>=72;
                positionKeep(end+1,:)=positionAll{n,m}(stimN,:);
            end
        elseif stimN==2;%the second stim: needs to be away from the one before and after
            if dis(stimN)>=positionDisThresh&dis(stimN+1)>=positionDisThresh;
               positionKeep(end+1,:)=positionAll{n,m}(stimN,:);
            end
        else %the 3rd stim: needs to be away from the one before
            if dis(stimN)>=positionDisThresh&dis(stimN+1)>=positionDisThresh;
                positionKeep(end+1,:)=positionAll{n,m}(stimN,:);
            end

        end
    end
    positionAll{n,m}=positionKeep;
    end
end

startBinStim={};%each row is a session/run: %run 1-10 for session1, then run 1-10 for session2... the three numbers are the start point of the three stim
for session=1:size(positionAll,1)
    for run=1:size(positionAll,2);
        if ~isempty(positionAll{session,run})

        startPoint=positionAll{session,run}(:,1);
        startBinStim{end+1}=ceil(startPoint);
        else
            startBinStim{end+1}=[];

    end
    end
end

%reorganize the RBR speed
RBRSpeed=[]; %run 1-10 for session1, then run 1-10 for session2...
S=allSDurRBRUseRun{N};
for session=1:length(mouseIdx); %starting from each session
    for run=1:length(S); %then each run
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
stimStart=15+beforeAfterZone;
% stimEnd=310-beforeAfterZone-15-exp+1;
stimEnd=340;

speed=[];
for n=1:length(startBinStim)
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


i=randperm(size(speed,1));
if length(i)>=61;
i=i(1:61);
end
speed=speed(i,:);
meanSpeed=meanSpeed(i,:);
meanSpeedAll{N}=meanSpeed;
%normalize
% speed=speed./meanSpeed(:,1);
subplot(4,2,2*(N-1)+1)
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

if N==1 | N==2;
ylim([22 35])
else
    ylim([28 44]);
end



subplot(4,2,N*2);
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

saveas(gcf,'aligniningStimZoneLowerLickFarApart50cmExpDownsample_combineStimDelay5.fig')
print -painters -depsc aligniningStimZoneLowerLickFarApart50cmExpDownsample_combineStimDelay5.eps
save('oneShuffle5.mat','meanSpeedAll')

%%  whether this is because the cs data had more points
load('meanSpeedAll5.mat')

NShuffle=1000;
NTake=size(meanSpeedAll{2},1);

allDownsamplePValues={};
meanCompare={};

idx=1;%downsample CHR2
NTotal=size(meanSpeedAll{idx},1);
pValuesShuffle=[];
meanCompareShuffle=[];
for n=1:NShuffle;
    i=randperm(NTotal);
    i=i(1:NTake);
    meanSpeedUse=meanSpeedAll{idx}(i,:);
    [~,pValuesShuffle(n,1)]=ttest(meanSpeedUse(:,1),meanSpeedUse(:,2)); %before vs before exp
[~,pValuesShuffle(n,2)]=ttest(meanSpeedUse(:,1),meanSpeedUse(:,3)); %before vs dur stim
[~,pValuesShuffle(n,3)]=ttest(meanSpeedUse(:,1),meanSpeedUse(:,4)); %before vs after exp
[~,pValuesShuffle(n,4)]=ttest(meanSpeedUse(:,1),meanSpeedUse(:,5)); %before vs after

meanCompareShuffle(n,1)=nanmean(meanSpeedUse(:,2),1)-nanmean(meanSpeedUse(:,1),1);
meanCompareShuffle(n,2)=nanmean(meanSpeedUse(:,3),1)-nanmean(meanSpeedUse(:,1),1);
meanCompareShuffle(n,3)=nanmean(meanSpeedUse(:,4),1)-nanmean(meanSpeedUse(:,1),1);
meanCompareShuffle(n,4)=nanmean(meanSpeedUse(:,5),1)-nanmean(meanSpeedUse(:,1),1);

end


allDownSamplePValues{idx}=pValuesShuffle;
meanCompare{idx}=meanCompareShuffle;


RSLevelBeforeVSStim=nanmean(meanSpeedAll{2}(:,3))-nanmean(meanSpeedAll{2}(:,1));
RSLevelBeforeVSBeforeDelay=nanmean(meanSpeedAll{2}(:,2))-nanmean(meanSpeedAll{2}(:,1));
RSLevelBeforeVSAfterDelay=nanmean(meanSpeedAll{2}(:,4))-nanmean(meanSpeedAll{2}(:,1));
   
save('allDownSamplePValuesCombineStimDelay5.mat','allDownSamplePValues','meanCompare')

%
%plot p value on a log scale

NShuffle=1000;

figure
 subplot(131)
p=allDownSamplePValues{1}(:,1);
m=meanCompare{1}(:,1);
r=RSLevelBeforeVSBeforeDelay;
semilogy(m,p, 'k.');  % semilogy uses log scale on the y-axis
hold on
semilogy([min(m) max(m)],[0.05 0.05], 'r-');  % semilogy uses log scale on the y-axis
hold on
semilogy([0 0],[min(p) max(p)], 'r-');  % semilogy uses log scale on the y-axis

hold on
semilogy([r r],[min(p) max(p)], 'g-');  % semilogy uses log scale on the y-axis

p1=length(find(p<=0.05))/1000;
p2=length(find(m<0)/NShuffle);
p3=length(find(m<r)/NShuffle);
title(['before',num2str([p1 p2 p3])]);
axis square

 subplot(132)
p=allDownSamplePValues{1}(:,2);
m=meanCompare{1}(:,2);
r=RSLevelBeforeVSStim;
semilogy(m,p, 'k.');  % semilogy uses log scale on the y-axis
hold on
semilogy([min(m) max(m)],[0.05 0.05], 'r-');  % semilogy uses log scale on the y-axis
hold on
semilogy([0 0],[min(p) max(p)], 'r-');  % semilogy uses log scale on the y-axis

hold on
semilogy([r r],[min(p) max(p)], 'g-');  % semilogy uses log scale on the y-axis

p1=length(find(p<=0.05))/1000;
p2=length(find(m<0)/NShuffle);
p3=length(find(m<r)/NShuffle);
title(['stim',num2str([p1 p2 p3])]);
axis square

 subplot(133)
p=allDownSamplePValues{1}(:,3);
m=meanCompare{1}(:,3);
r=RSLevelBeforeVSAfterDelay;
semilogy(m,p, 'k.');  % semilogy uses log scale on the y-axis
hold on
semilogy([min(m) max(m)],[0.05 0.05], 'r-');  % semilogy uses log scale on the y-axis
hold on
semilogy([0 0],[min(p) max(p)], 'r-');  % semilogy uses log scale on the y-axis

hold on
semilogy([r r],[min(p) max(p)], 'g-');  % semilogy uses log scale on the y-axis

p1=length(find(p<=0.05))/1000;
p2=length(find(m<0)/NShuffle);
p3=length(find(m<r)/NShuffle);
title(['after',num2str([p1 p2 p3])]);
axis square

 saveas(gcf,'pValuesForDownsampleCHRRCSCombineStimDelay_both5.fig')
print -painters -depsc pValuesForDownsampleCHRRCSCombineStimDelay_both5.eps