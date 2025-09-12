
%this stretagy did not work: cannot find a common baseline.

%Will just use 800 here
idxCHR2RS=[1:15];
idxCHR2CS=[27:44];
idxGFPRS=[16:26];
idxGFPCS=[45:56];


%%

%use 10 if RBR is using A/mean(A) (normalize to 1)
beforeDistance=20;
% beforeDistance=15;
afterDistance=34;
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\lickSegAllData.mat');
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\rewardSegAllData.mat');
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\positionSegAllData.mat');
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\threshold.mat');
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\lickSegAllDataBefore.mat');
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\rewardSegAllDataBefore.mat');
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\positionSegAllDataBefore.mat');
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\threshold.mat');
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\lickSegAllDataAfter.mat');
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\rewardSegAllDataAfter.mat');
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\positionSegAllDataAfter.mat');
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\threshold.mat');


edges=linspace(0,366,81);
width=nanmean(diff(edges));
binWindow=5;
cues=[63 85;163 185;279 301;374 396];
stimLoc=[];
stimLoc(1,:)=[366-15 366];
binWidth=1;
rewardLoc=366/binWidth;


load('Z:\papers\2025RewardOptogenetics\4mData\250128_data\cueTemp_4m.mat')
countBeforeAllCHR2CS=[];
countStimAllCHR2CS=[];
idx=idxCHR2CS;

for n=1:length(idx);
    [lickUsePositionStimAll,~,~,~] = absPreRewardLicking_RBR(lickSegAllData(idx(n)),rewardSegAllData(idx(n)),positionSegAllData(idx(n)),threshold(idx(n)),beforeDistance,afterDistance);
[lickUsePositionBeforeAll,~,~,~] = absPreRewardLicking_RBR(lickSegAllDataBefore(idx(n)),rewardSegAllDataBefore(idx(n)),positionSegAllDataBefore(idx(n)),threshold(idx(n)),beforeDistance,afterDistance);
lickUsePositionStimAll=lickUsePositionStimAll(lickUsePositionStimAll<366);
lickUsePositionBeforeAll=lickUsePositionBeforeAll(lickUsePositionBeforeAll<366);
h=histogram(lickUsePositionBeforeAll,edges);
countBeforeAllCHR2CS(n,:)=h.Values;

h=histogram(lickUsePositionStimAll,edges);
countStimAllCHR2CS(n,:)=h.Values;


end

countBeforeAllCHR2RS=[];
countStimAllCHR2RS=[];
idx=idxCHR2RS;

for n=1:length(idx);
    [lickUsePositionStimAll,~,~,~] = absPreRewardLicking_RBR(lickSegAllData(idx(n)),rewardSegAllData(idx(n)),positionSegAllData(idx(n)),threshold(idx(n)),beforeDistance,afterDistance);
[lickUsePositionBeforeAll,~,~,~] = absPreRewardLicking_RBR(lickSegAllDataBefore(idx(n)),rewardSegAllDataBefore(idx(n)),positionSegAllDataBefore(idx(n)),threshold(idx(n)),beforeDistance,afterDistance);
lickUsePositionStimAll=lickUsePositionStimAll(lickUsePositionStimAll<366);
lickUsePositionBeforeAll=lickUsePositionBeforeAll(lickUsePositionBeforeAll<366);
h=histogram(lickUsePositionBeforeAll,edges);
countBeforeAllCHR2RS(n,:)=h.Values;

h=histogram(lickUsePositionStimAll,edges);
countStimAllCHR2RS(n,:)=h.Values;


end

countBeforeAllGFPCS=[];
countStimAllGFPCS=[];
idx=idxGFPCS;

for n=1:length(idx);
    [lickUsePositionStimAll,~,~,~] = absPreRewardLicking_RBR(lickSegAllData(idx(n)),rewardSegAllData(idx(n)),positionSegAllData(idx(n)),threshold(idx(n)),beforeDistance,afterDistance);
[lickUsePositionBeforeAll,~,~,~] = absPreRewardLicking_RBR(lickSegAllDataBefore(idx(n)),rewardSegAllDataBefore(idx(n)),positionSegAllDataBefore(idx(n)),threshold(idx(n)),beforeDistance,afterDistance);
lickUsePositionStimAll=lickUsePositionStimAll(lickUsePositionStimAll<366);
lickUsePositionBeforeAll=lickUsePositionBeforeAll(lickUsePositionBeforeAll<366);
h=histogram(lickUsePositionBeforeAll,edges);
countBeforeAllGFPCS(n,:)=h.Values;

h=histogram(lickUsePositionStimAll,edges);
countStimAllGFPCS(n,:)=h.Values;



end

countBeforeAllGFPRS=[];
countStimAllGFPRS=[];
idx=idxGFPRS;

for n=1:length(idx);
    [lickUsePositionStimAll,~,~,~] = absPreRewardLicking_RBR(lickSegAllData(idx(n)),rewardSegAllData(idx(n)),positionSegAllData(idx(n)),threshold(idx(n)),beforeDistance,afterDistance);
[lickUsePositionBeforeAll,~,~,~] = absPreRewardLicking_RBR(lickSegAllDataBefore(idx(n)),rewardSegAllDataBefore(idx(n)),positionSegAllDataBefore(idx(n)),threshold(idx(n)),beforeDistance,afterDistance);
lickUsePositionStimAll=lickUsePositionStimAll(lickUsePositionStimAll<366);
lickUsePositionBeforeAll=lickUsePositionBeforeAll(lickUsePositionBeforeAll<366);
h=histogram(lickUsePositionBeforeAll,edges);
countBeforeAllGFPRS(n,:)=h.Values;

h=histogram(lickUsePositionStimAll,edges);
countStimAllGFPRS(n,:)=h.Values;

end

%%

a1=[countBeforeAllCHR2CS countStimAllCHR2CS];
a2=[countBeforeAllCHR2RS countStimAllCHR2RS];
a3=[countBeforeAllGFPCS countStimAllGFPCS];
a4=[countBeforeAllGFPRS countStimAllGFPRS];
% a1=reshape(a1,[size(a1,1)*size(a1,2) 1]);
% a2=reshape(a2,[size(a2,1)*size(a2,2) 1]);
% a3=reshape(a3,[size(a3,1)*size(a3,2) 1]);

a=[a1;a2;a3;a4];
a=reshape(a,[size(a,1)*size(a,2) 1]);
p=prctile(a,99.9);

idxRemove={};
idxRemove{1}=[];
for n=1:size(a1,1);
    if ~isempty(find(a1(n,:)>=p));
        idxRemove{1}(end+1)=n;
    end
end

idxRemove{2}=[];
for n=1:size(a2,1);
    if ~isempty(find(a2(n,:)>=p));
        idxRemove{2}(end+1)=n;
    end
end

idxRemove{3}=[];
for n=1:size(a3,1);
    if ~isempty(find(a3(n,:)>=p));
        idxRemove{3}(end+1)=n;
    end
end

idxRemove{4}=[];
for n=1:size(a4,1);
    if ~isempty(find(a4(n,:)>=p));
        idxRemove{4}(end+1)=n;
    end
end