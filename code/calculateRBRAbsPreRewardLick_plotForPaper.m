load('idxCHR2CS');
load('idxCHR2RS');
load('idxGFPCS');
load('idxGFPRS');

% %this is based on these from script: used all mice, not just good ones or the ones with good licking. This the good script to use: calculateRBRAbsPreRewardLick.m
%%%%%%%%%%%%%%%
%note: here: good performers are used for CHR2RS because they together
%showed better licking that are comparable with those for good performers.

% also removed the cells, which in some sessions (before, stim, or after)
% has very high licks and it is > 500 licks per bin

% idxCHR2RS=[1:15];
% idxCHR2CS=[27:44];
% idxGFPRS=[16:26];
% idxGFPCS=[45:56];
% allIdxAll={};%REMOVE CELLS WITH EXTREAMLY HIGH LICKS
% allIdxAll{1}=setdiff([1:1:length(idxCHR2RS)],[8 9 13]);%CHR2RS
% allIdxAll{2}=setdiff([1:1:length(idxCHR2CS)],[12 14]); %chr2 CS
% allIdxAll{3}=setdiff([1:1:length(idxGFPRS)],[]);%GFPRS
% allIdxAll{4}=setdiff([1:1:length(idxGFPCS)],[]); %GFP CS
% idxCHR2CS=idxCHR2CS(allIdxAll{2});
% idxCHR2RS=idxCHR2RS(allIdxAll{1});
% idxGFPCS=idxGFPCS(allIdxAll{4});
% idxGFPRS=idxGFPRS(allIdxAll{3});
load('threshold.mat')
beforeDistance=15;
afterDistance=30;
edges=linspace(0,366,81);
width=nanmean(diff(edges));
xvalues=edges(1:end-1)+width/2;
cues=[63 85;163 185;279 301;374 396];
stimLoc=[];
stimLoc(1,:)=[63+3.5 63+3.5+15];
stimLoc(2,:)=[163+3.5 163+3.5+15];
stimLoc(3,:)=[279+3.5 279+3.5+15];

binWidth=1;
rewardLoc=366/binWidth;
% stimExp=7.5; %expanding the stimulation zone to both sides
% stimExp=5; %expanding the stimulation zone to both sides
stimExp=6; %expanding the stimulation zone to both sides
samplingSpacing=1;%how fine to sample out stim area
%% chr2CS
idx=idxCHR2CS;

%%%during stim
load('lickSegAllData.mat')
load('rewardSegAllData.mat')
load('positionSegAllData.mat')
[lickUsePositionStimRBR,~,~,~] = absPreRewardLicking_RBRRBRLick(lickSegAllData(idx),rewardSegAllData(idx),positionSegAllData(idx),threshold(idx),beforeDistance,afterDistance);
for session=1:length(lickUsePositionStimRBR);
    for lap=1:length(lickUsePositionStimRBR{session});
lickUsePositionStimRBR{session}{lap}=lickUsePositionStimRBR{session}{lap}(lickUsePositionStimRBR{session}{lap}<366);
    end
end
countStimRBR=[];
allStimRBR={};
for lap=1:10;
    allStimRBR{lap}=[];
    for session=1:length(lickUsePositionStimRBR);
        allStimRBR{lap}=[allStimRBR{lap};lickUsePositionStimRBR{session}{lap}];
    end
    h=histogram(allStimRBR{lap},edges);
    countStimRBR(lap,:)=h.Values;
end

%%%before
load('lickSegAllDataBefore.mat')
load('rewardSegAllDataBefore.mat')
load('positionSegAllDataBefore.mat')
[lickUsePositionBeforeRBR,~,~,~] = absPreRewardLicking_RBRRBRLick(lickSegAllDataBefore(idx),rewardSegAllDataBefore(idx),positionSegAllDataBefore(idx),threshold(idx),beforeDistance,afterDistance);
for session=1:length(lickUsePositionBeforeRBR);
    for lap=1:length(lickUsePositionBeforeRBR{session});
lickUsePositionBeforeRBR{session}{lap}=lickUsePositionBeforeRBR{session}{lap}(lickUsePositionBeforeRBR{session}{lap}<366);
    end
end

% group licks across runs
countRBRBefore=[];
allBeforeRBR={};
for lap=1:10;
    allBeforeRBR{lap}=[];
    for session=1:length(lickUsePositionBeforeRBR);
        allBeforeRBR{lap}=[allBeforeRBR{lap};lickUsePositionBeforeRBR{session}{lap}];
    end
    h=histogram(allBeforeRBR{lap},edges);
    countRBRBefore(lap,:)=h.Values;
end

%%%after
load('lickSegAllDataAfter.mat')
load('rewardSegAllDataAfter.mat')
load('positionSegAllDataAfter.mat')
[lickUsePositionAfterRBR,~,~,~] = absPreRewardLicking_RBRRBRLick(lickSegAllDataAfter(idx),rewardSegAllDataAfter(idx),positionSegAllDataAfter(idx),threshold(idx),beforeDistance,afterDistance);
for session=1:length(lickUsePositionAfterRBR);
    for lap=1:length(lickUsePositionAfterRBR{session});
lickUsePositionAfterRBR{session}{lap}=lickUsePositionAfterRBR{session}{lap}(lickUsePositionAfterRBR{session}{lap}<366);
    end
end

%caculate percentile
[lickPercentileBefore,lickPercentileStim,lickPercentileAfter] = lickPercCueStim(lickUsePositionBeforeRBR,lickUsePositionStimRBR,lickUsePositionAfterRBR,stimLoc,stimExp,samplingSpacing);
save('LickPercentNearStimCHR2CS.mat','lickPercentileBefore','lickPercentileStim','lickPercentileAfter')
save('lickPositionCHR2CS.mat','lickUsePositionBeforeRBR','lickUsePositionStimRBR','lickUsePositionAfterRBR')

% group licks across runs
countRBRAfter=[];
allAfterRBR={};
for lap=1:10;
    allAfterRBR{lap}=[];
    for session=1:length(lickUsePositionAfterRBR);
        allAfterRBR{lap}=[allAfterRBR{lap};lickUsePositionAfterRBR{session}{lap}];
    end
    h=histogram(allAfterRBR{lap},edges);
    countRBRAfter(lap,:)=h.Values;
end

countAll=[countRBRBefore;countStimRBR;countRBRAfter];
%
% plot
figure
for n=1:size(countAll,1);
    subplot(size(countAll,1),2,1+2*(n-1))
    h=max(countAll(n,:));

    for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 h h 0],'y')
    end
if ~isempty(intersect(n,[11:1:20]))
     for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 h h 0],'b')
     end
end
hold on
line([rewardLoc rewardLoc],[0 h],'Color','r')
    hold on
    bar(xvalues,countAll(n,:),'FaceColor','m','EdgeColor','none');
        set(gca,'xticklabel',[])
    % set(gca,'yticklabel',[])
    xlim([0 400])
end




% every 5
% nDay=[0 3 3 4 3 3 4 3 3 4];
nDay=[0 1 3 3 3 1 3 3 3 1 3 3 3];
countAllAvg=[];
for n=1:size(nDay,2)-1;
    idx=[sum(nDay(1:n),2)+1:1:sum(nDay(1:n+1),2)];
    countAllAvg(n,:)=nanmean(countAll(idx,:),1);
end

allLicksPosition=[allBeforeRBR allStimRBR allAfterRBR];

lickPositionGroup={};
for n=1:size(nDay,2)-1;
    idx=[sum(nDay(1:n),2)+1:1:sum(nDay(1:n+1),2)];
    lickPositionGroup{n}=cell2mat(allLicksPosition(idx)');
end



% plot

for n=1:size(countAllAvg,1);
    % subplot(size(countAll,1),2,[(3*(n-1)+1)*2:2:(3*(n-1)+3)*2])
    subplot(size(countAll,1),2,[(2*(n-1)+1)*2:2:(2*(n-1)+2)*2])
    h=max(countAllAvg(n,:));

    for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 h h 0],'y')
end
if ~isempty(intersect(n,[5:1:8]))
     for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 h h 0],'b')
     end
end
hold on
line([rewardLoc rewardLoc],[0 h],'Color','r')
    hold on
    i=find(countAllAvg(n,:)>0);
    bar(xvalues,countAllAvg(n,:),'FaceColor','m','EdgeColor','none');
    hold on
    [p,x]=ksdensity(lickPositionGroup{n},'width',12);
    plot(x,p/max(p)*h)
       set(gca,'xticklabel',[])
    % set(gca,'yticklabel',[])
    xlim([0 400])
end
tightfig
    
saveas(gcf,'RBRLickDistri_CHR2CS.fig')
print -painters -depsc RBRLickDistri_CHR2CS.eps

% %% correlation between lick position and stim position
% stimPosition=zeros(1,366);
% for n=1:size(stimLoc,1);
%     stimPosition(stimLoc(n,1):stimLoc(n,2))=1;
% end
% 
% corrStimPosition=[];
% for n=1:length(lickPositionGroup2)
%     lick=lickPositionGroup2{n};
%     [p,x]=ksdensity(lick,[1:1:366]);
%     corrStimPosition(n)=corr(stimPosition',p');
% end
% figure,plot(corrStimPosition)

%% chr2RS
idx=idxCHR2RS;

%%%during stim
load('lickSegAllData.mat')
load('rewardSegAllData.mat')
load('positionSegAllData.mat')
[lickUsePositionStimRBR,~,~,~] = absPreRewardLicking_RBRRBRLick(lickSegAllData(idx),rewardSegAllData(idx),positionSegAllData(idx),threshold(idx),beforeDistance,afterDistance);
for session=1:length(lickUsePositionStimRBR);
    for lap=1:length(lickUsePositionStimRBR{session});
lickUsePositionStimRBR{session}{lap}=lickUsePositionStimRBR{session}{lap}(lickUsePositionStimRBR{session}{lap}<366);
    end
end
countStimRBR=[];
allStimRBR={};
for lap=1:10;
    allStimRBR{lap}=[];
    for session=1:length(lickUsePositionStimRBR);
        allStimRBR{lap}=[allStimRBR{lap};lickUsePositionStimRBR{session}{lap}];
    end
    h=histogram(allStimRBR{lap},edges);
    countStimRBR(lap,:)=h.Values;
end

%%%before
load('lickSegAllDataBefore.mat')
load('rewardSegAllDataBefore.mat')
load('positionSegAllDataBefore.mat')
[lickUsePositionBeforeRBR,~,~,~] = absPreRewardLicking_RBRRBRLick(lickSegAllDataBefore(idx),rewardSegAllDataBefore(idx),positionSegAllDataBefore(idx),threshold(idx),beforeDistance,afterDistance);
for session=1:length(lickUsePositionBeforeRBR);
    for lap=1:length(lickUsePositionBeforeRBR{session});
lickUsePositionBeforeRBR{session}{lap}=lickUsePositionBeforeRBR{session}{lap}(lickUsePositionBeforeRBR{session}{lap}<366);
    end
end

% group licks across runs
countRBRBefore=[];
allBeforeRBR={};
for lap=1:10;
    allBeforeRBR{lap}=[];
    for session=1:length(lickUsePositionBeforeRBR);
        allBeforeRBR{lap}=[allBeforeRBR{lap};lickUsePositionBeforeRBR{session}{lap}];
    end
    h=histogram(allBeforeRBR{lap},edges);
    countRBRBefore(lap,:)=h.Values;
end

%%%after
load('lickSegAllDataAfter.mat')
load('rewardSegAllDataAfter.mat')
load('positionSegAllDataAfter.mat')
[lickUsePositionAfterRBR,~,~,~] = absPreRewardLicking_RBRRBRLick(lickSegAllDataAfter(idx),rewardSegAllDataAfter(idx),positionSegAllDataAfter(idx),threshold(idx),beforeDistance,afterDistance);
for session=1:length(lickUsePositionAfterRBR);
    for lap=1:length(lickUsePositionAfterRBR{session});
lickUsePositionAfterRBR{session}{lap}=lickUsePositionAfterRBR{session}{lap}(lickUsePositionAfterRBR{session}{lap}<366);
    end
end

%caculate percentile
[lickPercentileBefore,lickPercentileStim,lickPercentileAfter] = lickPercCueStim(lickUsePositionBeforeRBR,lickUsePositionStimRBR,lickUsePositionAfterRBR,stimLoc,stimExp,samplingSpacing);
save('LickPercentNearStimCHR2RS.mat','lickPercentileBefore','lickPercentileStim','lickPercentileAfter')
save('lickPositionCHR2RS.mat','lickUsePositionBeforeRBR','lickUsePositionStimRBR','lickUsePositionAfterRBR')

% group licks across runs
countRBRAfter=[];
allAfterRBR={};
for lap=1:10;
    allAfterRBR{lap}=[];
    for session=1:length(lickUsePositionAfterRBR);
        allAfterRBR{lap}=[allAfterRBR{lap};lickUsePositionAfterRBR{session}{lap}];
    end
    h=histogram(allAfterRBR{lap},edges);
    countRBRAfter(lap,:)=h.Values;
end

countAll=[countRBRBefore;countStimRBR;countRBRAfter];
%
% plot
figure
for n=1:size(countAll,1);
    subplot(size(countAll,1),2,1+2*(n-1))
    h=max(countAll(n,:));

    for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 h h 0],'y')
    end
if ~isempty(intersect(n,[11:1:20]))
     for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 h h 0],'b')
     end
end
hold on
line([rewardLoc rewardLoc],[0 h],'Color','r')
    hold on
    bar(xvalues,countAll(n,:),'FaceColor',[0.5 0.5 0.5],'EdgeColor','none');
     set(gca,'xticklabel',[])
    % set(gca,'yticklabel',[])
    xlim([0 400])
end




% every 5
% nDay=[0 3 3 4 3 3 4 3 3 4];
nDay=[0 1 3 3 3 1 3 3 3 1 3 3 3];
countAllAvg=[];
for n=1:size(nDay,2)-1;
    idx=[sum(nDay(1:n),2)+1:1:sum(nDay(1:n+1),2)];
    countAllAvg(n,:)=nanmean(countAll(idx,:),1);
end

allLicksPosition=[allBeforeRBR allStimRBR allAfterRBR];

lickPositionGroup={};
for n=1:size(nDay,2)-1;
    idx=[sum(nDay(1:n),2)+1:1:sum(nDay(1:n+1),2)];
    lickPositionGroup{n}=cell2mat(allLicksPosition(idx)');
end



% plot

for n=1:size(countAllAvg,1);
    % subplot(size(countAll,1),2,[(3*(n-1)+1)*2:2:(3*(n-1)+3)*2])
    subplot(size(countAll,1),2,[(2*(n-1)+1)*2:2:(2*(n-1)+2)*2])
    h=max(countAllAvg(n,:));

    for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 h h 0],'y')
end
if ~isempty(intersect(n,[5:1:8]))
     for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 h h 0],'b')
     end
end
hold on
line([rewardLoc rewardLoc],[0 h],'Color',[0.5 0.5 0.5])
    hold on
    i=find(countAllAvg(n,:)>0);
    bar(xvalues,countAllAvg(n,:),'FaceColor',[0.5 0.5 0.5],'EdgeColor','none');
    hold on
    [p,x]=ksdensity(lickPositionGroup{n},'width',12);
    plot(x,p/max(p)*h)
        set(gca,'xticklabel',[])
    % set(gca,'yticklabel',[])
   
    xlim([0 400])
end
tightfig
    
saveas(gcf,'RBRLickDistri_CHR2RS.fig')
   print -painters -depsc RBRLickDistri_CHR2RS.eps 

%% gfpCS
idx=idxGFPCS;

%%%during stim
load('lickSegAllData.mat')
load('rewardSegAllData.mat')
load('positionSegAllData.mat')
[lickUsePositionStimRBR,~,~,~] = absPreRewardLicking_RBRRBRLick(lickSegAllData(idx),rewardSegAllData(idx),positionSegAllData(idx),threshold(idx),beforeDistance,afterDistance);
for session=1:length(lickUsePositionStimRBR);
    for lap=1:length(lickUsePositionStimRBR{session});
lickUsePositionStimRBR{session}{lap}=lickUsePositionStimRBR{session}{lap}(lickUsePositionStimRBR{session}{lap}<366);
    end
end
countStimRBR=[];
allStimRBR={};
for lap=1:10;
    allStimRBR{lap}=[];
    for session=1:length(lickUsePositionStimRBR);
        allStimRBR{lap}=[allStimRBR{lap};lickUsePositionStimRBR{session}{lap}];
    end
    h=histogram(allStimRBR{lap},edges);
    countStimRBR(lap,:)=h.Values;
end

%%%before
load('lickSegAllDataBefore.mat')
load('rewardSegAllDataBefore.mat')
load('positionSegAllDataBefore.mat')
[lickUsePositionBeforeRBR,~,~,~] = absPreRewardLicking_RBRRBRLick(lickSegAllDataBefore(idx),rewardSegAllDataBefore(idx),positionSegAllDataBefore(idx),threshold(idx),beforeDistance,afterDistance);
for session=1:length(lickUsePositionBeforeRBR);
    for lap=1:length(lickUsePositionBeforeRBR{session});
lickUsePositionBeforeRBR{session}{lap}=lickUsePositionBeforeRBR{session}{lap}(lickUsePositionBeforeRBR{session}{lap}<366);
    end
end

% group licks across runs
countRBRBefore=[];
allBeforeRBR={};
for lap=1:10;
    allBeforeRBR{lap}=[];
    for session=1:length(lickUsePositionBeforeRBR);
        allBeforeRBR{lap}=[allBeforeRBR{lap};lickUsePositionBeforeRBR{session}{lap}];
    end
    h=histogram(allBeforeRBR{lap},edges);
    countRBRBefore(lap,:)=h.Values;
end

%%%after
load('lickSegAllDataAfter.mat')
load('rewardSegAllDataAfter.mat')
load('positionSegAllDataAfter.mat')
[lickUsePositionAfterRBR,~,~,~] = absPreRewardLicking_RBRRBRLick(lickSegAllDataAfter(idx),rewardSegAllDataAfter(idx),positionSegAllDataAfter(idx),threshold(idx),beforeDistance,afterDistance);
for session=1:length(lickUsePositionAfterRBR);
    for lap=1:length(lickUsePositionAfterRBR{session});
lickUsePositionAfterRBR{session}{lap}=lickUsePositionAfterRBR{session}{lap}(lickUsePositionAfterRBR{session}{lap}<366);
    end
end

%caculate percentile
[lickPercentileBefore,lickPercentileStim,lickPercentileAfter] = lickPercCueStim(lickUsePositionBeforeRBR,lickUsePositionStimRBR,lickUsePositionAfterRBR,stimLoc,stimExp,samplingSpacing);
save('LickPercentNearStimGFPCS.mat','lickPercentileBefore','lickPercentileStim','lickPercentileAfter')
save('lickPositionGFPCS.mat','lickUsePositionBeforeRBR','lickUsePositionStimRBR','lickUsePositionAfterRBR')


% group licks across runs
countRBRAfter=[];
allAfterRBR={};
for lap=1:10;
    allAfterRBR{lap}=[];
    for session=1:length(lickUsePositionAfterRBR);
        allAfterRBR{lap}=[allAfterRBR{lap};lickUsePositionAfterRBR{session}{lap}];
    end
    h=histogram(allAfterRBR{lap},edges);
    countRBRAfter(lap,:)=h.Values;
end

countAll=[countRBRBefore;countStimRBR;countRBRAfter];
%
% plot
figure
for n=1:size(countAll,1);
    subplot(size(countAll,1),2,1+2*(n-1))
    h=max(countAll(n,:));

    for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 h h 0],'y')
    end
if ~isempty(intersect(n,[11:1:20]))
     for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 h h 0],'b')
     end
end
hold on
line([rewardLoc rewardLoc],[0 h],'Color','r')
    hold on
    bar(xvalues,countAll(n,:),'FaceColor','m','EdgeColor','none');
        set(gca,'xticklabel',[])
    % set(gca,'yticklabel',[])
    xlim([0 400])
end




% every 5
% nDay=[0 3 3 4 3 3 4 3 3 4];
nDay=[0 1 3 3 3 1 3 3 3 1 3 3 3];
countAllAvg=[];
for n=1:size(nDay,2)-1;
    idx=[sum(nDay(1:n),2)+1:1:sum(nDay(1:n+1),2)];
    countAllAvg(n,:)=nanmean(countAll(idx,:),1);
end

allLicksPosition=[allBeforeRBR allStimRBR allAfterRBR];

lickPositionGroup={};
for n=1:size(nDay,2)-1;
    idx=[sum(nDay(1:n),2)+1:1:sum(nDay(1:n+1),2)];
    lickPositionGroup{n}=cell2mat(allLicksPosition(idx)');
end



% plot

for n=1:size(countAllAvg,1);
    % subplot(size(countAll,1),2,[(3*(n-1)+1)*2:2:(3*(n-1)+3)*2])
    subplot(size(countAll,1),2,[(2*(n-1)+1)*2:2:(2*(n-1)+2)*2])
    h=max(countAllAvg(n,:));

    for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 h h 0],'y')
end
if ~isempty(intersect(n,[5:1:8]))
     for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 h h 0],'b')
     end
end
hold on
line([rewardLoc rewardLoc],[0 h],'Color','r')
    hold on
    i=find(countAllAvg(n,:)>0);
    bar(xvalues,countAllAvg(n,:),'FaceColor','m','EdgeColor','none');
    hold on
    [p,x]=ksdensity(lickPositionGroup{n},'width',12);
    plot(x,p/max(p)*h)
      set(gca,'xticklabel',[])
    % set(gca,'yticklabel',[])
    xlim([0 400])
end
tightfig
    
saveas(gcf,'RBRLickDistri_GFPCS.fig')
print -painters -depsc RBRLickDistri_GFPCS.eps

%
% %% correlation between lick position and stim position
% stimPosition=zeros(1,366);
% for n=1:size(stimLoc,1);
%     stimPosition(stimLoc(n,1):stimLoc(n,2))=1;
% end
% 
% corrStimPosition=[];
% for n=1:length(lickPositionGroup2)
%     lick=lickPositionGroup2{n};
%     [p,x]=ksdensity(lick,[1:1:366]);
%     corrStimPosition(n)=corr(stimPosition',p');
% end
% figure,plot(corrStimPosition)

%% gfpRS
idx=idxGFPRS;

%%%during stim
load('lickSegAllData.mat')
load('rewardSegAllData.mat')
load('positionSegAllData.mat')
[lickUsePositionStimRBR,~,~,~] = absPreRewardLicking_RBRRBRLick(lickSegAllData(idx),rewardSegAllData(idx),positionSegAllData(idx),threshold(idx),beforeDistance,afterDistance);
for session=1:length(lickUsePositionStimRBR);
    for lap=1:length(lickUsePositionStimRBR{session});
lickUsePositionStimRBR{session}{lap}=lickUsePositionStimRBR{session}{lap}(lickUsePositionStimRBR{session}{lap}<366);
    end
end
countStimRBR=[];
allStimRBR={};
for lap=1:10;
    allStimRBR{lap}=[];
    for session=1:length(lickUsePositionStimRBR);
        allStimRBR{lap}=[allStimRBR{lap};lickUsePositionStimRBR{session}{lap}];
    end
    h=histogram(allStimRBR{lap},edges);
    countStimRBR(lap,:)=h.Values;
end

%%%before
load('lickSegAllDataBefore.mat')
load('rewardSegAllDataBefore.mat')
load('positionSegAllDataBefore.mat')
[lickUsePositionBeforeRBR,~,~,~] = absPreRewardLicking_RBRRBRLick(lickSegAllDataBefore(idx),rewardSegAllDataBefore(idx),positionSegAllDataBefore(idx),threshold(idx),beforeDistance,afterDistance);
for session=1:length(lickUsePositionBeforeRBR);
    for lap=1:length(lickUsePositionBeforeRBR{session});
lickUsePositionBeforeRBR{session}{lap}=lickUsePositionBeforeRBR{session}{lap}(lickUsePositionBeforeRBR{session}{lap}<366);
    end
end

% group licks across runs
countRBRBefore=[];
allBeforeRBR={};
for lap=1:10;
    allBeforeRBR{lap}=[];
    for session=1:length(lickUsePositionBeforeRBR);
        allBeforeRBR{lap}=[allBeforeRBR{lap};lickUsePositionBeforeRBR{session}{lap}];
    end
    h=histogram(allBeforeRBR{lap},edges);
    countRBRBefore(lap,:)=h.Values;
end

%%%after
load('lickSegAllDataAfter.mat')
load('rewardSegAllDataAfter.mat')
load('positionSegAllDataAfter.mat')
[lickUsePositionAfterRBR,~,~,~] = absPreRewardLicking_RBRRBRLick(lickSegAllDataAfter(idx),rewardSegAllDataAfter(idx),positionSegAllDataAfter(idx),threshold(idx),beforeDistance,afterDistance);
for session=1:length(lickUsePositionAfterRBR);
    for lap=1:length(lickUsePositionAfterRBR{session});
lickUsePositionAfterRBR{session}{lap}=lickUsePositionAfterRBR{session}{lap}(lickUsePositionAfterRBR{session}{lap}<366);
    end
end

%caculate percentile
[lickPercentileBefore,lickPercentileStim,lickPercentileAfter] = lickPercCueStim(lickUsePositionBeforeRBR,lickUsePositionStimRBR,lickUsePositionAfterRBR,stimLoc,stimExp,samplingSpacing);
save('LickPercentNearStimGFPRS.mat','lickPercentileBefore','lickPercentileStim','lickPercentileAfter')
save('lickPositionGFPRS.mat','lickUsePositionBeforeRBR','lickUsePositionStimRBR','lickUsePositionAfterRBR')


% group licks across runs
countRBRAfter=[];
allAfterRBR={};
for lap=1:10;
    allAfterRBR{lap}=[];
    for session=1:length(lickUsePositionAfterRBR);
        allAfterRBR{lap}=[allAfterRBR{lap};lickUsePositionAfterRBR{session}{lap}];
    end
    h=histogram(allAfterRBR{lap},edges);
    countRBRAfter(lap,:)=h.Values;
end

countAll=[countRBRBefore;countStimRBR;countRBRAfter];
%
% plot
figure
for n=1:size(countAll,1);
    subplot(size(countAll,1),2,1+2*(n-1))
    h=max(countAll(n,:));

    for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 h h 0],'y')
    end
if ~isempty(intersect(n,[11:1:20]))
     for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 h h 0],'b')
     end
end
hold on
line([rewardLoc rewardLoc],[0 h],'Color','r')
    hold on
    bar(xvalues,countAll(n,:),'FaceColor',[0.5 0.5 0.5],'EdgeColor','none');
        set(gca,'xticklabel',[])
    % set(gca,'yticklabel',[])
    xlim([0 400])
end




% every 5
% nDay=[0 3 3 4 3 3 4 3 3 4];
nDay=[0 1 3 3 3 1 3 3 3 1 3 3 3];
countAllAvg=[];
for n=1:size(nDay,2)-1;
    idx=[sum(nDay(1:n),2)+1:1:sum(nDay(1:n+1),2)];
    countAllAvg(n,:)=nanmean(countAll(idx,:),1);
end

allLicksPosition=[allBeforeRBR allStimRBR allAfterRBR];

lickPositionGroup={};
for n=1:size(nDay,2)-1;
    idx=[sum(nDay(1:n),2)+1:1:sum(nDay(1:n+1),2)];
    lickPositionGroup{n}=cell2mat(allLicksPosition(idx)');
end



% plot

for n=1:size(countAllAvg,1);
    % subplot(size(countAll,1),2,[(3*(n-1)+1)*2:2:(3*(n-1)+3)*2])
    subplot(size(countAll,1),2,[(2*(n-1)+1)*2:2:(2*(n-1)+2)*2])
    h=max(countAllAvg(n,:));

    for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 h h 0],'y')
end
if ~isempty(intersect(n,[5:1:8]))
     for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 h h 0],'b')
     end
end
hold on
line([rewardLoc rewardLoc],[0 h],'Color','r')
    hold on
    i=find(countAllAvg(n,:)>0);
    bar(xvalues,countAllAvg(n,:),'FaceColor',[0.5 0.5 0.5],'EdgeColor','none');
    hold on
    [p,x]=ksdensity(lickPositionGroup{n},'width',12);
    plot(x,p/max(p)*h)
        set(gca,'xticklabel',[])
    % set(gca,'yticklabel',[])
    xlim([0 400])
end
tightfig
    
saveas(gcf,'RBRLickDistri_GFPRS.fig')
    print -painters -depsc RBRLickDistri_GFPRS.eps
%% plot per session
sessionDataAll={};
load('LickPercentNearStimCHR2CS.mat')
sessionDataAll{1}=[lickPercentileBefore' lickPercentileStim' lickPercentileAfter'];
load('LickPercentNearStimCHR2RS.mat')
sessionDataAll{2}=[lickPercentileBefore' lickPercentileStim' lickPercentileAfter'];
load('LickPercentNearStimGFPCS.mat')
sessionDataAll{3}=[lickPercentileBefore' lickPercentileStim' lickPercentileAfter'];
load('LickPercentNearStimGFPRS.mat')
sessionDataAll{4}=[lickPercentileBefore' lickPercentileStim' lickPercentileAfter'];


figure
for t=1:length(sessionDataAll)
allRBR=sessionDataAll{t};

allSession=[];
allSession(:,1)=nanmean(allRBR(:,[1:10]),2);
allSession(:,2)=nanmean(allRBR(:,[11:20]),2);
allSession(:,3)=nanmean(allRBR(:,[21:30]),2);
s=sum(allSession,2);
idxUse=find(~isnan(s));
allSession=allSession(idxUse,:);
p=[];
[~,p(1)]=ttest(allSession(:,1),allSession(:,2));
[~,p(2)]=ttest(allSession(:,1),allSession(:,3));
[~,p(3)]=ttest(allSession(:,2),allSession(:,3));
[isSignificant,adjusted_pvals,~]= bonferroni_holm(p,0.05); % adjusted_pvals are the corrected p values

subplot(1,4,t)
for i=1:size(allSession,1);
    hold on
    plot([1:1:3],allSession(i,:),'m-')
end

bar([1:1:3],nanmean(allSession,1),'FaceColor','w');
hold on
errorbar([1:1:3],nanmean(allSession,1),nansem(allSession,1),'k.');
ylim([0 1])
if t==1;
title(['CS',num2str(adjusted_pvals)])
elseif t==2;
    title(['RS',num2str(adjusted_pvals)])
elseif t==3;
    title(['CSg',num2str(adjusted_pvals)])
else
        title(['RSg',num2str(adjusted_pvals)])
end

end

saveas(gcf,'LickPercentNearStimPerSession.fig')
%% RBR normalized to 1


RBRCSadj=sessionDataAll{1};

RBRadj=sessionDataAll{2};

RBRCSGFPadj=sessionDataAll{3};

RBRRSGFPadj=sessionDataAll{4};

RBR=RBRadj;
RBRCS=RBRCSadj;
RBRCSGFP=RBRCSGFPadj;
RBRRSGFP=RBRRSGFPadj;
N=30;

figure
subplot(221)
A=RBR(:,1:N);
% for n=1:size(A,1);
    ma=nanmean(A(:,[1:10]),'all');
%     A(n,:)=A(n,:)/ma;
% end
% A(isinf(A))=nan;
A=A/ma;

B=RBRCS(:,1:N);
% for n=1:size(B,1);
    ma=nanmean(B(:,[1:10]),'all');
%     B(n,:)=B(n,:)/ma;
% end
% B(isinf(B))=nan;
B=B/ma;
hold on
errorbar([1:1:30],nanmean(A,1),nansem(A,1),'k')
hold on
errorbar([1:1:30],nanmean(B,1),nansem(B,1),'m')
% semshade(A,0.3,'k',[1:N])
% hold on
% semshade(B,0.3,'m',[1:N])
%

 % [pAnova,pMC,pLabel] = anovaRM2W_full_BH(A(:,11:20),B(:,11:20),1);
 % [pAnova2,pMC,pLabel] = anovaRM2W_full_BH(A(:,1:10),B(:,1:10),1);
 % [pAnova3,pMC,pLabel] = anovaRM2W_full_BH(A(:,21:30),B(:,21:30),1);
[~,p1]=ttest2(reshape(A(:,[11:20]),[1 size(A,1)*10]),reshape(B(:,[11:20]),[1 size(B,1)*10]));
[~,p2]=ttest2(reshape(A(:,[1:10]),[1 size(A,1)*10]),reshape(B(:,[1:10]),[1 size(B,1)*10]));
[~,p3]=ttest2(reshape(A(:,[21:30]),[1 size(A,1)*10]),reshape(B(:,[21:30]),[1 size(B,1)*10]));

% ylim([0 4])
% title(['chr2 CS RS,P=',num2str([pAnova(1) pAnova2(1) pAnova3(1)])])
title(['chr2 CS RS',num2str([p1 p2 p3])])

xlim([0 N])
p1=[];
for n=1:N;
    [r,p1(n)]= ttest2(A(:,n),B(:,n));
end


for n=1:N;
    if p1(n)<=0.05;
        hold on
        plot(n,0.5,'k*')
    end
end

subplot(222)
A=RBRRSGFP(:,1:N);
% for n=1:size(A,1);
    ma=nanmean(A(:,[1:10]),'all');
%     A(n,:)=A(n,:)/ma;
% end
A=A/ma;

B=RBRCSGFP(:,1:N);
% for n=1:size(B,1);
    ma=nanmean(B(:,[1:10]),'all');
%     B(n,:)=B(n,:)/ma;
% end
B=B/ma;
hold on
errorbar([1:1:30],nanmean(A,1),nansem(A,1),'k')
hold on
errorbar([1:1:30],nanmean(B,1),nansem(B,1),'m')
 % [pAnova,pMC,pLabel] = anovaRM2W_full_BH(A(:,11:20),B(:,11:20),1);
 % [pAnova2,pMC,pLabel] = anovaRM2W_full_BH(A(:,1:10),B(:,1:10),1);
 % [pAnova3,pMC,pLabel] = anovaRM2W_full_BH(A(:,21:30),B(:,21:30),1);
% ylim([0 4])
[~,p1]=ttest2(reshape(A(:,[11:20]),[1 size(A,1)*10]),reshape(B(:,[11:20]),[1 size(B,1)*10]));
[~,p2]=ttest2(reshape(A(:,[1:10]),[1 size(A,1)*10]),reshape(B(:,[1:10]),[1 size(B,1)*10]));
[~,p3]=ttest2(reshape(A(:,[21:30]),[1 size(A,1)*10]),reshape(B(:,[21:30]),[1 size(B,1)*10]));


% ylim([0.6 1.3])
% title(['gfp CS RS, P=',num2str([pAnova(1) pAnova2(1) pAnova3(1)])])
title(['gfp CS RS',num2str([p1 p2 p3])])

xlim([0 N])
p2=[];
for n=1:N;
    [r,p2(n)]= ttest2(A(:,n),B(:,n));
end
for n=1:N;
    if p2(n)<=0.05;
        hold on
        plot(n,0.5,'k*')
    end
end

%

subplot(224)
A=RBRCS(:,1:N);
% for n=1:size(A,1);
    ma=nanmean(A(:,[1:10]),'all');
%     A(n,:)=A(n,:)/ma;
% end

A=A/ma;

B=RBRCSGFP(:,1:N);
% for n=1:size(B,1);
    ma=nanmean(B(:,[1:10]),'all');
%     B(n,:)=B(n,:)/ma;
% end
B=B/ma;
hold on
errorbar([1:1:30],nanmean(A,1),nansem(A,1),'k')
hold on
errorbar([1:1:30],nanmean(B,1),nansem(B,1),'m')
 % [pAnova,pMC,pLabel] = anovaRM2W_full_BH(A(:,11:20),B(:,11:20),1);
 % [pAnova2,pMC,pLabel] = anovaRM2W_full_BH(A(:,1:10),B(:,1:10),1);
 % [pAnova3,pMC,pLabel] = anovaRM2W_full_BH(A(:,21:30),B(:,21:30),1);

[~,p1]=ttest2(reshape(A(:,[11:20]),[1 size(A,1)*10]),reshape(B(:,[11:20]),[1 size(B,1)*10]));
[~,p2]=ttest2(reshape(A(:,[1:10]),[1 size(A,1)*10]),reshape(B(:,[1:10]),[1 size(B,1)*10]));
[~,p3]=ttest2(reshape(A(:,[21:30]),[1 size(A,1)*10]),reshape(B(:,[21:30]),[1 size(B,1)*10]));


% ylim([0.6 1.3])
% title(['CS CHR2 GFP, P=',num2str([pAnova(1) pAnova2(1) pAnova3(1)])])
title(['CS CHR2 GFP',num2str([p1 p2 p3])])

xlim([0 N])
p3=[];
for n=1:N;
    [r,p3(n)]= ttest2(A(:,n),B(:,n));
end
for n=1:N;
    if p3(n)<=0.05;
        hold on
        plot(n,0.5,'k*')
    end
end



subplot(223)
A=RBR(:,1:N);
% for n=1:size(A,1);
    ma=nanmean(A(:,[1:10]),'all');
%     A(n,:)=A(n,:)/ma;
% end
A=A/ma;
B=RBRRSGFP(:,1:N);
% for n=1:size(B,1);
    ma=nanmean(B(:,[1:10]),'all');
%     B(n,:)=B(n,:)/ma;
% end
B=B/ma;
hold on
errorbar([1:1:30],nanmean(A,1),nansem(A,1),'k')
hold on
errorbar([1:1:30],nanmean(B,1),nansem(B,1),'m')
 % [pAnova,pMC,pLabel] = anovaRM2W_full_BH(A(:,11:20),B(:,11:20),1);
 %  [pAnova2,pMC,pLabel] = anovaRM2W_full_BH(A(:,1:10),B(:,1:10),1);
 % [pAnova3,pMC,pLabel] = anovaRM2W_full_BH(A(:,21:30),B(:,21:30),1);

[~,p1]=ttest2(reshape(A(:,[11:20]),[1 size(A,1)*10]),reshape(B(:,[11:20]),[1 size(B,1)*10]));
[~,p2]=ttest2(reshape(A(:,[1:10]),[1 size(A,1)*10]),reshape(B(:,[1:10]),[1 size(B,1)*10]));
[~,p3]=ttest2(reshape(A(:,[21:30]),[1 size(A,1)*10]),reshape(B(:,[21:30]),[1 size(B,1)*10]));
% ylim([0.6 1.3])
% title(['RS CHR2 GFP, P=',num2str([pAnova(1) pAnova2(1) pAnova3(1)])])
title(['RS CHR2 GFP',num2str([p1 p2 p3])])
xlim([0 N])
p4=[];
for n=1:N;
    [r,p4(n)]= ttest2(A(:,n),B(:,n));
end
for n=1:N;
    if p4(n)<=0.05;
        hold on
        plot(n,0.5,'k*')
    end
end
saveas(gcf,'LickPercentNearStimRBR.fig')

%% RBR normalized to 0


RBRCSadj=sessionDataAll{1};

RBRadj=sessionDataAll{2};

RBRCSGFPadj=sessionDataAll{3};

RBRRSGFPadj=sessionDataAll{4};

RBR=RBRadj;
RBRCS=RBRCSadj;
RBRCSGFP=RBRCSGFPadj;
RBRRSGFP=RBRRSGFPadj;
N=30;

figure
subplot(221)
A=RBR(:,1:N);
% for n=1:size(A,1);
    ma=nanmean(A(:,[1:10]),'all');
%     A(n,:)=A(n,:)-ma;
% end
% A(isinf(A))=nan;
A=A-ma;

B=RBRCS(:,1:N);
% for n=1:size(B,1);
    ma=nanmean(B(:,[1:10]),'all');
%     B(n,:)=B(n,:)-ma;
% end
% B(isinf(B))=nan;
B=B-ma;
hold on
errorbar([1:1:30],nanmean(A,1),nansem(A,1),'k')
hold on
errorbar([1:1:30],nanmean(B,1),nansem(B,1),'m')
% semshade(A,0.3,'k',[1:N])
% hold on
% semshade(B,0.3,'m',[1:N])
%

 % [pAnova,pMC,pLabel] = anovaRM2W_full_BH(A(:,11:20),B(:,11:20),1);
 % [pAnova2,pMC,pLabel] = anovaRM2W_full_BH(A(:,1:10),B(:,1:10),1);
 % [pAnova3,pMC,pLabel] = anovaRM2W_full_BH(A(:,21:30),B(:,21:30),1);
[~,p1]=ttest2(reshape(A(:,[11:20]),[1 size(A,1)*10]),reshape(B(:,[11:20]),[1 size(B,1)*10]));
[~,p2]=ttest2(reshape(A(:,[1:10]),[1 size(A,1)*10]),reshape(B(:,[1:10]),[1 size(B,1)*10]));
[~,p3]=ttest2(reshape(A(:,[21:30]),[1 size(A,1)*10]),reshape(B(:,[21:30]),[1 size(B,1)*10]));
A2=reshape(A(:,[11:20]),[1 size(A,1)*10]);
A1=reshape(A(:,[1:10]),[1 size(A,1)*10]);

B2=reshape(B(:,[11:20]),[1 size(B,1)*10]);
B1=reshape(B(:,[1:10]),[1 size(B,1)*10]);

[~,pPreDurRS]=ttest2(A2,A1)
[~,pPreDurCS]=ttest2(B2,B1)
% ylim([0 4])
% title(['chr2 CS RS,P=',num2str([pAnova(1) pAnova2(1) pAnova3(1)])])
title(['chr2 CS RS',num2str([p1 p2 p3 pPreDurCS pPreDurRS])])

xlim([0 N])
p1=[];
for n=1:N;
    [r,p1(n)]= ttest2(A(:,n),B(:,n));
end


for n=1:N;
    if p1(n)<=0.05;
        hold on
        plot(n,0.5,'k*')
    end
end

subplot(222)
A=RBRRSGFP(:,1:N);
% for n=1:size(A,1);
    ma=nanmean(A(:,[1:10]),'all');
%     A(n,:)=A(n,:)-ma;
% end
A=A-ma;

B=RBRCSGFP(:,1:N);
% for n=1:size(B,1);
    ma=nanmean(B(:,[1:10]),'all');
%     B(n,:)=B(n,:)-ma;
% end
B=B-ma;
hold on
errorbar([1:1:30],nanmean(A,1),nansem(A,1),'k')
hold on
errorbar([1:1:30],nanmean(B,1),nansem(B,1),'m')
 % [pAnova,pMC,pLabel] = anovaRM2W_full_BH(A(:,11:20),B(:,11:20),1);
 % [pAnova2,pMC,pLabel] = anovaRM2W_full_BH(A(:,1:10),B(:,1:10),1);
 % [pAnova3,pMC,pLabel] = anovaRM2W_full_BH(A(:,21:30),B(:,21:30),1);
% ylim([0 4])
[~,p1]=ttest2(reshape(A(:,[11:20]),[1 size(A,1)*10]),reshape(B(:,[11:20]),[1 size(B,1)*10]));
[~,p2]=ttest2(reshape(A(:,[1:10]),[1 size(A,1)*10]),reshape(B(:,[1:10]),[1 size(B,1)*10]));
[~,p3]=ttest2(reshape(A(:,[21:30]),[1 size(A,1)*10]),reshape(B(:,[21:30]),[1 size(B,1)*10]));
A2=reshape(A(:,[11:20]),[1 size(A,1)*10]);
A1=reshape(A(:,[1:10]),[1 size(A,1)*10]);

B2=reshape(B(:,[11:20]),[1 size(B,1)*10]);
B1=reshape(B(:,[1:10]),[1 size(B,1)*10]);

[~,pPreDurRS]=ttest2(A2,A1)
[~,pPreDurCS]=ttest2(B2,B1)

% ylim([0.6 1.3])
% title(['gfp CS RS, P=',num2str([pAnova(1) pAnova2(1) pAnova3(1)])])
title(['gfp CS RS',num2str([p1 p2 p3 pPreDurCS pPreDurRS])])

xlim([0 N])
p2=[];
for n=1:N;
    [r,p2(n)]= ttest2(A(:,n),B(:,n));
end
for n=1:N;
    if p2(n)<=0.05;
        hold on
        plot(n,0.5,'k*')
    end
end

%

subplot(224)
A=RBRCS(:,1:N);
% for n=1:size(A,1);
    ma=nanmean(A(:,[1:10]),'all');
%     A(n,:)=A(n,:)-ma;
% end

A=A-ma;

B=RBRCSGFP(:,1:N);
% for n=1:size(B,1);
    ma=nanmean(B(:,[1:10]),'all');
%     B(n,:)=B(n,:)-ma;
% end
B=B-ma;
hold on
errorbar([1:1:30],nanmean(A,1),nansem(A,1),'k')
hold on
errorbar([1:1:30],nanmean(B,1),nansem(B,1),'m')
 % [pAnova,pMC,pLabel] = anovaRM2W_full_BH(A(:,11:20),B(:,11:20),1);
 % [pAnova2,pMC,pLabel] = anovaRM2W_full_BH(A(:,1:10),B(:,1:10),1);
 % [pAnova3,pMC,pLabel] = anovaRM2W_full_BH(A(:,21:30),B(:,21:30),1);

[~,p1]=ttest2(reshape(A(:,[11:20]),[1 size(A,1)*10]),reshape(B(:,[11:20]),[1 size(B,1)*10]));
[~,p2]=ttest2(reshape(A(:,[1:10]),[1 size(A,1)*10]),reshape(B(:,[1:10]),[1 size(B,1)*10]));
[~,p3]=ttest2(reshape(A(:,[21:30]),[1 size(A,1)*10]),reshape(B(:,[21:30]),[1 size(B,1)*10]));


% ylim([0.6 1.3])
% title(['CS CHR2 GFP, P=',num2str([pAnova(1) pAnova2(1) pAnova3(1)])])
title(['CS CHR2 GFP',num2str([p1 p2 p3])])

xlim([0 N])
p3=[];
for n=1:N;
    [r,p3(n)]= ttest2(A(:,n),B(:,n));
end
for n=1:N;
    if p3(n)<=0.05;
        hold on
        plot(n,0.5,'k*')
    end
end



subplot(223)
A=RBR(:,1:N);
% for n=1:size(A,1);
    ma=nanmean(A(:,[1:10]),'all');
%     A(n,:)=A(n,:)-ma;
% end
A=A-ma;
B=RBRRSGFP(:,1:N);
% for n=1:size(B,1);
    ma=nanmean(B(:,[1:10]),'all');
%     B(n,:)=B(n,:)-ma;
% end
B=B-ma;
hold on
errorbar([1:1:30],nanmean(A,1),nansem(A,1),'k')
hold on
errorbar([1:1:30],nanmean(B,1),nansem(B,1),'m')
 % [pAnova,pMC,pLabel] = anovaRM2W_full_BH(A(:,11:20),B(:,11:20),1);
 %  [pAnova2,pMC,pLabel] = anovaRM2W_full_BH(A(:,1:10),B(:,1:10),1);
 % [pAnova3,pMC,pLabel] = anovaRM2W_full_BH(A(:,21:30),B(:,21:30),1);

[~,p1]=ttest2(reshape(A(:,[11:20]),[1 size(A,1)*10]),reshape(B(:,[11:20]),[1 size(B,1)*10]));
[~,p2]=ttest2(reshape(A(:,[1:10]),[1 size(A,1)*10]),reshape(B(:,[1:10]),[1 size(B,1)*10]));
[~,p3]=ttest2(reshape(A(:,[21:30]),[1 size(A,1)*10]),reshape(B(:,[21:30]),[1 size(B,1)*10]));
% ylim([0.6 1.3])
% title(['RS CHR2 GFP, P=',num2str([pAnova(1) pAnova2(1) pAnova3(1)])])
title(['RS CHR2 GFP',num2str([p1 p2 p3])])
xlim([0 N])
p4=[];
for n=1:N;
    [r,p4(n)]= ttest2(A(:,n),B(:,n));
end
for n=1:N;
    if p4(n)<=0.05;
        hold on
        plot(n,0.5,'k*')
    end
end
saveas(gcf,'LickPercentNearStimRBRZero.fig')

%% absolute lick at stim zone
%expand stim zone
stimLoc=[];
stimLoc(1,:)=[63+3.5 63+3.5+15];
stimLoc(2,:)=[163+3.5 163+3.5+15];
stimLoc(3,:)=[279+3.5 279+3.5+15];
stimExpThis=50;
stimLoc(:,1)=stimLoc(:,1)-stimExpThis;
stimLoc(:,2)=stimLoc(:,2)+stimExpThis;
load('lickPositionCHR2CS.mat');
%number of licks in stim zone
nLickStimZoneBefore=[];%each row is a session, each column is a run
nLickStimZoneStim=[];
nLickStimZoneAfter=[];

for session=1:length(lickUsePositionBeforeRBR)
    for lap=1:length(lickUsePositionBeforeRBR{session});
        runloc=lickUsePositionBeforeRBR{session}{lap};
        nRun=[];
        for s=1:size(stimLoc,1);
            nRun(s)=length(find(runloc>=stimLoc(s,1)&runloc<=stimLoc(s,2)));
        end
        nRun=nansum(nRun);
        nLickStimZoneBefore(session,lap)=nRun;

        runloc=lickUsePositionStimRBR{session}{lap};
        nRun=[];
        for s=1:size(stimLoc,1);
            nRun(s)=length(find(runloc>=stimLoc(s,1)&runloc<=stimLoc(s,2)));
        end
        nRun=nansum(nRun);
        nLickStimZoneStim(session,lap)=nRun;

        runloc=lickUsePositionAfterRBR{session}{lap};
        nRun=[];
        for s=1:size(stimLoc,1);
            nRun(s)=length(find(runloc>=stimLoc(s,1)&runloc<=stimLoc(s,2)));
        end
        nRun=nansum(nRun);
        nLickStimZoneAfter(session,lap)=nRun;
    end
end
        
RBRCS=[nLickStimZoneBefore nLickStimZoneStim nLickStimZoneAfter];
%%
load('lickPositionCHR2RS.mat');
%number of licks in stim zone
nLickStimZoneBefore=[];%each row is a session, each column is a run
nLickStimZoneStim=[];
nLickStimZoneAfter=[];

for session=1:length(lickUsePositionBeforeRBR)
    for lap=1:length(lickUsePositionBeforeRBR{session});
        runloc=lickUsePositionBeforeRBR{session}{lap};
        nRun=[];
        for s=1:size(stimLoc,1);
            nRun(s)=length(find(runloc>=stimLoc(s,1)&runloc<=stimLoc(s,2)));
        end
        nRun=nansum(nRun);
        nLickStimZoneBefore(session,lap)=nRun;

        runloc=lickUsePositionStimRBR{session}{lap};
        nRun=[];
        for s=1:size(stimLoc,1);
            nRun(s)=length(find(runloc>=stimLoc(s,1)&runloc<=stimLoc(s,2)));
        end
        nRun=nansum(nRun);
        nLickStimZoneStim(session,lap)=nRun;

        runloc=lickUsePositionAfterRBR{session}{lap};
        nRun=[];
        for s=1:size(stimLoc,1);
            nRun(s)=length(find(runloc>=stimLoc(s,1)&runloc<=stimLoc(s,2)));
        end
        nRun=nansum(nRun);
        nLickStimZoneAfter(session,lap)=nRun;
    end
end
        
RBR=[nLickStimZoneBefore nLickStimZoneStim nLickStimZoneAfter];

%%
load('lickPositionGFPCS.mat');
%number of licks in stim zone
nLickStimZoneBefore=[];%each row is a session, each column is a run
nLickStimZoneStim=[];
nLickStimZoneAfter=[];

for session=1:length(lickUsePositionBeforeRBR)
    for lap=1:length(lickUsePositionBeforeRBR{session});
        runloc=lickUsePositionBeforeRBR{session}{lap};
        nRun=[];
        for s=1:size(stimLoc,1);
            nRun(s)=length(find(runloc>=stimLoc(s,1)&runloc<=stimLoc(s,2)));
        end
        nRun=nansum(nRun);
        nLickStimZoneBefore(session,lap)=nRun;

        runloc=lickUsePositionStimRBR{session}{lap};
        nRun=[];
        for s=1:size(stimLoc,1);
            nRun(s)=length(find(runloc>=stimLoc(s,1)&runloc<=stimLoc(s,2)));
        end
        nRun=nansum(nRun);
        nLickStimZoneStim(session,lap)=nRun;

        runloc=lickUsePositionAfterRBR{session}{lap};
        nRun=[];
        for s=1:size(stimLoc,1);
            nRun(s)=length(find(runloc>=stimLoc(s,1)&runloc<=stimLoc(s,2)));
        end
        nRun=nansum(nRun);
        nLickStimZoneAfter(session,lap)=nRun;
    end
end
        
RBRCSGFP=[nLickStimZoneBefore nLickStimZoneStim nLickStimZoneAfter];

%%
load('lickPositionGFPRS.mat');
%number of licks in stim zone
nLickStimZoneBefore=[];%each row is a session, each column is a run
nLickStimZoneStim=[];
nLickStimZoneAfter=[];

for session=1:length(lickUsePositionBeforeRBR)
    for lap=1:length(lickUsePositionBeforeRBR{session});
        runloc=lickUsePositionBeforeRBR{session}{lap};
        nRun=[];
        for s=1:size(stimLoc,1);
            nRun(s)=length(find(runloc>=stimLoc(s,1)&runloc<=stimLoc(s,2)));
        end
        nRun=nansum(nRun);
        nLickStimZoneBefore(session,lap)=nRun;

        runloc=lickUsePositionStimRBR{session}{lap};
        nRun=[];
        for s=1:size(stimLoc,1);
            nRun(s)=length(find(runloc>=stimLoc(s,1)&runloc<=stimLoc(s,2)));
        end
        nRun=nansum(nRun);
        nLickStimZoneStim(session,lap)=nRun;

        runloc=lickUsePositionAfterRBR{session}{lap};
        nRun=[];
        for s=1:size(stimLoc,1);
            nRun(s)=length(find(runloc>=stimLoc(s,1)&runloc<=stimLoc(s,2)));
        end
        nRun=nansum(nRun);
        nLickStimZoneAfter(session,lap)=nRun;
    end
end
        
RBRRSGFP=[nLickStimZoneBefore nLickStimZoneStim nLickStimZoneAfter];

%% %absolute lick at stim zone

N=30;

figure
subplot(221)
A=RBR(:,1:N);
% for n=1:size(A,1);
    ma=nanmean(A(:,[1:10]),'all');
%     A(n,:)=A(n,:)/ma;
% end
% A(isinf(A))=nan;
A=A/ma;

B=RBRCS(:,1:N);
% for n=1:size(B,1);
    ma=nanmean(B(:,[1:10]),'all');
%     B(n,:)=B(n,:)/ma;
% end
% B(isinf(B))=nan;
B=B/ma;
hold on
errorbar([1:1:30],nanmean(A,1),nansem(A,1),'k')
hold on
errorbar([1:1:30],nanmean(B,1),nansem(B,1),'m')
% semshade(A,0.3,'k',[1:N])
% hold on
% semshade(B,0.3,'m',[1:N])
%

 % [pAnova,pMC,pLabel] = anovaRM2W_full_BH(A(:,11:20),B(:,11:20),1);
 % [pAnova2,pMC,pLabel] = anovaRM2W_full_BH(A(:,1:10),B(:,1:10),1);
 % [pAnova3,pMC,pLabel] = anovaRM2W_full_BH(A(:,21:30),B(:,21:30),1);


% ylim([0 4])
% title(['chr2 CS RS,P=',num2str([pAnova(1) pAnova2(1) pAnova3(1)])])
title(['chr2 CS RS'])

xlim([0 N])
p1=[];
for n=1:N;
    [r,p1(n)]= ttest2(A(:,n),B(:,n));
end


for n=1:N;
    if p1(n)<=0.05;
        hold on
        plot(n,0.5,'k*')
    end
end

subplot(222)
A=RBRRSGFP(:,1:N);
% for n=1:size(A,1);
    ma=nanmean(A(:,[1:10]),'all');
%     A(n,:)=A(n,:)/ma;
% end
A=A/ma;

B=RBRCSGFP(:,1:N);
% for n=1:size(B,1);
    ma=nanmean(B(:,[1:10]),'all');
%     B(n,:)=B(n,:)/ma;
% end
B=B/ma;
hold on
errorbar([1:1:30],nanmean(A,1),nansem(A,1),'k')
hold on
errorbar([1:1:30],nanmean(B,1),nansem(B,1),'m')
 % [pAnova,pMC,pLabel] = anovaRM2W_full_BH(A(:,11:20),B(:,11:20),1);
 % [pAnova2,pMC,pLabel] = anovaRM2W_full_BH(A(:,1:10),B(:,1:10),1);
 % [pAnova3,pMC,pLabel] = anovaRM2W_full_BH(A(:,21:30),B(:,21:30),1);
% ylim([0 4])

% ylim([0.6 1.3])
% title(['gfp CS RS, P=',num2str([pAnova(1) pAnova2(1) pAnova3(1)])])
title(['gfp CS RS']);

xlim([0 N])
p2=[];
for n=1:N;
    [r,p2(n)]= ttest2(A(:,n),B(:,n));
end
for n=1:N;
    if p2(n)<=0.05;
        hold on
        plot(n,0.5,'k*')
    end
end

%

subplot(224)
A=RBRCS(:,1:N);
% for n=1:size(A,1);
    ma=nanmean(A(:,[1:10]),'all');
%     A(n,:)=A(n,:)/ma;
% end

A=A/ma;

B=RBRCSGFP(:,1:N);
% for n=1:size(B,1);
    ma=nanmean(B(:,[1:10]),'all');
%     B(n,:)=B(n,:)/ma;
% end
B=B/ma;
hold on
errorbar([1:1:30],nanmean(A,1),nansem(A,1),'k')
hold on
errorbar([1:1:30],nanmean(B,1),nansem(B,1),'m')
 % [pAnova,pMC,pLabel] = anovaRM2W_full_BH(A(:,11:20),B(:,11:20),1);
 % [pAnova2,pMC,pLabel] = anovaRM2W_full_BH(A(:,1:10),B(:,1:10),1);
 % [pAnova3,pMC,pLabel] = anovaRM2W_full_BH(A(:,21:30),B(:,21:30),1);


% ylim([0.6 1.3])
% title(['CS CHR2 GFP, P=',num2str([pAnova(1) pAnova2(1) pAnova3(1)])])
title(['CS CHR2 GFP'])

xlim([0 N])
p3=[];
for n=1:N;
    [r,p3(n)]= ttest2(A(:,n),B(:,n));
end
for n=1:N;
    if p3(n)<=0.05;
        hold on
        plot(n,0.5,'k*')
    end
end



subplot(223)
A=RBR(:,1:N);
% for n=1:size(A,1);
    ma=nanmean(A(:,[1:10]),'all');
%     A(n,:)=A(n,:)/ma;
% end
A=A/ma;
B=RBRRSGFP(:,1:N);
% for n=1:size(B,1);
    ma=nanmean(B(:,[1:10]),'all');
%     B(n,:)=B(n,:)/ma;
% end
B=B/ma;
hold on
errorbar([1:1:30],nanmean(A,1),nansem(A,1),'k')
hold on
errorbar([1:1:30],nanmean(B,1),nansem(B,1),'m')
 % [pAnova,pMC,pLabel] = anovaRM2W_full_BH(A(:,11:20),B(:,11:20),1);
 %  [pAnova2,pMC,pLabel] = anovaRM2W_full_BH(A(:,1:10),B(:,1:10),1);
 % [pAnova3,pMC,pLabel] = anovaRM2W_full_BH(A(:,21:30),B(:,21:30),1);


% ylim([0.6 1.3])
% title(['RS CHR2 GFP, P=',num2str([pAnova(1) pAnova2(1) pAnova3(1)])])
title(['RS CHR2 GFP'])
xlim([0 N])
p4=[];
for n=1:N;
    [r,p4(n)]= ttest2(A(:,n),B(:,n));
end
for n=1:N;
    if p4(n)<=0.05;
        hold on
        plot(n,0.5,'k*')
    end
end
saveas(gcf,'LickAbsoluteWithinStimRBR.fig')

%% %development of pattern with stim pattern
track=zeros(400,1);
stimLoc=[];
stimLoc(1,:)=[63+3.5 63+3.5+15];
stimLoc(2,:)=[163+3.5 163+3.5+15];
stimLoc(3,:)=[279+3.5 279+3.5+15];
% stimExp1=0;
% stimExp2=10;
% stimExp1=10;
% stimExp1=5;
% stimExp2=50;
%use if bin width =5; track is 1:400/binWidth
% stimExp1=5;
% stimExp2=50;
%use if bin width =1; track is 16:400/binWidth
stimExp1=0;
stimExp2=25;

% stimExp1=0;
% stimExp2=20;
stimLoc(:,1)=stimLoc(:,1)-stimExp1;
stimLoc(:,2)=stimLoc(:,2)+stimExp2;

for n=1:size(stimLoc,1)
    track(ceil(stimLoc(n,1)):ceil(stimLoc(n,2)))=1;
end
%bin the track
% binWidth=5; %this could work too but cannot speficically do 16cm to 400
% to be consistent with here:
% Z:\papers\2025RewardOptogenetics\speedExamples\cue3_15\plotForPaper: the
% beginnign of the track: the mice made stops to lick for reward
binWidth=1;
% binWidth=5;
binEdges=[0:binWidth:length(track)];
binCenters=[binEdges(1:end-1)+binWidth/2];
%group cells in every bin
binnedTrack=[];
for n=1:length(track)/binWidth;
    binnedTrack(n,1)=nanmean(track(binWidth*(n-1)+1:binWidth*n));
end

% useRegion=[5:400/binWidth];
% useRegion=[1:400/binWidth]; 
useRegion=[16:400/binWidth];
% useRegion=[1:400/binWidth];
rollingWindow=4;


% chr2 cs
CHR2CSCorr=[];%EACH row is a run, each column is a session
clear lickUsePositionBeforeRBR
clear lickUsePositionAfterRBR
clear lickUsePositionStimRBR
load('lickPositionCHR2CS.mat');


beforeCorr=[];
B=[];
A=lickUsePositionBeforeRBR;
for session=1:length(A);
    for lap=1:length(A{session})-rollingWindow+1; %rolling average every 2
        lick=cell2mat(A{session}(lap:lap+rollingWindow-1)');
        if isempty(lick);
            B(lap,session)=nan;
        else
            [p,x]=ksdensity(lick,binCenters);
            B(lap,session)=corr(p(useRegion)',binnedTrack(useRegion));
        end
    end
end
beforeCorr=B;

stimCorr=[];
B=[];
A=lickUsePositionStimRBR;
for session=1:length(A);
    for lap=1:length(A{session})-rollingWindow+1; %rolling average every 2
        lick=cell2mat(A{session}(lap:lap+rollingWindow-1)');
        if isempty(lick);
            B(lap,session)=nan;
        else
            [p,x]=ksdensity(lick,binCenters);
            B(lap,session)=corr(p(useRegion)',binnedTrack(useRegion));
        end
    end
end
stimCorr=B;

afterCorr=[];
B=[];
A=lickUsePositionAfterRBR;
for session=1:length(A);
    for lap=1:length(A{session})-rollingWindow+1; %rolling average every 2
        lick=cell2mat(A{session}(lap:lap+rollingWindow-1)');
        if isempty(lick);
            B(lap,session)=nan;
        else
            [p,x]=ksdensity(lick,binCenters);
            B(lap,session)=corr(p(useRegion)',binnedTrack(useRegion));
        end
    end
end
afterCorr=B;

CHR2CSCorr=[beforeCorr;stimCorr;afterCorr];

% chr2 rs
CHR2RSCorr=[];%EACH row is a run, each column is a session
clear lickUsePositionBeforeRBR
clear lickUsePositionAfterRBR
clear lickUsePositionStimRBR

load('lickPositionCHR2RS.mat');
beforeCorr=[];
B=[];
A=lickUsePositionBeforeRBR;
for session=1:length(A);
    for lap=1:length(A{session})-rollingWindow+1; %rolling average every 2
        lick=cell2mat(A{session}(lap:lap+rollingWindow-1)');
        if isempty(lick);
            B(lap,session)=nan;
        else
            [p,x]=ksdensity(lick,binCenters);
            B(lap,session)=corr(p(useRegion)',binnedTrack(useRegion));
        end
    end
end
beforeCorr=B;

stimCorr=[];
B=[];
A=lickUsePositionStimRBR;
for session=1:length(A);
    for lap=1:length(A{session})-rollingWindow+1; %rolling average every 2
        lick=cell2mat(A{session}(lap:lap+rollingWindow-1)');
        if isempty(lick);
            B(lap,session)=nan;
        else
            [p,x]=ksdensity(lick,binCenters);
            B(lap,session)=corr(p(useRegion)',binnedTrack(useRegion));
        end
    end
end
stimCorr=B;

afterCorr=[];
B=[];
A=lickUsePositionAfterRBR;
for session=1:length(A);
    for lap=1:length(A{session})-rollingWindow+1; %rolling average every 2
        lick=cell2mat(A{session}(lap:lap+rollingWindow-1)');
        if isempty(lick);
            B(lap,session)=nan;
        else
            [p,x]=ksdensity(lick,binCenters);
            B(lap,session)=corr(p(useRegion)',binnedTrack(useRegion));
        end
    end
end
afterCorr=B;

CHR2RSCorr=[beforeCorr;stimCorr;afterCorr];

% GFP cs
GFPCSCorr=[];%EACH row is a run, each column is a session
clear lickUsePositionBeforeRBR
clear lickUsePositionAfterRBR
clear lickUsePositionStimRBR

load('lickPositionGFPCS.mat');
beforeCorr=[];
B=[];
A=lickUsePositionBeforeRBR;
for session=1:length(A);
    for lap=1:length(A{session})-rollingWindow+1; %rolling average every 2
        lick=cell2mat(A{session}(lap:lap+rollingWindow-1)');
        if isempty(lick);
            B(lap,session)=nan;
        else
            [p,x]=ksdensity(lick,binCenters);
            B(lap,session)=corr(p(useRegion)',binnedTrack(useRegion));
        end
    end
end
beforeCorr=B;

stimCorr=[];
B=[];
A=lickUsePositionStimRBR;
for session=1:length(A);
    for lap=1:length(A{session})-rollingWindow+1; %rolling average every 2
        lick=cell2mat(A{session}(lap:lap+rollingWindow-1)');
        if isempty(lick);
            B(lap,session)=nan;
        else
            [p,x]=ksdensity(lick,binCenters);
            B(lap,session)=corr(p(useRegion)',binnedTrack(useRegion));
        end
    end
end
stimCorr=B;

afterCorr=[];
B=[];
A=lickUsePositionAfterRBR;
for session=1:length(A);
    for lap=1:length(A{session})-rollingWindow+1; %rolling average every 2
        lick=cell2mat(A{session}(lap:lap+rollingWindow-1)');
        if isempty(lick);
            B(lap,session)=nan;
        else
            [p,x]=ksdensity(lick,binCenters);
            B(lap,session)=corr(p(useRegion)',binnedTrack(useRegion));
        end
    end
end
afterCorr=B;

GFPCSCorr=[beforeCorr;stimCorr;afterCorr];

% chr2 rs
GFPRSCorr=[];%EACH row is a run, each column is a session
clear lickUsePositionBeforeRBR
clear lickUsePositionAfterRBR
clear lickUsePositionStimRBR

load('lickPositionGFPRS.mat');
beforeCorr=[];
B=[];
A=lickUsePositionBeforeRBR;
for session=1:length(A);
    for lap=1:length(A{session})-rollingWindow+1; %rolling average every 2
        lick=cell2mat(A{session}(lap:lap+rollingWindow-1)');
        if isempty(lick);
            B(lap,session)=nan;
        else
            [p,x]=ksdensity(lick,binCenters);
            B(lap,session)=corr(p(useRegion)',binnedTrack(useRegion));
        end
    end
end
beforeCorr=B;

stimCorr=[];
B=[];
A=lickUsePositionStimRBR;
for session=1:length(A);
    for lap=1:length(A{session})-rollingWindow+1; %rolling average every 2
        lick=cell2mat(A{session}(lap:lap+rollingWindow-1)');
        if isempty(lick);
            B(lap,session)=nan;
        else
            [p,x]=ksdensity(lick,binCenters);
            B(lap,session)=corr(p(useRegion)',binnedTrack(useRegion));
        end
    end
end
stimCorr=B;

afterCorr=[];
B=[];
A=lickUsePositionAfterRBR;
for session=1:length(A);
    for lap=1:length(A{session})-rollingWindow+1; %rolling average every 2
        lick=cell2mat(A{session}(lap:lap+rollingWindow-1)');
        if isempty(lick);
            B(lap,session)=nan;
        else
            [p,x]=ksdensity(lick,binCenters);
            B(lap,session)=corr(p(useRegion)',binnedTrack(useRegion));
        end
    end
end
afterCorr=B;

GFPRSCorr=[beforeCorr;stimCorr;afterCorr];

corrAll={};
corrAll{1}=CHR2CSCorr;
corrAll{2}=CHR2RSCorr;
corrAll{3}=GFPCSCorr;
corrAll{4}=GFPRSCorr;
%
figure
for n=1:length(corrAll);
    subplot(3,2,n);

    if n==1 | n==3;

    errorbar([1:1:size(corrAll{1},1)],nanmean(corrAll{n},2),nansem(corrAll{n},2),'m');
    else
        errorbar([1:1:size(corrAll{1},1)],nanmean(corrAll{n},2),nansem(corrAll{n},2),'k');
    end
hold on
line([1 size(corrAll{1},1)],[0 0])
    pValue=[];

    for lap=1:size(corrAll{n},1);
        [~,pValue(lap)]=ttest(corrAll{n}(lap,:),0);
    end

    i=find(pValue<=0.05);
    hold on
    plot(i,ones(1,length(i))*0.35,'r*')    
    if n==1 | n==2;
        ylim([-0.2 0.4])
    else
        ylim([-0.25 0.4])
    end

   
   if n==1;
    title('chr2cs')
   elseif n==2;
       title('chr2rs')
   elseif n==3;
title('gfpcs')
   else
       title('gfprs');
   end
end

%plot example
subplot(3,2,5)

load('lickPositionCHR2CS.mat');
A=lickUsePositionStimRBR;
for session=4;
    for lap=3; %rolling average every 2
        lick=cell2mat(A{session}(lap:lap+rollingWindow-1)');
        if isempty(lick);
            B(lap,session)=nan;
        else
       
            [p,x]=ksdensity(lick,binCenters);
        end
    end
end

h=histogram(lick,20);
a=h.Values;
x=h.BinEdges(1:length(a))+h.BinWidth/2;
a=a/max(a)*max(p);

plot(binnedTrack(useRegion)*max(p))
hold on
plot(p(useRegion),'m')
hold on
bar(x,a,'FaceColor','m')

% n=1;
%  pValue=[];
% 
%     for lap=1:size(corrAll{n},1);
%         [~,pValue(lap)]=ttest(corrAll{n}(lap,:),0);
%     end
%     pValue(7:12)

    saveas(gcf,'lickPatternCorrStimTemp.fig')


%% pattern correlation: include RS patterns


%Get stim zones

%stim zones of RS data
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\stimZoneIdxAllDataFix.mat');
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\positionSegAllData.mat');
%these data order: CHR2R2, GRPRS, CHR2CS, GRPCS
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\idxCHR2CS.mat')
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\idxCHR2RS.mat')
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\idxGFPCS.mat')
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\idxGFPRS.mat')
%THESE idx oder: CHR2CS, CHR2RS, GFPCS, GFPRS
allIdxAll={};
allIdxAll{1}=idxCHR2CS;
allIdxAll{2}=idxCHR2RS;
allIdxAll{3}=idxGFPCS;
allIdxAll{4}=idxGFPRS;

positionAll={};
for N=1:length(allIdxAll);
    positionAll{N}={};
mouseIdx=allIdxAll{N};

position=positionSegAllData(mouseIdx);
stim=stimZoneIdxAllDataFix(mouseIdx);

for n=1:length(position); %per session
    for m=1:length(position{n}) %per run
        positionAll{N}{n,m}=position{n}{m}(stim{n}{m});
    end
end
end

%

track=zeros(400,1);
% stimLoc=[];
% stimLoc(1,:)=[63+3.5 63+3.5+15];
% stimLoc(2,:)=[163+3.5 163+3.5+15];
% stimLoc(3,:)=[279+3.5 279+3.5+15];
% stimExp1=0;
% stimExp2=10;
% stimExp1=10;
% stimExp1=5;
% stimExp2=50;
%use if bin width =5; track is 1:400/binWidth
% stimExp1=5;
% stimExp2=50;
%use if bin width =1; track is 16:400/binWidth
stimExp1=12;
stimExp2=25;
% stimExp1=12;
% stimExp2=12;
% stimExp1=0;
% stimExp2=25;
% stimLoc(:,1)=stimLoc(:,1)-stimExp1;
% stimLoc(:,2)=stimLoc(:,2)+stimExp2;

% for n=1:size(stimLoc,1)
%     track(ceil(stimLoc(n,1)):ceil(stimLoc(n,2)))=1;
% end
%bin the track
% binWidth=5; %this could work too but cannot speficically do 16cm to 400
% to be consistent with here:
% Z:\papers\2025RewardOptogenetics\speedExamples\cue3_15\plotForPaper: the
% beginnign of the track: the mice made stops to lick for reward
binWidth=1;
% binWidth=5;
binEdges=[0:binWidth:length(track)];
binCenters=[binEdges(1:end-1)+binWidth/2];
%group cells in every bin
% binnedTrack=[];
% for n=1:length(track)/binWidth;
%     binnedTrack(n,1)=nanmean(track(binWidth*(n-1)+1:binWidth*n));
% end

% useRegion=[5:400/binWidth];
% useRegion=[1:400/binWidth]; 
useRegion=[16:400/binWidth];
% useRegion=[1:340/binWidth];
% useRegion=[1:400/binWidth];
rollingWindow=1;

% chr2 cs
CHR2CSCorr=[];%EACH row is a run, each column is a session
% clear lickUsePositionBeforeRBR
% clear lickUsePositionAfterRBR
clear lickUsePositionStimRBR
load('lickPositionCHR2CS.mat');
position=positionAll{1};
stimCorr=[];
B=[];
A=lickUsePositionStimRBR;
for session=1:length(A);
    for lap=1:length(A{session})-rollingWindow+1; %rolling average every 2
        lick=cell2mat(A{session}(lap:lap+rollingWindow-1)');
        positionThis=position(session,(lap:lap+rollingWindow-1))';
        positionThis=cell2mat(positionThis);
        positionThis(:,1)=positionThis(:,1)-stimExp1;
        positionThis(:,2)=positionThis(:,2)+stimExp2;
positionThis(positionThis<=0)=1;
         trackThis=track;
            for i=1:size(positionThis,1);
                trackThis(ceil(positionThis(i,1)):ceil(positionThis(i,2)))=1;
            end

            binnedTrack=[];
for n=1:length(trackThis)/binWidth;
    binnedTrack(n,1)=nanmean(trackThis(binWidth*(n-1)+1:binWidth*n));
end


        if isempty(lick);
            B(lap,session)=nan;
        else
            [p,x]=ksdensity(lick,binCenters);
            B(lap,session)=corr(p(useRegion)',binnedTrack(useRegion));
        end
    end
end
stimCorr=B;

CHR2CSCorr=[stimCorr];

% chr2 rs
CHR2RSCorr=[];%EACH row is a run, each column is a session
clear lickUsePositionStimRBR

load('lickPositionCHR2RS.mat');

position=positionAll{2};
stimCorr=[];
B=[];
A=lickUsePositionStimRBR;
for session=1:length(A);
    for lap=1:length(A{session})-rollingWindow+1; %rolling average every 2
        lick=cell2mat(A{session}(lap:lap+rollingWindow-1)');

                positionThis=position(session,(lap:lap+rollingWindow-1))';
        positionThis=cell2mat(positionThis);
        positionThis(:,1)=positionThis(:,1)-stimExp1;
        positionThis(:,2)=positionThis(:,2)+stimExp2;
        positionThis(positionThis<=0)=1;

         trackThis=track;
            for i=1:size(positionThis,1);
                trackThis(ceil(positionThis(i,1)):ceil(positionThis(i,2)))=1;
            end

            binnedTrack=[];
for n=1:length(trackThis)/binWidth;
    binnedTrack(n,1)=nanmean(trackThis(binWidth*(n-1)+1:binWidth*n));
end
        if isempty(lick);
            B(lap,session)=nan;
        else
            [p,x]=ksdensity(lick,binCenters);
            B(lap,session)=corr(p(useRegion)',binnedTrack(useRegion));
        end
    end
end
stimCorr=B;


CHR2RSCorr=[stimCorr];

% chr2 cs
GFPCSCorr=[];%EACH row is a run, each column is a session

clear lickUsePositionStimRBR

load('lickPositionGFPCS.mat');
position=positionAll{3};
stimCorr=[];
B=[];
A=lickUsePositionStimRBR;
for session=1:length(A);
    for lap=1:length(A{session})-rollingWindow+1; %rolling average every 2
        lick=cell2mat(A{session}(lap:lap+rollingWindow-1)');
                positionThis=position(session,(lap:lap+rollingWindow-1))';
        positionThis=cell2mat(positionThis);
        positionThis(:,1)=positionThis(:,1)-stimExp1;
        positionThis(:,2)=positionThis(:,2)+stimExp2;
positionThis(positionThis<=0)=1;
         trackThis=track;
            for i=1:size(positionThis,1);
                trackThis(ceil(positionThis(i,1)):ceil(positionThis(i,2)))=1;
            end

            binnedTrack=[];
for n=1:length(trackThis)/binWidth;
    binnedTrack(n,1)=nanmean(trackThis(binWidth*(n-1)+1:binWidth*n));
end
        if isempty(lick);
            B(lap,session)=nan;
        else
            [p,x]=ksdensity(lick,binCenters);
            B(lap,session)=corr(p(useRegion)',binnedTrack(useRegion));
        end
    end
end
stimCorr=B;

GFPCSCorr=[stimCorr];

% chr2 rs
GFPRSCorr=[];%EACH row is a run, each column is a session
clear lickUsePositionStimRBR

load('lickPositionGFPRS.mat');
position=positionAll{4};
stimCorr=[];
B=[];
A=lickUsePositionStimRBR;
for session=1:length(A);
    for lap=1:length(A{session})-rollingWindow+1; %rolling average every 2
        lick=cell2mat(A{session}(lap:lap+rollingWindow-1)');
                positionThis=position(session,(lap:lap+rollingWindow-1))';
        positionThis=cell2mat(positionThis);
        positionThis(:,1)=positionThis(:,1)-stimExp1;
        positionThis(:,2)=positionThis(:,2)+stimExp2;
positionThis(positionThis<=0)=1;
         trackThis=track;
            for i=1:size(positionThis,1);
                trackThis(ceil(positionThis(i,1)):ceil(positionThis(i,2)))=1;
            end

            binnedTrack=[];
for n=1:length(trackThis)/binWidth;
    binnedTrack(n,1)=nanmean(trackThis(binWidth*(n-1)+1:binWidth*n));
end
        if isempty(lick);
            B(lap,session)=nan;
        else
            [p,x]=ksdensity(lick,binCenters);
            B(lap,session)=corr(p(useRegion)',binnedTrack(useRegion));
        end
    end
end
stimCorr=B;

GFPRSCorr=[stimCorr];

corrAll={};
corrAll{1}=CHR2CSCorr;
corrAll{2}=CHR2RSCorr;
corrAll{3}=GFPCSCorr;
corrAll{4}=GFPRSCorr;
%
figure
for n=1:length(corrAll);
    subplot(4,2,n);

    if n==1 | n==3;

    errorbar([1:1:size(corrAll{1},1)],nanmean(corrAll{n},2),nansem(corrAll{n},2),'m');
    else
        errorbar([1:1:size(corrAll{1},1)],nanmean(corrAll{n},2),nansem(corrAll{n},2),'k');
    end
hold on
line([1 size(corrAll{1},1)],[0 0])
    pValue=[];

    for lap=1:size(corrAll{n},1);
        [~,pValue(lap)]=ttest(corrAll{n}(lap,:),0);
    end

    i=find(pValue<=0.05);
    hold on
    plot(i,ones(1,length(i))*0.25,'r*')    
    if n==1 | n==2;
        ylim([-0.25 0.35])
    else
        ylim([-0.25 0.35])
    end

   
   if n==1;
    title('chr2cs')
   elseif n==2;
       title('chr2rs')
   elseif n==3;
title('gfpcs')
   else
       title('gfprs');
   end
end

%merge CS RS
subplot(4,2,5)
    n=1 ;
    errorbar([1:1:size(corrAll{1},1)],nanmean(corrAll{n},2),nansem(corrAll{n},2),'m');
 hold on
 n=2;
        errorbar([1:1:size(corrAll{1},1)],nanmean(corrAll{n},2),nansem(corrAll{n},2),'k');
hold on
line([1 size(corrAll{1},1)],[0 0])
    pValue=[];

    for lap=1:size(corrAll{1},1);
        [~,pValue(lap)]=ttest2(corrAll{1}(lap,:),corrAll{2}(lap,:));
    end

    i=find(pValue<=0.05);
    hold on
    plot(i,ones(1,length(i))*0.25,'r*') 
    startRun=4;
    n=1;
    c1=corrAll{n}(startRun:end,:);
    c1=reshape(c1,[size(c1,1)*size(c1,2) 1]);
 n=2;
    c2=corrAll{n}(startRun:end,:);
    c2=reshape(c2,[size(c2,1)*size(c2,2) 1]);
      [~,p]=ttest2(c1,c2);
        ylim([-0.2 0.3])

    title(['chr2cs vs rs',num2str(p)])

    %merge CS RS
subplot(4,2,6)
    n=3 ;
    errorbar([1:1:size(corrAll{1},1)],nanmean(corrAll{n},2),nansem(corrAll{n},2),'m');
 hold on
 n=4;
        errorbar([1:1:size(corrAll{1},1)],nanmean(corrAll{n},2),nansem(corrAll{n},2),'k');
hold on
line([1 size(corrAll{1},1)],[0 0])
    pValue=[];

    for lap=1:size(corrAll{1},1);
        [~,pValue(lap)]=ttest2(corrAll{3}(lap,:),corrAll{4}(lap,:));
    end

    i=find(pValue<=0.05);
    hold on
    plot(i,ones(1,length(i))*0.25,'r*')  
    n=3;
    c1=corrAll{n}(startRun:end,:);
    c1=reshape(c1,[size(c1,1)*size(c1,2) 1]);
 n=4;
    c2=corrAll{n}(startRun:end,:);
    c2=reshape(c2,[size(c2,1)*size(c2,2) 1]);
    [~,p]=ttest2(c1,c2);
  
        ylim([-0.2 0.3])

    title(['gfpcs vs rs',num2str(p)])

    %

%plot example
subplot(4,2,7)

load('lickPositionCHR2RS.mat');
A=lickUsePositionStimRBR;
position=positionAll{2};
for session=4;
    for lap=10; 
        lick=cell2mat(A{session}(lap:lap+rollingWindow-1)');
        positionThis=position(session,(lap:lap+rollingWindow-1))';
        positionThis=cell2mat(positionThis);
        positionThis(:,1)=positionThis(:,1)-stimExp1;
        positionThis(:,2)=positionThis(:,2)+stimExp2;
positionThis(positionThis<=0)=1;
         trackThis=track;
            for i=1:size(positionThis,1);
                trackThis(ceil(positionThis(i,1)):ceil(positionThis(i,2)))=1;
            end

            binnedTrack=[];
for n=1:length(trackThis)/binWidth;
    binnedTrack(n,1)=nanmean(trackThis(binWidth*(n-1)+1:binWidth*n));
end
        if isempty(lick);
            B(lap,session)=nan;
        else
       
            [p,x]=ksdensity(lick,binCenters);
        end
    end
end

h=histogram(lick,20);
a=h.Values;
x=h.BinEdges(1:length(a))+h.BinWidth/2;
a=a/max(a)*max(p);

plot(binnedTrack(useRegion)*max(p))
hold on
plot(p(useRegion),'m')
hold on
bar(x,a,'FaceColor','m')

subplot(4,2,8)

load('lickPositionCHR2RS.mat');
A=lickUsePositionStimRBR;
position=positionAll{2};
for session=4;
    for lap=7; 
        lick=cell2mat(A{session}(lap:lap+rollingWindow-1)');
        positionThis=position(session,(lap:lap+rollingWindow-1))';
        positionThis=cell2mat(positionThis);
        positionThis(:,1)=positionThis(:,1)-stimExp1;
        positionThis(:,2)=positionThis(:,2)+stimExp2;

         trackThis=track;
            for i=1:size(positionThis,1);
                trackThis(ceil(positionThis(i,1)):ceil(positionThis(i,2)))=1;
            end

            binnedTrack=[];
for n=1:length(trackThis)/binWidth;
    binnedTrack(n,1)=nanmean(trackThis(binWidth*(n-1)+1:binWidth*n));
end
        if isempty(lick);
            B(lap,session)=nan;
        else
       
            [p,x]=ksdensity(lick,binCenters);
        end
    end
end

h=histogram(lick,20);
a=h.Values;
x=h.BinEdges(1:length(a))+h.BinWidth/2;
a=a/max(a)*max(p);

plot(binnedTrack(useRegion)*max(p))
hold on
plot(p(useRegion),'m')
hold on
bar(x,a,'FaceColor','m')


% n=1;
%  pValue=[];
% 
%     for lap=1:size(corrAll{n},1);
%         [~,pValue(lap)]=ttest(corrAll{n}(lap,:),0);
%     end
%     pValue(7:12)

    saveas(gcf,'lickPatternCorrStimTemp_forRS.fig')
print -painters -depsc lickPatternCorrStimTemp_forRS.eps
%% alinging the stimulation zone to check licks before and after stim, mostly to check RS
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\stimZoneIdxAllDataFix.mat');
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\positionSegAllData.mat');
%these data order: CHR2R2, GRPRS, CHR2CS, GRPCS
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\idxCHR2CS.mat')
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\idxCHR2RS.mat')
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\idxGFPCS.mat')
load('Z:\papers\2025RewardOptogenetics\stopInOutStimZone\cue15cmStim\idxGFPRS.mat')
%THESE idx oder: CHR2CS, CHR2RS, GFPCS, GFPRS
allIdxAll={};
allIdxAll{1}=idxCHR2CS;
allIdxAll{2}=idxCHR2RS;
allIdxAll{3}=idxGFPCS;
allIdxAll{4}=idxGFPRS;

lickUsePositionStimRBRAll={};
nLick={};
% chr2 CS
clear lickUsePositionStimRBR
load('lickPositionCHR2CS.mat');
lickUsePositionStimRBRAll{1}=lickUsePositionStimRBR;
% chr2 RS
clear lickUsePositionStimRBR
load('lickPositionCHR2RS.mat');
lickUsePositionStimRBRAll{2}=lickUsePositionStimRBR;
% GFP CS
clear lickUsePositionStimRBR
load('lickPositionGFPCS.mat');
lickUsePositionStimRBRAll{3}=lickUsePositionStimRBR;
% GFP RS
clear lickUsePositionStimRBR
load('lickPositionGFPRS.mat');
lickUsePositionStimRBRAll{4}=lickUsePositionStimRBR;

beforeAfterZone=50;
% beforeAfterZone=15;
stimStart=15+beforeAfterZone;
stimEnd=310-beforeAfterZone-15+1;
nLine=11;
% nLine=21;
edgeBefore=linspace(0,beforeAfterZone,nLine);
edgeStim=linspace(beforeAfterZone,beforeAfterZone+15,ceil(nLine/beforeAfterZone*15));
edgeAfter=linspace(beforeAfterZone+15,2*beforeAfterZone+15,nLine);

data={};

for N=1:length(allIdxAll);
mouseIdx=allIdxAll{N};

position=positionSegAllData(mouseIdx);
stim=stimZoneIdxAllDataFix(mouseIdx);
positionAll={}; %each row is one session one run: session 1 run 1-10, session2 run 1-10....
%6 numbers: first stim start, first stim end, second stim start, second
%stim end, third stim start, third stim end.

for n=1:length(position); %per session
    for m=1:length(position{n}) %per run
        positionThis=position{n}{m}(stim{n}{m});
        positionAll{end+1,1}=positionThis;
    end
end

A=lickUsePositionStimRBRAll{N};
lickAll={}; %same as position all: %each cell is one session one run: session 1 run 1-10, session2 run 1-10....
for n=1:length(A); %per session
    for m=1:length(A{n}) %per run
        lickAll{end+1,1}=A{n}{m};
    end
end

lickAllInRange={};%all licks that within the stim range
for n=1:size(positionAll,1);
    thisLick=lickAll{n};
    for m=1:size(positionAll{n},1);
        thisStim=positionAll{n}(m,:);
        if thisStim(1)>stimStart&thisStim(2)<stimEnd;
            %in before zone
            lickAllInRange{end+1,1}=thisLick(thisLick<thisStim(1)&thisLick>thisStim(1)-beforeAfterZone)-(thisStim(1)-beforeAfterZone);
            %in stim zone
            lickAllInRange{end,2}=thisLick(thisLick>thisStim(1)&thisLick<thisStim(2))-(thisStim(1)-beforeAfterZone);
            %in after zone
            lickAllInRange{end,3}=thisLick(thisLick>thisStim(2)&thisLick<thisStim(2)+beforeAfterZone)-(thisStim(1)-beforeAfterZone);
        end
    end
end

nLick{N}=[];
for n=1:size(lickAllInRange,1)
    nLick{N}(n,1)=length(lickAllInRange{n,1})/beforeAfterZone;
    nLick{N}(n,2)=length(lickAllInRange{n,2})/15;
    nLick{N}(n,3)=length(lickAllInRange{n,3})/beforeAfterZone;
end

before=cell2mat(lickAllInRange(:,1));
h=histogram(before,edgeBefore);
countBefore=h.Values/length(lickAllInRange);
xBefore=h.BinEdges(1:end-1)+h.BinWidth/2;

stim=cell2mat(lickAllInRange(:,2));
h=histogram(stim,edgeStim);
countStim=h.Values/length(lickAllInRange);
xStim=h.BinEdges(1:end-1)+h.BinWidth/2;

after=cell2mat(lickAllInRange(:,3));
h=histogram(after,edgeAfter);
countAfter=h.Values/length(lickAllInRange);
xAfter=h.BinEdges(1:end-1)+h.BinWidth/2;

%SAVE THE DATA AND PLOT LATER
data{N}{1}=xBefore;
data{N}{2}=countBefore;
data{N}{3}=xStim;
data{N}{4}=countStim;
data{N}{5}=xAfter;
data{N}{6}=countAfter;

end
close all

figure,
for N=1:length(data);
    subplot(4,2,2*(N-1)+1);
    for m=1:3;
    hold on
    if m==2;
    bar(data{N}{2*(m-1)+1},data{N}{2*m},'FaceColor','b');
    else
        bar(data{N}{2*(m-1)+1},data{N}{2*m},'FaceColor','w');
    end
    ylim([0 1.6])
    end

    subplot(4,2,2*N);
    nData=nLick{N};
    M=nanmean(nData,1);
    E=nansem(nData,1);
    bar([1 2 3],M,'FaceColor','w');
    hold on
    errorbar([1 2 3],M,E,'k.');
    ylim([0 0.2])
    [~,p1]=ttest(nData(:,1),nData(:,2));
[~,p2]=ttest(nData(:,2),nData(:,3));
[~,p3]=ttest(nData(:,1),nData(:,3));
title(num2str([p1 p2 p3]))
end
saveas(gcf,'aligniningStimZone.fig')
print -painters -depsc aligniningStimZone.eps
    






