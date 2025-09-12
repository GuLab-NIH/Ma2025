%this code uses only a few samples chr2RS which had higher lick baseline
%so other samples: chr2cs, gfp cs, gfp rs, are all the same as before. so
%in the 'allSRBRUseRunLowerLick.mat' below, just the second cell (for
%chr2rs) is different
load('Z:\papers\2025RewardOptogenetics\optRBR\cue3_15\allIdxAll.mat')
load('Z:\papers\2025RewardOptogenetics\optRBR\cue3_15\allIdxAllLowerLick.mat');
[~,i,iCHR2CSNew]=intersect(allIdxAllLowerLick{1},allIdxAll{1});
load('allSRBRUseRun.mat')
idx=1;
for n=1:length(allSPreRBRUseRun{idx})
allSPreRBRUseRun{idx}{n}=allSPreRBRUseRun{idx}{n}(iCHR2CSNew,:);
allSPostRBRUseRun{idx}{n}=allSPostRBRUseRun{idx}{n}(iCHR2CSNew,:);
allSDurRBRUseRun{idx}{n}=allSDurRBRUseRun{idx}{n}(iCHR2CSNew,:);
end

save('allSRBRUseRunLowerLick.mat','allSPreRBRUseRun','allSPostRBRUseRun','allSDurRBRUseRun')

load('allSpeedYUse.mat')
allSDurUse{idx}=allSDurUse{idx}(iCHR2CSNew);
allSPostUse{idx}=allSPostUse{idx}(iCHR2CSNew);
allSPreUse{idx}=allSPreUse{idx}(iCHR2CSNew);

allYDurUse{idx}=allYDurUse{idx}(iCHR2CSNew);
allYPostUse{idx}=allYPostUse{idx}(iCHR2CSNew);
allYPreUse{idx}=allYPreUse{idx}(iCHR2CSNew);
save('allSpeedYUseLowerLick.mat','allSPostUse','allYPostUse','allSDurUse','allYDurUse','allSPreUse','allYPreUse');

%%
cues=[63 85;163 185;279 301;374 396];
stimLoc=[];
stimLoc(1,:)=[63+3.5 63+3.5+15];
stimLoc(2,:)=[163+3.5 163+3.5+15];
stimLoc(3,:)=[279+3.5 279+3.5+15];
stimExp1=0;
stimExp2=0;
samplingSpacing=1;


binWidth=1;
rewardLoc=366/binWidth;
xvalues=[1:1:400];

%% 
load('allSRBRUseRunLowerLick.mat')
idx=1; % no need to replot chr2cs
pre=allSPreRBRUseRun{idx};
dur=allSDurRBRUseRun{idx};
post=allSPostRBRUseRun{idx};

speedAll=[];%each row is the mean speed of all sessions
for n=1:length(pre);
    speedAll=[speedAll;nanmean(pre{n},1)];
end
for n=1:length(dur);
    speedAll=[speedAll;nanmean(dur{n},1)];
end
for n=1:length(post);
    speedAll=[speedAll;nanmean(post{n},1)];
end

speedAllSession=[pre';dur';post'];

%
% plot
figure
for n=1:size(speedAll,1);
    subplot(size(speedAll,1),2,1+2*(n-1))
    h=max(speedAll(n,:));

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
    % plot(xvalues,speedAll(n,:),'g');
    semshade(speedAllSession{n},0.3,'g',xvalues);

    set(gca,'xticklabel',[])
    % set(gca,'yticklabel',[])

    xlim([0 400])
end

%
% every 5
nDay=[0 1 3 3 3 1 3 3 3 1 3 3 3];
speedAllAvg=[];
for n=1:size(nDay,2)-1;
    idx=[sum(nDay(1:n),2)+1:1:sum(nDay(1:n+1),2)];
    speedAllAvg(n,:)=nanmean(speedAll(idx,:),1);
end

speedAllAvgSession={};
for n=1:size(nDay,2)-1;
    idx=[sum(nDay(1:n),2)+1:1:sum(nDay(1:n+1),2)];
    speedAllAvgSession{n}=cell2mat(speedAllSession(idx));
end



% plot

for n=1:size(speedAllAvg,1);
    subplot(size(speedAll,1),2,[(2*(n-1)+1)*2:2:(2*(n-1)+2)*2])
    h=max(speedAllAvg(n,:));

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
    i=find(speedAllAvg(n,:)>0);
    % plot(xvalues,speedAllAvg(n,:),'g');
   semshade(speedAllAvgSession{n},0.3,'g',xvalues);
    xlim([0 400])
        set(gca,'xticklabel',[])
    % set(gca,'yticklabel',[])

end
tightfig


saveas(gcf,'RBRDeceDistri_CHR2CSLowerLick.fig')
print -painters -depsc RBRDeceDistri_CHR2CSLowerLick.eps
%
% idx=2;
% pre=allSPreRBRUseRun{idx};
% dur=allSDurRBRUseRun{idx};
% post=allSPostRBRUseRun{idx};
% 
% speedAll=[];%each row is the mean speed of all sessions
% for n=1:length(pre);
%     speedAll=[speedAll;nanmean(pre{n},1)];
% end
% for n=1:length(dur);
%     speedAll=[speedAll;nanmean(dur{n},1)];
% end
% for n=1:length(post);
%     speedAll=[speedAll;nanmean(post{n},1)];
% end
% 
% speedAllSession=[pre';dur';post'];
% 
% %
% % plot
% figure
% for n=1:size(speedAll,1);
%     subplot(size(speedAll,1),2,1+2*(n-1))
%     h=max(speedAll(n,:));
% 
%     for i=1:size(cues,1)
%     hold on
%     patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 h h 0],'y')
%     end
% if ~isempty(intersect(n,[11:1:20]))
%      for i=1:size(stimLoc,1)
% hold on
%         patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 h h 0],'b')
%      end
% end
% hold on
% line([rewardLoc rewardLoc],[0 h],'Color','r')
%     hold on
%     % plot(xvalues,speedAll(n,:),'g');
%     semshade(speedAllSession{n},0.3,[0.5 0.5 0.5],xvalues);
%     % axis off
%     set(gca,'xticklabel',[])
%     % set(gca,'yticklabel',[])
% 
%     xlim([0 400])
% end
% 
% %
% % every 5
% nDay=[0 1 3 3 3 1 3 3 3 1 3 3 3];
% speedAllAvg=[];
% for n=1:size(nDay,2)-1;
%     idx=[sum(nDay(1:n),2)+1:1:sum(nDay(1:n+1),2)];
%     speedAllAvg(n,:)=nanmean(speedAll(idx,:),1);
% end
% 
% speedAllAvgSession={};
% for n=1:size(nDay,2)-1;
%     idx=[sum(nDay(1:n),2)+1:1:sum(nDay(1:n+1),2)];
%     speedAllAvgSession{n}=cell2mat(speedAllSession(idx));
% end
% 
% 
% 
% % plot
% 
% for n=1:size(speedAllAvg,1);
%     subplot(size(speedAll,1),2,[(2*(n-1)+1)*2:2:(2*(n-1)+2)*2])
%     h=max(speedAllAvg(n,:));
% 
%     for i=1:size(cues,1)
%     hold on
%     patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 h h 0],'y')
% end
% if ~isempty(intersect(n,[5:1:8]))
%      for i=1:size(stimLoc,1)
% hold on
%         patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 h h 0],'b')
%      end
% end
% hold on
% line([rewardLoc rewardLoc],[0 h],'Color','r')
%     hold on
%     i=find(speedAllAvg(n,:)>0);
%     % plot(xvalues,speedAllAvg(n,:),'g');
%    semshade(speedAllAvgSession{n},0.3,[0.5 0.5 0.5],xvalues);
%     xlim([0 400])
%         % axis off
%     set(gca,'xticklabel',[])
%     % set(gca,'yticklabel',[])
% end
% tightfig
% 
% saveas(gcf,'RBRDeceDistri_CHR2RS.fig')
% print -painters -depsc RBRDeceDistri_CHR2RS.eps

% %% no need to replot gfp cs
% idx=3;
% pre=allSPreRBRUseRun{idx};
% dur=allSDurRBRUseRun{idx};
% post=allSPostRBRUseRun{idx};
% 
% speedAll=[];%each row is the mean speed of all sessions
% for n=1:length(pre);
%     speedAll=[speedAll;nanmean(pre{n},1)];
% end
% for n=1:length(dur);
%     speedAll=[speedAll;nanmean(dur{n},1)];
% end
% for n=1:length(post);
%     speedAll=[speedAll;nanmean(post{n},1)];
% end
% 
% speedAllSession=[pre';dur';post'];
% 
% %
% % plot
% figure
% for n=1:size(speedAll,1);
%     subplot(size(speedAll,1),2,1+2*(n-1))
%     h=max(speedAll(n,:));
% 
%     for i=1:size(cues,1)
%     hold on
%     patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 h h 0],'y')
%     end
% if ~isempty(intersect(n,[11:1:20]))
%      for i=1:size(stimLoc,1)
% hold on
%         patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 h h 0],'b')
%      end
% end
% hold on
% line([rewardLoc rewardLoc],[0 h],'Color','r')
%     hold on
%     % plot(xvalues,speedAll(n,:),'g');
%     semshade(speedAllSession{n},0.3,'g',xvalues);
%         % axis off
%     set(gca,'xticklabel',[])
%     % set(gca,'yticklabel',[])
%     xlim([0 400])
% end
% 
% %
% % every 5
% nDay=[0 1 3 3 3 1 3 3 3 1 3 3 3];
% speedAllAvg=[];
% for n=1:size(nDay,2)-1;
%     idx=[sum(nDay(1:n),2)+1:1:sum(nDay(1:n+1),2)];
%     speedAllAvg(n,:)=nanmean(speedAll(idx,:),1);
% end
% 
% speedAllAvgSession={};
% for n=1:size(nDay,2)-1;
%     idx=[sum(nDay(1:n),2)+1:1:sum(nDay(1:n+1),2)];
%     speedAllAvgSession{n}=cell2mat(speedAllSession(idx));
% end
% 
% 
% 
% % plot
% 
% for n=1:size(speedAllAvg,1);
%     subplot(size(speedAll,1),2,[(2*(n-1)+1)*2:2:(2*(n-1)+2)*2])
%     h=max(speedAllAvg(n,:));
% 
%     for i=1:size(cues,1)
%     hold on
%     patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 h h 0],'y')
% end
% if ~isempty(intersect(n,[5:1:8]))
%      for i=1:size(stimLoc,1)
% hold on
%         patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 h h 0],'b')
%      end
% end
% hold on
% line([rewardLoc rewardLoc],[0 h],'Color','r')
%     hold on
%     i=find(speedAllAvg(n,:)>0);
%     % plot(xvalues,speedAllAvg(n,:),'g');
%    semshade(speedAllAvgSession{n},0.3,'g',xvalues);
%     xlim([0 400])
%         % axis off
%     set(gca,'xticklabel',[])
%     % set(gca,'yticklabel',[])
% end
% tightfig
% 
% saveas(gcf,'RBRDeceDistri_GFPCS.fig')
% print -painters -depsc RBRDeceDistri_GFPCS.eps
% %% no need to replot GFP RS
% idx=4;
% pre=allSPreRBRUseRun{idx};
% dur=allSDurRBRUseRun{idx};
% post=allSPostRBRUseRun{idx};
% 
% speedAll=[];%each row is the mean speed of all sessions
% for n=1:length(pre);
%     speedAll=[speedAll;nanmean(pre{n},1)];
% end
% for n=1:length(dur);
%     speedAll=[speedAll;nanmean(dur{n},1)];
% end
% for n=1:length(post);
%     speedAll=[speedAll;nanmean(post{n},1)];
% end
% 
% speedAllSession=[pre';dur';post'];
% 
% %
% % plot
% figure
% for n=1:size(speedAll,1);
%     subplot(size(speedAll,1),2,1+2*(n-1))
%     h=max(speedAll(n,:));
% 
%     for i=1:size(cues,1)
%     hold on
%     patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 h h 0],'y')
%     end
% if ~isempty(intersect(n,[11:1:20]))
%      for i=1:size(stimLoc,1)
% hold on
%         patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 h h 0],'b')
%      end
% end
% hold on
% line([rewardLoc rewardLoc],[0 h],'Color','r')
%     hold on
%     % plot(xvalues,speedAll(n,:),'g');
%     semshade(speedAllSession{n},0.3,[0.5 0.5 0.5],xvalues);
%         % axis off
%     set(gca,'xticklabel',[])
%     % set(gca,'yticklabel',[])
%     xlim([0 400])
% end
% 
% %
% % every 5
% nDay=[0 1 3 3 3 1 3 3 3 1 3 3 3];
% speedAllAvg=[];
% for n=1:size(nDay,2)-1;
%     idx=[sum(nDay(1:n),2)+1:1:sum(nDay(1:n+1),2)];
%     speedAllAvg(n,:)=nanmean(speedAll(idx,:),1);
% end
% 
% speedAllAvgSession={};
% for n=1:size(nDay,2)-1;
%     idx=[sum(nDay(1:n),2)+1:1:sum(nDay(1:n+1),2)];
%     speedAllAvgSession{n}=cell2mat(speedAllSession(idx));
% end
% 
% 
% 
% % plot
% 
% for n=1:size(speedAllAvg,1);
%     subplot(size(speedAll,1),2,[(2*(n-1)+1)*2:2:(2*(n-1)+2)*2])
%     h=max(speedAllAvg(n,:));
% 
%     for i=1:size(cues,1)
%     hold on
%     patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 h h 0],'y')
% end
% if ~isempty(intersect(n,[5:1:8]))
%      for i=1:size(stimLoc,1)
% hold on
%         patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 h h 0],'b')
%      end
% end
% hold on
% line([rewardLoc rewardLoc],[0 h],'Color','r')
%     hold on
%     i=find(speedAllAvg(n,:)>0);
%     % plot(xvalues,speedAllAvg(n,:),'g');
%    semshade(speedAllAvgSession{n},0.3,[0.5 0.5 0.5],xvalues);
%     xlim([0 400])
%         % axis off
%     set(gca,'xticklabel',[])
%     % set(gca,'yticklabel',[])
% end
% tightfig
% 
% saveas(gcf,'RBRDeceDistri_GFPRS.fig')
% print -painters -depsc RBRDeceDistri_GFPRS.eps
%%
load('allSpeedYUseLowerLick.mat')
idx=1; %chr2cs
speedBefore=allSPreUse{idx};
speedStim=allSDurUse{idx};
speedAfter=allSPostUse{idx};
yBefore=allYPreUse{idx};
yStim=allYDurUse{idx};
yAfter=allYPostUse{idx};
[decePercentileBefore,decePercentileStim,decePercentileAfter] = decelerationPercCueStim(speedBefore,speedStim,speedAfter,yBefore,yStim,yAfter,stimLoc,stimExp1,stimExp2,samplingSpacing);
save('decePercentNearStimCHR2CSLowerLick.mat','decePercentileBefore','decePercentileStim','decePercentileAfter')

% idx=2; %chr2rs
% speedBefore=allSPreUse{idx};
% speedStim=allSDurUse{idx};
% speedAfter=allSPostUse{idx};
% yBefore=allYPreUse{idx};
% yStim=allYDurUse{idx};
% yAfter=allYPostUse{idx};
% [decePercentileBefore,decePercentileStim,decePercentileAfter] = decelerationPercCueStim(speedBefore,speedStim,speedAfter,yBefore,yStim,yAfter,stimLoc,stimExp1,stimExp2,samplingSpacing);
% save('decePercentNearStimCHR2RSLowerLick.mat','decePercentileBefore','decePercentileStim','decePercentileAfter')

% idx=3; %GFPcs
% speedBefore=allSPreUse{idx};
% speedStim=allSDurUse{idx};
% speedAfter=allSPostUse{idx};
% yBefore=allYPreUse{idx};
% yStim=allYDurUse{idx};
% yAfter=allYPostUse{idx};
% [decePercentileBefore,decePercentileStim,decePercentileAfter] = decelerationPercCueStim(speedBefore,speedStim,speedAfter,yBefore,yStim,yAfter,stimLoc,stimExp1,stimExp2,samplingSpacing);
% save('decePercentNearStimGFPCS.mat','decePercentileBefore','decePercentileStim','decePercentileAfter')
% 
% idx=4; %GFPcs
% speedBefore=allSPreUse{idx};
% speedStim=allSDurUse{idx};
% speedAfter=allSPostUse{idx};
% yBefore=allYPreUse{idx};
% yStim=allYDurUse{idx};
% yAfter=allYPostUse{idx};
% [decePercentileBefore,decePercentileStim,decePercentileAfter] = decelerationPercCueStim(speedBefore,speedStim,speedAfter,yBefore,yStim,yAfter,stimLoc,stimExp1,stimExp2,samplingSpacing);
% save('decePercentNearStimGFPRS.mat','decePercentileBefore','decePercentileStim','decePercentileAfter')

%%

sessionDataAll={};
load('decePercentNearStimCHR2CSLowerLick.mat')
sessionDataAll{1}=[decePercentileBefore' decePercentileStim' decePercentileAfter'];
load('decePercentNearStimCHR2RS.mat')
sessionDataAll{2}=[decePercentileBefore' decePercentileStim' decePercentileAfter'];
load('decePercentNearStimGFPCS.mat')
sessionDataAll{3}=[decePercentileBefore' decePercentileStim' decePercentileAfter'];
load('decePercentNearStimGFPRS.mat')
sessionDataAll{4}=[decePercentileBefore' decePercentileStim' decePercentileAfter'];

figure
for t=1:length(sessionDataAll)
allRBR=sessionDataAll{t};
% if t==2; %for chr2RS
% allRBR=allRBR(1:end~=9,:); %remove this mouse: has extremely low baseline and made the following RBR analysis weired
% end

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
    plot([1:1:3],allSession(i,:),'g-')
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

saveas(gcf,'decePercentNearStimPerSessionLowerLick.fig')

%% RBR


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
for n=1:size(A,1);
   ma=nanmean(A(n,[1:10]),2);
    A(n,:)=A(n,:)/ma;
end
A(isinf(A))=nan;
A=A(1:end~=8,:);%this mouse was eliminated because the basedline was very low compared to all other mice
% ma=nanmean(A(:,[1:10]),'all');
% A=A/ma;

B=RBRCS(:,1:N);
for n=1:size(B,1);
  ma=nanmean(B(n,[1:10]),2);
    B(n,:)=B(n,:)/ma;
end
B(isinf(B))=nan;
% ma=nanmean(B(:,[1:10]),'all');
% B=B/ma;

hold on
errorbar([1:1:30],nanmean(A,1),nansem(A,1),'k')
hold on
errorbar([1:1:30],nanmean(B,1),nansem(B,1),'g')
% semshade(A,0.3,'k',[1:N])
% hold on
% semshade(B,0.3,'m',[1:N])
%

 [pAnova,pMC,pLabel] = anovaRM2W_full_BH(A(:,11:20),B(:,11:20),1);
 [pAnova2,pMC,pLabel] = anovaRM2W_full_BH(A(:,1:10),B(:,1:10),1);
 [pAnova3,pMC,pLabel] = anovaRM2W_full_BH(A(:,21:30),B(:,21:30),1);


% ylim([0 4])
title(['chr2 CS RS,P=',num2str([pAnova(1) pAnova2(1) pAnova3(1)])])
% title(['chr2 CS RS'])

xlim([0 N])
p1=[];
for n=1:N;
    [r,p1(n)]= ttest2(A(:,n),B(:,n));
end


for n=1:N;
    if p1(n)<=0.05;
        hold on
        plot(n,1.2,'k*')
    end
end

subplot(222)
A=RBRRSGFP(:,1:N);
for n=1:size(A,1);
  ma=nanmean(A(n,[1:10]),2);
    A(n,:)=A(n,:)/ma;
end
A(isinf(A))=nan;
% ma=nanmean(A(:,[1:10]),'all');
% A=A/ma;

B=RBRCSGFP(:,1:N);
for n=1:size(B,1);
  ma=nanmean(B(n,[1:10]),2);
    B(n,:)=B(n,:)/ma;
end
B(isinf(B))=nan;
% ma=nanmean(B(:,[1:10]),'all');
% B=B/ma;

hold on
errorbar([1:1:30],nanmean(A,1),nansem(A,1),'k')
hold on
errorbar([1:1:30],nanmean(B,1),nansem(B,1),'g')
 [pAnova,pMC,pLabel] = anovaRM2W_full_BH(A(:,11:20),B(:,11:20),1);
 [pAnova2,pMC,pLabel] = anovaRM2W_full_BH(A(:,1:10),B(:,1:10),1);
 [pAnova3,pMC,pLabel] = anovaRM2W_full_BH(A(:,21:30),B(:,21:30),1);
% ylim([0 4])

% ylim([0.6 1.3])
title(['gfp CS RS, P=',num2str([pAnova(1) pAnova2(1) pAnova3(1)])])
% title(['gfp CS RS']);

xlim([0 N])
p2=[];
for n=1:N;
    [r,p2(n)]= ttest2(A(:,n),B(:,n));
end
for n=1:N;
    if p2(n)<=0.05;
        hold on
        plot(n,1.2,'k*')
    end
end

%

subplot(224)
A=RBRCS(:,1:N);
for n=1:size(A,1);
   ma=nanmean(A(n,[1:10]),2);
    A(n,:)=A(n,:)/ma;
end
A(isinf(A))=nan;
% ma=nanmean(A(:,[1:10]),'all');
% A=A/ma;
B=RBRCSGFP(:,1:N);
for n=1:size(B,1);
  ma=nanmean(B(n,[1:10]),2);
    B(n,:)=B(n,:)/ma;
end
B(isinf(B))=nan;
% ma=nanmean(B(:,[1:10]),'all');
% B=B/ma;
hold on
errorbar([1:1:30],nanmean(A,1),nansem(A,1),'k')
hold on
errorbar([1:1:30],nanmean(B,1),nansem(B,1),'g')
 [pAnova,pMC,pLabel] = anovaRM2W_full_BH(A(:,11:20),B(:,11:20),1);
 [pAnova2,pMC,pLabel] = anovaRM2W_full_BH(A(:,1:10),B(:,1:10),1);
 [pAnova3,pMC,pLabel] = anovaRM2W_full_BH(A(:,21:30),B(:,21:30),1);


% ylim([0.6 1.3])
title(['CS CHR2 GFP, P=',num2str([pAnova(1) pAnova2(1) pAnova3(1)])])
% title(['CS CHR2 GFP'])

xlim([0 N])
p3=[];
for n=1:N;
    [r,p3(n)]= ttest2(A(:,n),B(:,n));
end
for n=1:N;
    if p3(n)<=0.05;
        hold on
        plot(n,1.2,'k*')
    end
end



subplot(223)
A=RBR(:,1:N);
for n=1:size(A,1);
    ma=nanmean(A(n,[1:10]),2);
    A(n,:)=A(n,:)/ma;
end
A(isinf(A))=nan;
% ma=nanmean(A(:,[1:10]),'all');
% A=A/ma;
B=RBRRSGFP(:,1:N);
for n=1:size(B,1);
    ma=nanmean(B(n,[1:10]),2);
    B(n,:)=B(n,:)/ma;
end
B(isinf(B))=nan;
% ma=nanmean(B(:,[1:10]),'all');
% B=B/ma;
hold on
errorbar([1:1:30],nanmean(A,1),nansem(A,1),'k')
hold on
errorbar([1:1:30],nanmean(B,1),nansem(B,1),'g')
 [pAnova,pMC,pLabel] = anovaRM2W_full_BH(A(:,11:20),B(:,11:20),1);
  [pAnova2,pMC,pLabel] = anovaRM2W_full_BH(A(:,1:10),B(:,1:10),1);
 [pAnova3,pMC,pLabel] = anovaRM2W_full_BH(A(:,21:30),B(:,21:30),1);


% ylim([0.6 1.3])
title(['RS CHR2 GFP, P=',num2str([pAnova(1) pAnova2(1) pAnova3(1)])])
% title(['RS CHR2 GFP'])
xlim([0 N])
p4=[];
for n=1:N;
    [r,p4(n)]= ttest2(A(:,n),B(:,n));
end
for n=1:N;
    if p4(n)<=0.05;
        hold on
        plot(n,1.2,'k*')
    end
end
saveas(gcf,'decePercentNearStimRBRLowerLick.fig')

%% RBR: normalized to zero


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
for n=1:size(A,1);
   ma=nanmean(A(n,[1:10]),2);
    A(n,:)=A(n,:)-ma;
end
A(isinf(A))=nan;
A=A(1:end~=8,:);%this mouse was eliminated because the basedline was very low compared to all other mice
% ma=nanmean(A(:,[1:10]),'all');
% A=A-ma;

B=RBRCS(:,1:N);
for n=1:size(B,1);
  ma=nanmean(B(n,[1:10]),2);
    B(n,:)=B(n,:)-ma;
end
B(isinf(B))=nan;
% ma=nanmean(B(:,[1:10]),'all');
% B=B-ma;

hold on
errorbar([1:1:30],nanmean(A,1),nansem(A,1),'k')
hold on
errorbar([1:1:30],nanmean(B,1),nansem(B,1),'g')
% semshade(A,0.3,'k',[1:N])
% hold on
% semshade(B,0.3,'m',[1:N])
%

 [pAnova,pMC,pLabel] = anovaRM2W_full_BH(A(:,11:20),B(:,11:20),1);
 [pAnova2,pMC,pLabel] = anovaRM2W_full_BH(A(:,1:10),B(:,1:10),1);
 [pAnova3,pMC,pLabel] = anovaRM2W_full_BH(A(:,21:30),B(:,21:30),1);

A2=reshape(A(:,[11:20]),[1 size(A,1)*10]);
A1=reshape(A(:,[1:10]),[1 size(A,1)*10]);

B2=reshape(B(:,[11:20]),[1 size(B,1)*10]);
B1=reshape(B(:,[1:10]),[1 size(B,1)*10]);

[~,pPreDurRS]=ttest2(A2,A1)
[~,pPreDurCS]=ttest2(B2,B1)
% ylim([0 4])
title(['chr2 CS RS,P=',num2str([pAnova(1) pAnova2(1) pAnova3(1) pPreDurCS pPreDurRS])])
% title(['chr2 CS RS'])

xlim([0 N])
p1=[];
for n=1:N;
    [r,p1(n)]= ttest2(A(:,n),B(:,n));
end


for n=1:N;
    if p1(n)<=0.05;
        hold on
        plot(n,0.4,'k*')
    end
end

subplot(222)
A=RBRRSGFP(:,1:N);
for n=1:size(A,1);
  ma=nanmean(A(n,[1:10]),2);
    A(n,:)=A(n,:)-ma;
end
A(isinf(A))=nan;
% ma=nanmean(A(:,[1:10]),'all');
% A=A-ma;

B=RBRCSGFP(:,1:N);
for n=1:size(B,1);
  ma=nanmean(B(n,[1:10]),2);
    B(n,:)=B(n,:)-ma;
end
B(isinf(B))=nan;
% ma=nanmean(B(:,[1:10]),'all');
% B=B-ma;

hold on
errorbar([1:1:30],nanmean(A,1),nansem(A,1),'k')
hold on
errorbar([1:1:30],nanmean(B,1),nansem(B,1),'g')
 [pAnova,pMC,pLabel] = anovaRM2W_full_BH(A(:,11:20),B(:,11:20),1);
 [pAnova2,pMC,pLabel] = anovaRM2W_full_BH(A(:,1:10),B(:,1:10),1);
 [pAnova3,pMC,pLabel] = anovaRM2W_full_BH(A(:,21:30),B(:,21:30),1);

 A2=reshape(A(:,[11:20]),[1 size(A,1)*10]);
A1=reshape(A(:,[1:10]),[1 size(A,1)*10]);

B2=reshape(B(:,[11:20]),[1 size(B,1)*10]);
B1=reshape(B(:,[1:10]),[1 size(B,1)*10]);

[~,pPreDurRS]=ttest2(A2,A1)
[~,pPreDurCS]=ttest2(B2,B1)
% ylim([0 4])

% ylim([0.6 1.3])
title(['gfp CS RS, P=',num2str([pAnova(1) pAnova2(1) pAnova3(1) pPreDurCS pPreDurRS])])
% title(['gfp CS RS']);

xlim([0 N])
p2=[];
for n=1:N;
    [r,p2(n)]= ttest2(A(:,n),B(:,n));
end
for n=1:N;
    if p2(n)<=0.05;
        hold on
        plot(n,0.4,'k*')
    end
end

%

subplot(224)
A=RBRCS(:,1:N);
for n=1:size(A,1);
   ma=nanmean(A(n,[1:10]),2);
    A(n,:)=A(n,:)-ma;
end
A(isinf(A))=nan;
% ma=nanmean(A(:,[1:10]),'all');
% A=A-ma;
B=RBRCSGFP(:,1:N);
for n=1:size(B,1);
  ma=nanmean(B(n,[1:10]),2);
    B(n,:)=B(n,:)-ma;
end
B(isinf(B))=nan;
% ma=nanmean(B(:,[1:10]),'all');
% B=B-ma;
hold on
errorbar([1:1:30],nanmean(A,1),nansem(A,1),'k')
hold on
errorbar([1:1:30],nanmean(B,1),nansem(B,1),'g')
 [pAnova,pMC,pLabel] = anovaRM2W_full_BH(A(:,11:20),B(:,11:20),1);
 [pAnova2,pMC,pLabel] = anovaRM2W_full_BH(A(:,1:10),B(:,1:10),1);
 [pAnova3,pMC,pLabel] = anovaRM2W_full_BH(A(:,21:30),B(:,21:30),1);


% ylim([0.6 1.3])
title(['CS CHR2 GFP, P=',num2str([pAnova(1) pAnova2(1) pAnova3(1)])])
% title(['CS CHR2 GFP'])

xlim([0 N])
p3=[];
for n=1:N;
    [r,p3(n)]= ttest2(A(:,n),B(:,n));
end
for n=1:N;
    if p3(n)<=0.05;
        hold on
        plot(n,0.4,'k*')
    end
end



subplot(223)
A=RBR(:,1:N);
for n=1:size(A,1);
    ma=nanmean(A(n,[1:10]),2);
    A(n,:)=A(n,:)-ma;
end
A(isinf(A))=nan;
% ma=nanmean(A(:,[1:10]),'all');
% A=A-ma;
B=RBRRSGFP(:,1:N);
for n=1:size(B,1);
    ma=nanmean(B(n,[1:10]),2);
    B(n,:)=B(n,:)-ma;
end
B(isinf(B))=nan;
% ma=nanmean(B(:,[1:10]),'all');
% B=B-ma;
hold on
errorbar([1:1:30],nanmean(A,1),nansem(A,1),'k')
hold on
errorbar([1:1:30],nanmean(B,1),nansem(B,1),'g')
 [pAnova,pMC,pLabel] = anovaRM2W_full_BH(A(:,11:20),B(:,11:20),1);
  [pAnova2,pMC,pLabel] = anovaRM2W_full_BH(A(:,1:10),B(:,1:10),1);
 [pAnova3,pMC,pLabel] = anovaRM2W_full_BH(A(:,21:30),B(:,21:30),1);


% ylim([0.6 1.3])
title(['RS CHR2 GFP, P=',num2str([pAnova(1) pAnova2(1) pAnova3(1)])])
% title(['RS CHR2 GFP'])
xlim([0 N])
p4=[];
for n=1:N;
    [r,p4(n)]= ttest2(A(:,n),B(:,n));
end
for n=1:N;
    if p4(n)<=0.05;
        hold on
        plot(n,0.4,'k*')
    end
end
saveas(gcf,'decePercentNearStimRBRZeroLowerLick.fig')
%% 

track=zeros(400,1);
stimLoc=[];
stimLoc(1,:)=[63+3.5 63+3.5+15];
stimLoc(2,:)=[163+3.5 163+3.5+15];
stimLoc(3,:)=[279+3.5 279+3.5+15];
stimExp1=3;
stimExp2=3;
stimLoc(:,1)=stimLoc(:,1)-stimExp1;
stimLoc(:,2)=stimLoc(:,2)+stimExp2;

for n=1:size(stimLoc,1)
    track(ceil(stimLoc(n,1)):ceil(stimLoc(n,2)))=1;
end
useRegion=[16:400];%remove the first 15: many mice nade stop there
useRegion=[21:380];
load('allSRBRUseRunLowerLick.mat')
window=10;
deceDur={};
corrDur={};

for type=1:length(allSDurRBRUseRun);
    deceDur{type}={};
    for lap=1:length(allSDurRBRUseRun{type});
        deceDur{type}{lap}=[];
        for session=1:size(allSDurRBRUseRun{type}{lap},1);
            speed=allSDurRBRUseRun{type}{lap}(session,:);
            speed=smoothdata(speed,'movmean',window);
            dece=-[0 diff(speed)];
            deceDur{type}{lap}(session,:)=dece;
            corrDur{type}(lap,session)=corr(dece(useRegion)',track(useRegion));
        end
    end
end
%

decePost={};
corrPost={};

for type=1:length(allSPostRBRUseRun);
    decePost{type}={};
    for lap=1:length(allSPostRBRUseRun{type});
        decePost{type}{lap}=[];
        for session=1:size(allSPostRBRUseRun{type}{lap},1);
            speed=allSPostRBRUseRun{type}{lap}(session,:);
            speed=smoothdata(speed,'movmean',window);
            dece=-[0 diff(speed)];
            decePost{type}{lap}(session,:)=dece;
            corrPost{type}(lap,session)=corr(dece(useRegion)',track(useRegion));
        end
    end
end
%
decePre={};
corrPre={};

for type=1:length(allSPreRBRUseRun);
    decePre{type}={};
    for lap=1:length(allSPreRBRUseRun{type});
        decePre{type}{lap}=[];
        for session=1:size(allSPreRBRUseRun{type}{lap},1);
            speed=allSPreRBRUseRun{type}{lap}(session,:);
            speed=smoothdata(speed,'movmean',window);
            dece=-[0 diff(speed)];
            decePre{type}{lap}(session,:)=dece;
            corrPre{type}(lap,session)=corr(dece(useRegion)',track(useRegion));
        end
    end
end
%
corrAll={};
for n=1:4;
    corrAll{n}=[corrPre{n};corrDur{n};corrPost{n}];
end
%
figure
for n=1:length(corrAll)
    subplot(3,2,n);
     if n==1 | n==3;

    errorbar([1:1:size(corrAll{1},1)],nanmean(corrAll{n},2),nansem(corrAll{n},2),'g');
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
    ylim([-0.24 0.4])
   
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
plot(useRegion,track(useRegion))
hold on
plot(useRegion,deceDur{1}{1}(6,useRegion),'g')


    saveas(gcf,'corrDeceWithStimPatternLowerLick.fig')
print -painters -depsc corrDeceWithStimPatternLowerLick.eps
%% correlation between stim temp and dece: including RS pattern
load('allSRBRUseRunLowerLick.mat')
window=10;
% window=40;
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

track=zeros(400,1);
useRegion=[16:400];
useRegion=[21:380];
positionAll={};
for N=1:length(allIdxAll);
    positionAll{N}={};
mouseIdx=allIdxAll{N};

position=positionSegAllData(mouseIdx);
stim=stimZoneIdxAllDataFix(mouseIdx);
expRight=3;
expLeft=3;
for n=1:length(position); %per session
    for m=1:length(position{n}) %per run
        positionAll{N}{n,m}=position{n}{m}(stim{n}{m});
        positionAll{N}{n,m}(:,2)=positionAll{N}{n,m}(:,2)+expRight;
         positionAll{N}{n,m}(:,1)=positionAll{N}{n,m}(:,1)-expLeft;
         a=positionAll{N}{n,m}(:,1);
         a(a<1)=1;
         positionAll{N}{n,m}(:,1)=a;
    end
end
end
deceDur={};
corrDur={};
trackDur={};
for N=1:length(positionAll);
    corrDur{N}=[];
    deceDur{N}={};
    trackDur{N}={};
    for session=1:size(positionAll{N},1);
        for lap=1:size(positionAll{N},2);
            positionThis=positionAll{N}{session,lap};
            trackThis=track;
            for i=1:size(positionThis,1);
                trackThis(ceil(positionThis(i,1)):ceil(positionThis(i,2)))=1;
            end
            trackThis=trackThis(useRegion);
            speed=allSDurRBRUseRun{N}{lap}(session,:);
            speed=smoothdata(speed,'movmean',window);
            dece=-[0 diff(speed)];
            dece=dece(useRegion);
            deceDur{N}{session,lap}=dece;
            trackDur{N}{session,lap}=trackThis;
            corrDur{N}(session,lap)=corr(dece',trackThis);
        end
    end
end

%

figure
for N=1:length(corrDur)
    subplot(5,2,N)
    C=corrDur{N};
    errorbar([1:1:size(C,2)],nanmean(C,1),nansem(C,1),'k');

     [pAnova,pMC,pLabel] = anovaRM2W_full_BH(C,zeros(size(C,1),size(C,2)),1);
    ylim([-0.1 0.35])
    hold on
    line([0 size(C,2)],[0 0])
    pvalue=[];
    for n=1:size(C,2);
        [~,pvalue(n)]=ttest(C(:,n),0);
        if pvalue(n)<=0.05;
            plot(n,0.2,'r*')
        end
    end

    [~,p]=corr([1:1:10]',nanmean(C,1)');
    title([num2str([pAnova(1) p])])


end

subplot(5,2,5)
N=1;
 C1=corrDur{N};
    errorbar([1:1:size(C1,2)],nanmean(C1,1),nansem(C1,1),'g');
    hold on
    N=2;
   C2=corrDur{N};
   hold on
    errorbar([1:1:size(C2,2)],nanmean(C2,1),nansem(C2,1),'Color',[0.5 0.5 0.5]);  
      ylim([-0.1 0.35])
    hold on
    line([0 size(C1,2)],[0 0])

      pvalue=[];
    for n=1:size(C1,2);
        [~,pvalue(n)]=ttest2(C1(:,n),C2(:,n));
        if pvalue(n)<=0.05;
            plot(n,0.2,'r*')
        end
    end

     [pAnova,pMC,pLabel] = anovaRM2W_full_BH(C1(:,5:end),C2(:,5:end),1);
    title(['CHR2 CS vs RS, ANOVA',num2str(pAnova(1))])

    subplot(5,2,7)

    C1=C1(:,[5:end]);
    C1=reshape(C1,[size(C1,1)*size(C1,2) 1]);
    C2=C2(:,[5:end]);
    C2=reshape(C2,[size(C2,1)*size(C2,2) 1]);

    M=[];
    E=[];
    M(1)=nanmean(C1,1);
    M(2)=nanmean(C2,1);

    E(1)=nansem(C1,1);
    E(2)=nansem(C2,1);

    bar([1 2],M)
    hold on
    errorbar([1 2],M,E)

    [r,p]=ttest2(C1,C2);
    title([num2str(p)])

    subplot(5,2,6)
N=3;
 C1=corrDur{N};
    errorbar([1:1:size(C1,2)],nanmean(C1,1),nansem(C1,1),'g');
    hold on
    N=4;
   C2=corrDur{N};
   hold on
    errorbar([1:1:size(C2,2)],nanmean(C2,1),nansem(C2,1),'Color',[0.5 0.5 0.5]);  
      ylim([-0.15 0.35])
    hold on
    line([0 size(C1,2)],[0 0])

      pvalue=[];
    for n=1:size(C1,2);
        [~,pvalue(n)]=ttest2(C1(:,n),C2(:,n));
        if pvalue(n)<=0.05;
            plot(n,0.2,'r*')
        end
    end
         [pAnova,pMC,pLabel] = anovaRM2W_full_BH(C1(:,2:end),C2(:,2:end),1);

    title(['GFP CS vs RS, ANOVA',num2str(pAnova(1))])

     subplot(5,2,8)

    C1=C1(:,[2:end]);
    C1=reshape(C1,[size(C1,1)*size(C1,2) 1]);
    C2=C2(:,[2:end]);
    C2=reshape(C2,[size(C2,1)*size(C2,2) 1]);

    M=[];
    E=[];
    M(1)=nanmean(C1,1);
    M(2)=nanmean(C2,1);

    E(1)=nansem(C1,1);
    E(2)=nansem(C2,1);

    bar([1 2],M)
    hold on
    errorbar([1 2],M,E)

    [r,p]=ttest2(C1,C2);
    title([num2str(p)])

%plot example
type=2;
session=3;
lap=2;
subplot(5,2,9)
hold on
plot(useRegion,trackDur{N}{session,lap},'k');
hold on
plot(useRegion,deceDur{N}{session,lap},'g');

type=2;
session=3;
lap=7;
subplot(5,2,10)
hold on
plot(useRegion,trackDur{N}{session,lap},'k');
hold on
plot(useRegion,deceDur{N}{session,lap},'g');

    saveas(gcf,'corrDeceWithStimPatternForRSLowerLick.fig')
    print -painters -depsc corrDeceWithStimPatternForRSLowerLick.eps


           
%% aligning stim zone: more for RS
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


startBinStim=[];%each row is a session/run: %run 1-10 for session1, then run 1-10 for session2... the three numbers are the start point of the three stim
for session=1:size(positionAll,1)
    for run=1:size(positionAll,2);
        startPoint=positionAll{session,run}(:,1);
        startBinStim(end+1,:)=ceil(startPoint);
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

beforeAfterZone=15; %first 15min: many sessions show a large deceleration
stimStart=15+beforeAfterZone;
stimEnd=310-beforeAfterZone-15+1;



speed=[];
for n=1:size(startBinStim,1)
    for m=1:size(startBinStim,2)
        startThis=startBinStim(n,m);
        if startThis-stimStart>0 & startThis-stimEnd<0;
            speed(end+1,:)=RBRSpeed(n,startThis-beforeAfterZone:startThis+15+beforeAfterZone-1);
        end
    end
end

subplot(4,2,2*(N-1)+1)
semshade(speed,0.3,'g',[1:1:size(speed,2)]);
M=nanmean(speed,1);
E=nansem(speed,1);
y2=max(M)+max(E);
y1=min(M)-min(E);
hold on
line([beforeAfterZone beforeAfterZone],[y1 y2])
hold on
line([beforeAfterZone+15 beforeAfterZone+15],[y1 y2])
ylim([24 43])

zones={};
zones{1}=[1:1:beforeAfterZone];
zones{2}=[beforeAfterZone+1:1:beforeAfterZone+15];
zones{3}=[beforeAfterZone+16:1:size(speed,2)];

meanSpeed=[];
meanSpeed=[nanmean(speed(:,zones{1}),2) nanmean(speed(:,zones{2}),2) nanmean(speed(:,zones{3}),2)];

meanSpeedAll{N}=meanSpeed;
subplot(4,2,N*2);
bar([1:1:3],nanmean(meanSpeed,1),'FaceColor','w');
hold on
errorbar([1:1:3],nanmean(meanSpeed,1),nansem(meanSpeed,1),'k')
ylim([0 45])
[~,p1]=ttest(meanSpeed(:,1),meanSpeed(:,2));
[~,p2]=ttest(meanSpeed(:,2),meanSpeed(:,3));
[~,p3]=ttest(meanSpeed(:,1),meanSpeed(:,3));
title(num2str([p1 p2 p3]))
end

saveas(gcf,'aligniningStimZoneLowerLick.fig')
print -painters -depsc aligniningStimZoneLowerLick.eps

%plot speed change %NO DIFFERENCE
speedChangeCS=meanSpeedAll{1}(:,1)-meanSpeedAll{1}(:,2);
speedChangeRS=meanSpeedAll{2}(:,1)-meanSpeedAll{2}(:,2);
[R,P]=ttest2(speedChangeRS,speedChangeCS)

%% aligning stim zone: more for RS: only take the stim region with at least 100 cm from other loc (the beginning region can have less)
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

beforeAfterZone=15; %first 15min: many sessions show a large deceleration
stimStart=15+beforeAfterZone;
stimEnd=310-beforeAfterZone-15+1;

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
zones{1}=[1:1:beforeAfterZone];
zones{2}=[beforeAfterZone+1:1:beforeAfterZone+15];
zones{3}=[beforeAfterZone+16:1:size(speed,2)];


meanSpeed=[];
meanSpeed=[nanmean(speed(:,zones{1}),2) nanmean(speed(:,zones{2}),2) nanmean(speed(:,zones{3}),2)];

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
xlim([1 45])
if N==1 | N==2;
ylim([25 37])
else
    ylim([28 43]);
end


meanSpeedAll{N}=meanSpeed;
subplot(4,2,N*2);
bar([1:1:3],nanmean(meanSpeed,1),'FaceColor','w');
hold on
errorbar([1:1:3],nanmean(meanSpeed,1),nansem(meanSpeed,1),'k')
ylim([0 45])
[~,p1]=ttest(meanSpeed(:,1),meanSpeed(:,2));
[~,p2]=ttest(meanSpeed(:,2),meanSpeed(:,3));
[~,p3]=ttest(meanSpeed(:,1),meanSpeed(:,3));
title(num2str([p1 p2 p3]))
end

saveas(gcf,'aligniningStimZoneLowerLickFarApart.fig')
print -painters -depsc aligniningStimZoneLowerLickFarApart.eps

%plot speed change %NO DIFFERENCE
speedChangeCS=meanSpeedAll{1}(:,1)-meanSpeedAll{1}(:,2);
speedChangeRS=meanSpeedAll{2}(:,1)-meanSpeedAll{2}(:,2);
[R,P]=ttest2(speedChangeRS,speedChangeCS)


%% 50cm before and after: aligning stim zone: more for RS: only take the stim region with at least 100 cm from other loc (the beginning region can have less)
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

beforeAfterZone=50; %first 15min: many sessions show a large deceleration
stimStart=15+beforeAfterZone;
stimEnd=310-beforeAfterZone-15+1;

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
zones{1}=[1:1:beforeAfterZone];
zones{2}=[beforeAfterZone+1:1:beforeAfterZone+15];
zones{3}=[beforeAfterZone+16:1:size(speed,2)];


meanSpeed=[];
meanSpeed=[nanmean(speed(:,zones{1}),2) nanmean(speed(:,zones{2}),2) nanmean(speed(:,zones{3}),2)];

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

if N==1 | N==2;
ylim([25 42])
else
    ylim([28 43]);
end


meanSpeedAll{N}=meanSpeed;
subplot(4,2,N*2);
bar([1:1:3],nanmean(meanSpeed,1),'FaceColor','w');
hold on
errorbar([1:1:3],nanmean(meanSpeed,1),nansem(meanSpeed,1),'k')
ylim([0 45])
[~,p1]=ttest(meanSpeed(:,1),meanSpeed(:,2));
[~,p2]=ttest(meanSpeed(:,2),meanSpeed(:,3));
[~,p3]=ttest(meanSpeed(:,1),meanSpeed(:,3));
title(num2str([p1 p2 p3]))
end

saveas(gcf,'aligniningStimZoneLowerLickFarApart50cm.fig')
print -painters -depsc aligniningStimZoneLowerLickFarApart50cm.eps

%% 50 cm before after and expanded stim zone, aligning stim zone: more for RS: only take the stim region with at least 100 cm from other loc (the beginning region can have less)
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
exp=10;
beforeAfterZone=50; %here we use before after region and the before regiion we used 50cm, but after region, since there is a 10 cm expansion, we only included 40cm after that (so the total region after stim is still 50cm)
stimStart=15+beforeAfterZone;
stimEnd=310-beforeAfterZone-15-exp+1;

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
zones{1}=[1:1:beforeAfterZone];
zones{2}=[beforeAfterZone+1:1:beforeAfterZone+15];
zones{3}=[beforeAfterZone+16:1:beforeAfterZone+15+exp];
zones{4}=[beforeAfterZone+16+exp:1:size(speed,2)];
meanSpeed=[];
meanSpeed=[nanmean(speed(:,zones{1}),2) nanmean(speed(:,zones{2}),2) nanmean(speed(:,zones{3}),2) nanmean(speed(:,zones{4}),2)];
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

if N==1 | N==2;
ylim([25 42])
else
    ylim([28 44]);
end



subplot(4,2,N*2);
bar([1:1:4],nanmean(meanSpeed,1),'FaceColor','w');
hold on
errorbar([1:1:4],nanmean(meanSpeed,1),nansem(meanSpeed,1),'k')
ylim([0 45])
[~,p1]=ttest(meanSpeed(:,1),meanSpeed(:,2));
[~,p2]=ttest(meanSpeed(:,1),meanSpeed(:,3));
[~,p3]=ttest(meanSpeed(:,2),meanSpeed(:,3));
[~,p4]=ttest(meanSpeed(:,1),meanSpeed(:,4));
title(num2str([p1 p2 p3 p4]))
end

saveas(gcf,'aligniningStimZoneLowerLickFarApart50cmExp.fig')
print -painters -depsc aligniningStimZoneLowerLickFarApart50cmExp.eps

%% Downsample all to CHR2RS: 50 cm before after and expanded stim zone, aligning stim zone: more for RS: only take the stim region with at least 100 cm from other loc (the beginning region can have less)
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
exp=10;
beforeAfterZone=50; %first 15min: many sessions show a large deceleration
stimStart=15+beforeAfterZone;
stimEnd=310-beforeAfterZone-15-exp+1;

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
zones{1}=[1:1:beforeAfterZone];
zones{2}=[beforeAfterZone+1:1:beforeAfterZone+15];
zones{3}=[beforeAfterZone+16:1:beforeAfterZone+15+exp];
zones{4}=[beforeAfterZone+16+exp:1:size(speed,2)];
meanSpeed=[];
meanSpeed=[nanmean(speed(:,zones{1}),2) nanmean(speed(:,zones{2}),2) nanmean(speed(:,zones{3}),2) nanmean(speed(:,zones{4}),2)];
NPoints(N,1)=size(meanSpeed,1);
meanSpeedAll{N}=meanSpeed;

i=randperm(size(speed,1));
i=i(1:35);
speed=speed(i,:);
meanSpeed=meanSpeed(i,:);
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

if N==1 | N==2;
ylim([25 42])
else
    ylim([28 44]);
end



subplot(4,2,N*2);
bar([1:1:4],nanmean(meanSpeed,1),'FaceColor','w');
hold on
errorbar([1:1:4],nanmean(meanSpeed,1),nansem(meanSpeed,1),'k')
ylim([0 45])
[~,p1]=ttest(meanSpeed(:,1),meanSpeed(:,2));
[~,p2]=ttest(meanSpeed(:,1),meanSpeed(:,3));
[~,p3]=ttest(meanSpeed(:,2),meanSpeed(:,3));
[~,p4]=ttest(meanSpeed(:,1),meanSpeed(:,4));
title(num2str([p1 p2 p3 p4]))
end

saveas(gcf,'aligniningStimZoneLowerLickFarApart50cmExpDownsample.fig')
print -painters -depsc aligniningStimZoneLowerLickFarApart50cmExpDownsample.eps

%%  whether this is because the cs data had more points

NShuffle=1000;
NTake=size(meanSpeedAll{2},1);

allDownsamplePValues={};

idx=1;%downsample CHR2
NTotal=size(meanSpeedAll{idx},1);
pValuesShuffle=[];
for n=1:NShuffle;
    i=randperm(NTotal);
    i=i(1:NTake);
    meanSpeedUse=meanSpeedAll{idx}(i,:);
    [~,pValuesShuffle(n,1)]=ttest(meanSpeedUse(:,1),meanSpeedUse(:,2));
[~,pValuesShuffle(n,2)]=ttest(meanSpeedUse(:,1),meanSpeedUse(:,3));
[~,pValuesShuffle(n,3)]=ttest(meanSpeedUse(:,2),meanSpeedUse(:,3));
[~,pValuesShuffle(n,4)]=ttest(meanSpeedUse(:,1),meanSpeedUse(:,4));
end

length(find(pValuesShuffle(:,1)<=0.05))/NShuffle
length(find(pValuesShuffle(:,2)<=0.05))/NShuffle
length(find(pValuesShuffle(:,3)<=0.05))/NShuffle
length(find(pValuesShuffle(:,4)<=0.05))/NShuffle

allDownSamplePValues{idx}=pValuesShuffle;

idx=3;%downsample GFP CS
NTotal=size(meanSpeedAll{idx},1);
pValuesShuffle=[];
for n=1:NShuffle;
    i=randperm(NTotal);
    i=i(1:NTake);
    meanSpeedUse=meanSpeedAll{idx}(i,:);
    [~,pValuesShuffle(n,1)]=ttest(meanSpeedUse(:,1),meanSpeedUse(:,2));
[~,pValuesShuffle(n,2)]=ttest(meanSpeedUse(:,1),meanSpeedUse(:,3));
[~,pValuesShuffle(n,3)]=ttest(meanSpeedUse(:,2),meanSpeedUse(:,3));
[~,pValuesShuffle(n,4)]=ttest(meanSpeedUse(:,1),meanSpeedUse(:,4));
end

length(find(pValuesShuffle(:,1)<=0.05))/NShuffle
length(find(pValuesShuffle(:,2)<=0.05))/NShuffle
length(find(pValuesShuffle(:,3)<=0.05))/NShuffle
length(find(pValuesShuffle(:,4)<=0.05))/NShuffle

allDownSamplePValues{idx}=pValuesShuffle;

idx=4;%downsample GFP RS
NTotal=size(meanSpeedAll{idx},1);
pValuesShuffle=[];
for n=1:NShuffle;
    i=randperm(NTotal);
    i=i(1:NTake);
    meanSpeedUse=meanSpeedAll{idx}(i,:);
    [~,pValuesShuffle(n,1)]=ttest(meanSpeedUse(:,1),meanSpeedUse(:,2));
[~,pValuesShuffle(n,2)]=ttest(meanSpeedUse(:,1),meanSpeedUse(:,3));
[~,pValuesShuffle(n,3)]=ttest(meanSpeedUse(:,2),meanSpeedUse(:,3));
[~,pValuesShuffle(n,4)]=ttest(meanSpeedUse(:,1),meanSpeedUse(:,4));
end

length(find(pValuesShuffle(:,1)<=0.05))/NShuffle
length(find(pValuesShuffle(:,2)<=0.05))/NShuffle
length(find(pValuesShuffle(:,3)<=0.05))/NShuffle
length(find(pValuesShuffle(:,4)<=0.05))/NShuffle

allDownSamplePValues{idx}=pValuesShuffle;

   
save('allDownSamplePValues.mat','allDownSamplePValues')

%% plot p value on a log scale
load('allDownSamplePValues.mat')
NShuffle=1000;
figure,
subplot(121)
p=allDownSamplePValues{1}(:,1);
semilogy(p, 'k.');  % semilogy uses log scale on the y-axis
hold on
semilogy(0.05*ones(1,length(p)), 'r-');  % semilogy uses log scale on the y-axis

xlabel('Index');
ylabel('p-value (log scale)');

NP=length(find(p<=0.05))/NShuffle;

title(['P-values on log scale',num2str(NP)]);
grid on;


subplot(122)
p=allDownSamplePValues{1}(:,2);
semilogy(p, 'k.');  % semilogy uses log scale on the y-axis
hold on
semilogy(0.05*ones(1,length(p)), 'r-');  % semilogy uses log scale on the y-axis

xlabel('Index');
ylabel('p-value (log scale)');
NP=length(find(p<=0.05))/NShuffle;

title(['P-values on log scale',num2str(NP)]);
grid on;

saveas(gcf,'pValuesForDownsampleCHRRCS.fig')
print -painters -depsc pValuesForDownsampleCHRRCS.eps

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
exp=10;
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
zones{1}=[1:1:beforeAfterZone];
zones{2}=[beforeAfterZone+1:1:beforeAfterZone+15];
zones{3}=[beforeAfterZone+15+1:1:beforeAfterZone+15+exp];
zones{4}=[beforeAfterZone+15+exp+1:1:size(speed,2)];

meanSpeed=[];
meanSpeed=[nanmean(speed(:,zones{1}),2) nanmean(speed(:,zones{2}),2) nanmean(speed(:,zones{3}),2) nanmean(speed(:,zones{4}),2)];
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

if N==1 | N==2;
ylim([25 44])
else
    ylim([28 44]);
end



subplot(4,2,N*2);
bar([1:1:4],nanmean(meanSpeed,1),'FaceColor','w');
hold on
errorbar([1:1:4],nanmean(meanSpeed,1),nansem(meanSpeed,1),'k')
ylim([0 45])
[~,p1]=ttest(meanSpeed(:,1),meanSpeed(:,2));
[~,p2]=ttest(meanSpeed(:,1),meanSpeed(:,3));
[~,p3]=ttest(meanSpeed(:,2),meanSpeed(:,3));
[~,p4]=ttest(meanSpeed(:,1),meanSpeed(:,4));
title(num2str([p1 p2 p3 p4]))
end
save('meanSpeedAll.mat','meanSpeedAll')
saveas(gcf,'aligniningStimZoneLowerLickFarApart50cmExp_combineStimDelay.fig')
print -painters -depsc aligniningStimZoneLowerLickFarApart50cmExp_combineStimDelay.eps

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
exp=10;
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
zones{1}=[1:1:beforeAfterZone];
zones{2}=[beforeAfterZone+1:1:beforeAfterZone+15];
zones{3}=[beforeAfterZone+15+1:1:beforeAfterZone+15+exp];
zones{4}=[beforeAfterZone+15+exp+1:1:size(speed,2)];
% zones{4}=[(beforeAfterZone-exp+1):1:beforeAfterZone];
meanSpeed=[];
meanSpeed=[nanmean(speed(:,zones{1}),2) nanmean(speed(:,zones{2}),2) nanmean(speed(:,zones{3}),2) nanmean(speed(:,zones{4}),2)];
NPoints(N,1)=size(meanSpeed,1);
meanSpeedAll{N}=meanSpeed;

i=randperm(size(speed,1));
if length(i)>=61;
i=i(1:61);
end
speed=speed(i,:);
meanSpeed=meanSpeed(i,:);
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

if N==1 | N==2;
ylim([25 42])
else
    ylim([28 44]);
end



subplot(4,2,N*2);
bar([1:1:4],nanmean(meanSpeed,1),'FaceColor','w');
hold on
errorbar([1:1:4],nanmean(meanSpeed,1),nansem(meanSpeed,1),'k')
ylim([0 45])
[~,p1]=ttest(meanSpeed(:,1),meanSpeed(:,2));
[~,p2]=ttest(meanSpeed(:,1),meanSpeed(:,3));
[~,p3]=ttest(meanSpeed(:,2),meanSpeed(:,3));
[~,p4]=ttest(meanSpeed(:,1),meanSpeed(:,4));
title(num2str([p1 p2 p3 p4]))
end

saveas(gcf,'aligniningStimZoneLowerLickFarApart50cmExpDownsample_combineStimDelay.fig')
print -painters -depsc aligniningStimZoneLowerLickFarApart50cmExpDownsample_combineStimDelay.eps

%%  whether this is because the cs data had more points
load('meanSpeedAll.mat')

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
    [~,pValuesShuffle(n,1)]=ttest(meanSpeedUse(:,1),meanSpeedUse(:,2));
[~,pValuesShuffle(n,2)]=ttest(meanSpeedUse(:,1),meanSpeedUse(:,3));
[~,pValuesShuffle(n,3)]=ttest(meanSpeedUse(:,2),meanSpeedUse(:,3));
[~,pValuesShuffle(n,4)]=ttest(meanSpeedUse(:,1),meanSpeedUse(:,4));

meanCompareShuffle(n,1)=nanmean(meanSpeedUse(:,2),1)-nanmean(meanSpeedUse(:,1),1);
meanCompareShuffle(n,2)=nanmean(meanSpeedUse(:,3),1)-nanmean(meanSpeedUse(:,1),1);
meanCompareShuffle(n,3)=nanmean(meanSpeedUse(:,3),1)-nanmean(meanSpeedUse(:,2),1);
meanCompareShuffle(n,4)=nanmean(meanSpeedUse(:,4),1)-nanmean(meanSpeedUse(:,1),1);

end

% length(find(pValuesShuffle(:,1)<=0.05))/NShuffle
% length(find(pValuesShuffle(:,2)<=0.05))/NShuffle
% length(find(pValuesShuffle(:,3)<=0.05))/NShuffle
% length(find(pValuesShuffle(:,4)<=0.05))/NShuffle

allDownSamplePValues{idx}=pValuesShuffle;
meanCompare{idx}=meanCompareShuffle;

idx=3;%downsample GFP CS
NTotal=size(meanSpeedAll{idx},1);
pValuesShuffle=[];
for n=1:NShuffle;
    i=randperm(NTotal);
    i=i(1:NTake);
    meanSpeedUse=meanSpeedAll{idx}(i,:);
    [~,pValuesShuffle(n,1)]=ttest(meanSpeedUse(:,1),meanSpeedUse(:,2));
[~,pValuesShuffle(n,2)]=ttest(meanSpeedUse(:,1),meanSpeedUse(:,3));
[~,pValuesShuffle(n,3)]=ttest(meanSpeedUse(:,2),meanSpeedUse(:,3));
[~,pValuesShuffle(n,4)]=ttest(meanSpeedUse(:,1),meanSpeedUse(:,4));

meanCompareShuffle(n,1)=nanmean(meanSpeedUse(:,2),1)-nanmean(meanSpeedUse(:,1),1);
meanCompareShuffle(n,2)=nanmean(meanSpeedUse(:,3),1)-nanmean(meanSpeedUse(:,1),1);
meanCompareShuffle(n,3)=nanmean(meanSpeedUse(:,3),1)-nanmean(meanSpeedUse(:,2),1);
meanCompareShuffle(n,4)=nanmean(meanSpeedUse(:,4),1)-nanmean(meanSpeedUse(:,1),1);

end


% length(find(pValuesShuffle(:,1)<=0.05))/NShuffle
% length(find(pValuesShuffle(:,2)<=0.05))/NShuffle
% length(find(pValuesShuffle(:,3)<=0.05))/NShuffle
% length(find(pValuesShuffle(:,4)<=0.05))/NShuffle

allDownSamplePValues{idx}=pValuesShuffle;
meanCompare{idx}=meanCompareShuffle;

% idx=4;%downsample GFP RS
% NTotal=size(meanSpeedAll{idx},1);
% pValuesShuffle=[];
% for n=1:NShuffle;
%     i=randperm(NTotal);
%     i=i(1:NTake);
%     meanSpeedUse=meanSpeedAll{idx}(i,:);
%     [~,pValuesShuffle(n,1)]=ttest(meanSpeedUse(:,1),meanSpeedUse(:,2));
% [~,pValuesShuffle(n,2)]=ttest(meanSpeedUse(:,1),meanSpeedUse(:,3));
% [~,pValuesShuffle(n,3)]=ttest(meanSpeedUse(:,2),meanSpeedUse(:,3));
% [~,pValuesShuffle(n,4)]=ttest(meanSpeedUse(:,1),meanSpeedUse(:,4));
% meanCompareShuffle(n,1)=nanmean(meanSpeedUse(:,2),1)-nanmean(meanSpeedUse(:,1),1);
% meanCompareShuffle(n,2)=nanmean(meanSpeedUse(:,3),1)-nanmean(meanSpeedUse(:,1),1);
% meanCompareShuffle(n,3)=nanmean(meanSpeedUse(:,3),1)-nanmean(meanSpeedUse(:,2),1);
% meanCompareShuffle(n,4)=nanmean(meanSpeedUse(:,4),1)-nanmean(meanSpeedUse(:,1),1);
% end

% length(find(pValuesShuffle(:,1)<=0.05))/NShuffle
% length(find(pValuesShuffle(:,2)<=0.05))/NShuffle
% length(find(pValuesShuffle(:,3)<=0.05))/NShuffle
% length(find(pValuesShuffle(:,4)<=0.05))/NShuffle

allDownSamplePValues{idx}=pValuesShuffle;
meanCompare{idx}=meanCompareShuffle;

RSLevelBeforeVSStim=nanmean(meanSpeedAll{2}(:,2))-nanmean(meanSpeedAll{2}(:,1));
RSLevelBeforeVSDelay=nanmean(meanSpeedAll{2}(:,3))-nanmean(meanSpeedAll{2}(:,1));

   
save('allDownSamplePValuesCombineStimDelay.mat','allDownSamplePValues','meanCompare')

%
%plot p value on a log scale
load('allDownSamplePValuesCombineStimDelay.mat')
NShuffle=1000;
figure,
subplot(121)
p=allDownSamplePValues{1}(:,1);
semilogy(p, 'k.');  % semilogy uses log scale on the y-axis
hold on
semilogy(0.05*ones(1,length(p)), 'r-');  % semilogy uses log scale on the y-axis

xlabel('Index');
ylabel('p-value (log scale)');

NP=length(find(p<=0.05))/NShuffle;

title(['P-values on log scale',num2str(NP)]);
grid on;


subplot(122)
p=allDownSamplePValues{1}(:,2);
semilogy(p, 'k.');  % semilogy uses log scale on the y-axis
hold on
semilogy(0.05*ones(1,length(p)), 'r-');  % semilogy uses log scale on the y-axis

xlabel('Index');
ylabel('p-value (log scale)');
NP=length(find(p<=0.05))/NShuffle;

title(['P-values on log scale',num2str(NP)]);
grid on;

saveas(gcf,'pValuesForDownsampleCHRRCSCombineStimDelay.fig')
print -painters -depsc pValuesForDownsampleCHRRCSCombineStimDelay.eps

figure
 subplot(121)
p=allDownSamplePValues{1}(:,1);
m=meanCompare{1}(:,1);
semilogy(m,p, 'k.');  % semilogy uses log scale on the y-axis
hold on
semilogy([min(m) max(m)],[0.05 0.05], 'r-');  % semilogy uses log scale on the y-axis
hold on
semilogy([0 0],[min(p) max(p)], 'r-');  % semilogy uses log scale on the y-axis

hold on
semilogy([RSLevelBeforeVSStim RSLevelBeforeVSStim],[min(p) max(p)], 'g-');  % semilogy uses log scale on the y-axis


p1=length(find(p<=0.05))/1000;
p2=length(find(m<0)/NShuffle);
p3=length(find(m<RSLevelBeforeVSStim)/NShuffle);
title([num2str([p1 p2 p3])]);
axis square


 subplot(122)
p=allDownSamplePValues{1}(:,2);
m=meanCompare{1}(:,2);
semilogy(m,p, 'k.');  % semilogy uses log scale on the y-axis
hold on
semilogy([min(m) max(m)],[0.05 0.05], 'r-');  % semilogy uses log scale on the y-axis
hold on
semilogy([0 0],[min(p) max(p)], 'r-');  % semilogy uses log scale on the y-axis

hold on
semilogy([RSLevelBeforeVSDelay RSLevelBeforeVSDelay],[min(p) max(p)], 'g-');  % semilogy uses log scale on the y-axis


p1=length(find(p<=0.05))/1000;
p2=length(find(m<0)/NShuffle);
p3=length(find(m<RSLevelBeforeVSDelay)/NShuffle);

title([num2str([p1 p2 p3])]);
axis square
saveas(gcf,'pValuesForDownsampleCHRRCSCombineStimDelay_both.fig')
print -painters -depsc pValuesForDownsampleCHRRCSCombineStimDelay_both.eps