load('Z:\papers\2025RewardOptogenetics\speedExamples\cue3_15\stimPositionAll.mat')

load('threshold.mat')
beforeDistance=15;
afterDistance=30;
edges=linspace(0,366,81);
width=nanmean(diff(edges));
xvalues=edges(1:end-1)+width/2;
cues=[63 85;163 185;279 301;374 396];
binWidth=1;
rewardLoc=366/binWidth;
% stimExp=7.5; %expanding the stimulation zone to both sides
% stimExp=5; %expanding the stimulation zone to both sides
% stimExp=6; %good if expansion is symmetrical
stimExp1=12;
stimExp2=29;

% stimExp1=-1;
% stimExp2=-1;

samplingSpacing=1;%how fine to sample out stim area
load('idxCHR2CS');
load('idxCHR2RS');
load('idxGFPCS');
load('idxGFPRS');
%% chr2CS
idx=idxCHR2CS;
stimLocAll=stimPositionAll{1};
%%%during stim
load('lickSegAllData.mat')
load('rewardSegAllData.mat')
load('positionSegAllData.mat')
[lickUsePositionStimRBR,~,~,~] = absPreRewardLicking_RBRRBRLick(lickSegAllData(idx),rewardSegAllData(idx),positionSegAllData(idx),threshold(idx),beforeDistance,afterDistance);
%%%before
load('lickSegAllDataBefore.mat')
load('rewardSegAllDataBefore.mat')
load('positionSegAllDataBefore.mat')
[lickUsePositionBeforeRBR,~,~,~] = absPreRewardLicking_RBRRBRLick(lickSegAllDataBefore(idx),rewardSegAllDataBefore(idx),positionSegAllDataBefore(idx),threshold(idx),beforeDistance,afterDistance);
%%%after
load('lickSegAllDataAfter.mat')
load('rewardSegAllDataAfter.mat')
load('positionSegAllDataAfter.mat')
[lickUsePositionAfterRBR,~,~,~] = absPreRewardLicking_RBRRBRLick(lickSegAllDataAfter(idx),rewardSegAllDataAfter(idx),positionSegAllDataAfter(idx),threshold(idx),beforeDistance,afterDistance);
[lickPercentileBefore,lickPercentileStim,lickPercentileAfter] = lickPercCueStimRS(lickUsePositionBeforeRBR,lickUsePositionStimRBR,lickUsePositionAfterRBR,stimLocAll,stimExp1,stimExp2,samplingSpacing);
save('LickPercentNearStimCHR2CSRBRStimLoc.mat','lickPercentileBefore','lickPercentileStim','lickPercentileAfter')

%% chr2RS
idx=idxCHR2RS;
stimLocAll=stimPositionAll{2};
%%%during stim
load('lickSegAllData.mat')
load('rewardSegAllData.mat')
load('positionSegAllData.mat')
[lickUsePositionStimRBR,~,~,~] = absPreRewardLicking_RBRRBRLick(lickSegAllData(idx),rewardSegAllData(idx),positionSegAllData(idx),threshold(idx),beforeDistance,afterDistance);
%%%before
load('lickSegAllDataBefore.mat')
load('rewardSegAllDataBefore.mat')
load('positionSegAllDataBefore.mat')
[lickUsePositionBeforeRBR,~,~,~] = absPreRewardLicking_RBRRBRLick(lickSegAllDataBefore(idx),rewardSegAllDataBefore(idx),positionSegAllDataBefore(idx),threshold(idx),beforeDistance,afterDistance);
%%%after
load('lickSegAllDataAfter.mat')
load('rewardSegAllDataAfter.mat')
load('positionSegAllDataAfter.mat')
[lickUsePositionAfterRBR,~,~,~] = absPreRewardLicking_RBRRBRLick(lickSegAllDataAfter(idx),rewardSegAllDataAfter(idx),positionSegAllDataAfter(idx),threshold(idx),beforeDistance,afterDistance);
[lickPercentileBefore,lickPercentileStim,lickPercentileAfter] = lickPercCueStimRS(lickUsePositionBeforeRBR,lickUsePositionStimRBR,lickUsePositionAfterRBR,stimLocAll,stimExp1,stimExp2,samplingSpacing);
save('LickPercentNearStimCHR2RSRBRStimLoc.mat','lickPercentileBefore','lickPercentileStim','lickPercentileAfter')
%% GFPCS
idx=idxGFPCS;
stimLocAll=stimPositionAll{3};
%%%during stim
load('lickSegAllData.mat')
load('rewardSegAllData.mat')
load('positionSegAllData.mat')
[lickUsePositionStimRBR,~,~,~] = absPreRewardLicking_RBRRBRLick(lickSegAllData(idx),rewardSegAllData(idx),positionSegAllData(idx),threshold(idx),beforeDistance,afterDistance);
%%%before
load('lickSegAllDataBefore.mat')
load('rewardSegAllDataBefore.mat')
load('positionSegAllDataBefore.mat')
[lickUsePositionBeforeRBR,~,~,~] = absPreRewardLicking_RBRRBRLick(lickSegAllDataBefore(idx),rewardSegAllDataBefore(idx),positionSegAllDataBefore(idx),threshold(idx),beforeDistance,afterDistance);
%%%after
load('lickSegAllDataAfter.mat')
load('rewardSegAllDataAfter.mat')
load('positionSegAllDataAfter.mat')
[lickUsePositionAfterRBR,~,~,~] = absPreRewardLicking_RBRRBRLick(lickSegAllDataAfter(idx),rewardSegAllDataAfter(idx),positionSegAllDataAfter(idx),threshold(idx),beforeDistance,afterDistance);
[lickPercentileBefore,lickPercentileStim,lickPercentileAfter] = lickPercCueStimRS(lickUsePositionBeforeRBR,lickUsePositionStimRBR,lickUsePositionAfterRBR,stimLocAll,stimExp1,stimExp2,samplingSpacing);
save('LickPercentNearStimGFPCSRBRStimLoc.mat','lickPercentileBefore','lickPercentileStim','lickPercentileAfter')
%% GFPRS
idx=idxGFPRS;
stimLocAll=stimPositionAll{4};
%%%during stim
load('lickSegAllData.mat')
load('rewardSegAllData.mat')
load('positionSegAllData.mat')
[lickUsePositionStimRBR,~,~,~] = absPreRewardLicking_RBRRBRLick(lickSegAllData(idx),rewardSegAllData(idx),positionSegAllData(idx),threshold(idx),beforeDistance,afterDistance);
%%%before
load('lickSegAllDataBefore.mat')
load('rewardSegAllDataBefore.mat')
load('positionSegAllDataBefore.mat')
[lickUsePositionBeforeRBR,~,~,~] = absPreRewardLicking_RBRRBRLick(lickSegAllDataBefore(idx),rewardSegAllDataBefore(idx),positionSegAllDataBefore(idx),threshold(idx),beforeDistance,afterDistance);
%%%after
load('lickSegAllDataAfter.mat')
load('rewardSegAllDataAfter.mat')
load('positionSegAllDataAfter.mat')
[lickUsePositionAfterRBR,~,~,~] = absPreRewardLicking_RBRRBRLick(lickSegAllDataAfter(idx),rewardSegAllDataAfter(idx),positionSegAllDataAfter(idx),threshold(idx),beforeDistance,afterDistance);
[lickPercentileBefore,lickPercentileStim,lickPercentileAfter] = lickPercCueStimRS(lickUsePositionBeforeRBR,lickUsePositionStimRBR,lickUsePositionAfterRBR,stimLocAll,stimExp1,stimExp2,samplingSpacing);
save('LickPercentNearStimGFPRSRBRStimLoc.mat','lickPercentileBefore','lickPercentileStim','lickPercentileAfter')


%% plot per session
sessionDataAll={};
load('LickPercentNearStimCHR2CSRBRStimLoc.mat')
sessionDataAll{1}=[lickPercentileBefore' lickPercentileStim' lickPercentileAfter'];
load('LickPercentNearStimCHR2RSRBRStimLoc.mat')
sessionDataAll{2}=[lickPercentileBefore' lickPercentileStim' lickPercentileAfter'];
load('LickPercentNearStimGFPCSRBRStimLoc.mat')
sessionDataAll{3}=[lickPercentileBefore' lickPercentileStim' lickPercentileAfter'];
load('LickPercentNearStimGFPRSRBRStimLoc.mat')
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

saveas(gcf,'LickPercentNearStimPerSessionRBRStimLoc.fig')

%% 
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
saveas(gcf,'LickPercentNearStimRBR_RBSSTimLoc.fig')

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
saveas(gcf,'LickPercentNearStimRBRZero_RBSSTimLoc.fig')

