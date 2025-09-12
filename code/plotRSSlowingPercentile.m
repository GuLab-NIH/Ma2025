
%% GET ALL STIM POSITIONS
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

stimPositionAll={}; 
%1: chr2cs, 2. chr2 rs, 3, gfp cs, 4. gfp rs. PER MOUSE. WITHIN: ROW: PER SESSION, COLUMN: PER RUN
for N=1:length(allIdxAll);
mouseIdx=allIdxAll{N};

position=positionSegAllData(mouseIdx);
stim=stimZoneIdxAllDataFix(mouseIdx);
stimPositionAll{N}={};


for n=1:length(position); %per session
    for m=1:length(position{n}) %per run
        stimPositionAll{N}{n,m}=position{n}{m}(stim{n}{m});
    end
end
end
save('stimPositionAll.mat','stimPositionAll')

%%
cues=[63 85;163 185;279 301;374 396];
stimExp1=0;
stimExp2=0;
samplingSpacing=1;


binWidth=1;
rewardLoc=366/binWidth;
xvalues=[1:1:400];

load('allSpeedYUse.mat')
idx=1; %chr2cs
speedBefore=allSPreUse{idx};
speedStim=allSDurUse{idx};
speedAfter=allSPostUse{idx};
yBefore=allYPreUse{idx};
yStim=allYDurUse{idx};
yAfter=allYPostUse{idx};
stimLocAll=stimPositionAll{idx};
[decePercentileBefore,decePercentileStim,decePercentileAfter] = decelerationPercCueStimRS(speedBefore,speedStim,speedAfter,yBefore,yStim,yAfter,stimLocAll,stimExp1,stimExp2,samplingSpacing);
save('decePercentNearStimCHR2CSRBRStimLoc.mat','decePercentileBefore','decePercentileStim','decePercentileAfter')

idx=2; %chr2rs
speedBefore=allSPreUse{idx};
speedStim=allSDurUse{idx};
speedAfter=allSPostUse{idx};
yBefore=allYPreUse{idx};
yStim=allYDurUse{idx};
yAfter=allYPostUse{idx};
stimLocAll=stimPositionAll{idx};
[decePercentileBefore,decePercentileStim,decePercentileAfter] = decelerationPercCueStimRS(speedBefore,speedStim,speedAfter,yBefore,yStim,yAfter,stimLocAll,stimExp1,stimExp2,samplingSpacing);
save('decePercentNearStimCHR2RSRBRStimLoc.mat','decePercentileBefore','decePercentileStim','decePercentileAfter')

idx=3; %GFPcs
speedBefore=allSPreUse{idx};
speedStim=allSDurUse{idx};
speedAfter=allSPostUse{idx};
yBefore=allYPreUse{idx};
yStim=allYDurUse{idx};
yAfter=allYPostUse{idx};
stimLocAll=stimPositionAll{idx};
[decePercentileBefore,decePercentileStim,decePercentileAfter] = decelerationPercCueStimRS(speedBefore,speedStim,speedAfter,yBefore,yStim,yAfter,stimLocAll,stimExp1,stimExp2,samplingSpacing);
save('decePercentNearStimGFPCSRBRStimLoc.mat','decePercentileBefore','decePercentileStim','decePercentileAfter')

idx=4; %GFPrs
speedBefore=allSPreUse{idx};
speedStim=allSDurUse{idx};
speedAfter=allSPostUse{idx};
yBefore=allYPreUse{idx};
yStim=allYDurUse{idx};
yAfter=allYPostUse{idx};
stimLocAll=stimPositionAll{idx};
[decePercentileBefore,decePercentileStim,decePercentileAfter] = decelerationPercCueStimRS(speedBefore,speedStim,speedAfter,yBefore,yStim,yAfter,stimLocAll,stimExp1,stimExp2,samplingSpacing);
save('decePercentNearStimGFPRSRBRStimLoc.mat','decePercentileBefore','decePercentileStim','decePercentileAfter')
%%

sessionDataAll={};
load('decePercentNearStimCHR2CSRBRStimLoc.mat')
sessionDataAll{1}=[decePercentileBefore' decePercentileStim' decePercentileAfter'];
load('decePercentNearStimCHR2RSRBRStimLoc.mat')
sessionDataAll{2}=[decePercentileBefore' decePercentileStim' decePercentileAfter'];
load('decePercentNearStimGFPCSRBRStimLoc.mat')
sessionDataAll{3}=[decePercentileBefore' decePercentileStim' decePercentileAfter'];
load('decePercentNearStimGFPRSRBRStimLoc.mat')
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

saveas(gcf,'decePercentNearStimPerSessionRBRStimLoc.fig')

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
saveas(gcf,'decePercentNearStimRBR_RBSStimLoc.fig')

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
saveas(gcf,'decePercentNearStimRBRZero_RBSStimLoc.fig')
