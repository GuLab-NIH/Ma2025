beforeDistance=20;
afterDistance=30;
% 
% beforeDistance=10;
% afterDistance=25;
% load('Z:\papers\2025RewardOptogenetics\optRBR\dual\allIdxAll.mat')
%% stim
load('lickSegAllData.mat');
load('rewardSegAllData.mat');
load('positionSegAllData.mat');
load('threshold.mat');

[predLickStimRBR,predLickStimAllRuns] = predLickingPercentileRollingWindow_RBR(lickSegAllData,rewardSegAllData,positionSegAllData,threshold,beforeDistance,afterDistance);
save('predLickStimRBR.mat','predLickStimRBR');
save('predLickStimAllRuns.mat','predLickStimAllRuns');

%% 10 runs before
load('lickSegAllDataBefore.mat');
load('rewardSegAllDataBefore.mat');
load('positionSegAllDataBefore.mat');
load('threshold.mat');

[predLickBeforeRBR,predLickBeforeAllRuns] = predLickingPercentileRollingWindow_RBR(lickSegAllDataBefore,rewardSegAllDataBefore,positionSegAllDataBefore,threshold,beforeDistance,afterDistance);
save('predLickBeforeRBR.mat','predLickBeforeRBR');
save('predLickBeforeAllRuns.mat','predLickBeforeAllRuns');
%% 10 runs after
load('lickSegAllDataAfter.mat');
load('rewardSegAllDataAfter.mat');
load('positionSegAllDataAfter.mat');
load('threshold.mat');

[predLickAfterRBR,predLickAfterAllRuns] = predLickingPercentileRollingWindow_RBR(lickSegAllDataAfter,rewardSegAllDataAfter,positionSegAllDataAfter,threshold,beforeDistance,afterDistance);
save('predLickAfterRBR.mat','predLickAfterRBR');
save('predLickAfterAllRuns.mat','predLickAfterAllRuns');
%% 
% idxCHR2RS=[1:30];
% idxCHR2CS=[52:81];
% idxGFPRS=[31:51];
% idxGFPCS=[82:102];
% 
% idxUseCHR2CS=allIdxAll{1};
% idxUseCHR2RS=allIdxAll{2};
% idxUseGFPCS=allIdxAll{3};
% idxUseGFPRS=allIdxAll{4};
% 
% idxCHR2CS=idxCHR2CS(idxUseCHR2CS);
% idxCHR2RS=idxCHR2RS(idxUseCHR2RS);
% idxGFPCS=idxGFPCS(idxUseGFPCS);
% idxGFPRS=idxGFPRS(idxUseGFPRS);
load('idxCHR2RS.mat')
load('idxCHR2CS.mat')
load('idxGFPRS.mat')
load('idxGFPCS.mat')
% idxCHR2RS=[1:15];
% allIdxAll{1}=setdiff([1:1:length(idxCHR2RS)],[1 2 3 10 8 13]);%CHR2RS: remove 8 and 13: too many lick signals per bin
% idxCHR2RS=idxCHR2RS(allIdxAll{1});

%% CHR2 CS
allRBR=[];
idx=idxCHR2CS;
allRBR=[predLickBeforeRBR(idx,:) predLickStimRBR(idx,:) predLickAfterRBR(idx,:)];
allAll=[predLickBeforeAllRuns(idx,:) predLickStimAllRuns(idx,:) predLickAfterAllRuns(idx,:)];

figure
subplot(241)
errorbar([1:1:size(allRBR,2)],nanmean(allRBR,1),nansem(allRBR,1));
hold on
plot([11:1:20],nanmean(allRBR(:,[11:1:20]),1),'ro')
before=allRBR(:,1:10);
before=reshape(before,[size(before,1)*size(before,2),1]);
stim=allRBR(:,11:20);
stim=reshape(stim,[size(stim,1)*size(stim,2),1]);

[~,p]=ttest2(before,stim);
title(['chr2cs RBR',num2str(p)])

subplot(245)
errorbar([1:1:size(allAll,2)],nanmean(allAll,1),nansem(allAll,1));
hold on
plot(2,nanmean(allAll(:,2),1),'ro')
[~,p]=ttest(allAll(:,1),allAll(:,2));
title(['chr2cs all runs',num2str(p)])

%% CHR2 RS
allRBR=[];
idx=idxCHR2RS;
allRBR=[predLickBeforeRBR(idx,:) predLickStimRBR(idx,:) predLickAfterRBR(idx,:)];
allAll=[predLickBeforeAllRuns(idx,:) predLickStimAllRuns(idx,:) predLickAfterAllRuns(idx,:)];


subplot(242)
errorbar([1:1:size(allRBR,2)],nanmean(allRBR,1),nansem(allRBR,1));
hold on
plot([11:1:20],nanmean(allRBR(:,[11:1:20]),1),'ro')
before=allRBR(:,1:10);
before=reshape(before,[size(before,1)*size(before,2),1]);
stim=allRBR(:,11:20);
stim=reshape(stim,[size(stim,1)*size(stim,2),1]);

[~,p]=ttest2(before,stim);
title(['chr2rs RBR',num2str(p)])

subplot(246)
errorbar([1:1:size(allAll,2)],nanmean(allAll,1),nansem(allAll,1));
hold on
plot(2,nanmean(allAll(:,2),1),'ro')
[~,p]=ttest(allAll(:,1),allAll(:,2));
title(['chr2rs all runs',num2str(p)])

%% GFP CS
allRBR=[];
idx=idxGFPCS;
allRBR=[predLickBeforeRBR(idx,:) predLickStimRBR(idx,:) predLickAfterRBR(idx,:)];
allAll=[predLickBeforeAllRuns(idx,:) predLickStimAllRuns(idx,:) predLickAfterAllRuns(idx,:)];

subplot(243)
errorbar([1:1:size(allRBR,2)],nanmean(allRBR,1),nansem(allRBR,1));
hold on
plot([11:1:20],nanmean(allRBR(:,[11:1:20]),1),'ro')
before=allRBR(:,1:10);
before=reshape(before,[size(before,1)*size(before,2),1]);
stim=allRBR(:,11:20);
stim=reshape(stim,[size(stim,1)*size(stim,2),1]);

[~,p]=ttest2(before,stim);
title(['gfpcs RBR',num2str(p)])

subplot(247)
errorbar([1:1:size(allAll,2)],nanmean(allAll,1),nansem(allAll,1));
hold on
plot(2,nanmean(allAll(:,2),1),'ro')
[~,p]=ttest(allAll(:,1),allAll(:,2));
title(['gfpcs all runs',num2str(p)])

%% GFP RS
allRBR=[];
idx=idxGFPRS;
allRBR=[predLickBeforeRBR(idx,:) predLickStimRBR(idx,:) predLickAfterRBR(idx,:)];
allAll=[predLickBeforeAllRuns(idx,:) predLickStimAllRuns(idx,:) predLickAfterAllRuns(idx,:)];


subplot(244)
errorbar([1:1:size(allRBR,2)],nanmean(allRBR,1),nansem(allRBR,1));
hold on
plot([11:1:20],nanmean(allRBR(:,[11:1:20]),1),'ro')
before=allRBR(:,1:10);
before=reshape(before,[size(before,1)*size(before,2),1]);
stim=allRBR(:,11:20);
stim=reshape(stim,[size(stim,1)*size(stim,2),1]);

[~,p]=ttest2(before,stim);
title(['gfprs RBR',num2str(p)])


subplot(248)
errorbar([1:1:size(allAll,2)],nanmean(allAll,1),nansem(allAll,1));
hold on
plot(2,nanmean(allAll(:,2),1),'ro')
[~,p]=ttest(allAll(:,1),allAll(:,2));
title(['gfprs all runs',num2str(p)])

saveas(gcf,'predictiveLickingPercentileRollingWindow.fig')
%% bar graph per session
idxAll={};
idxAll{1}=idxCHR2CS;
idxAll{2}=idxCHR2RS;
idxAll{3}=idxGFPCS;
idxAll{4}=idxGFPRS;

figure
for t=1:length(idxAll)
allRBR=[];
idx=idxAll{t};

allRBR=[predLickBeforeRBR(idx,:) predLickStimRBR(idx,:) predLickAfterRBR(idx,:)];

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
saveas(gcf,'predictiveLickingPercentileRollingWindowPerSession.fig')
%% RBR shade

idx=idxCHR2CS;
RBRCSadj=[predLickBeforeRBR(idx,:) predLickStimRBR(idx,:) predLickAfterRBR(idx,:)];

idx=idxCHR2RS;
RBRadj=[predLickBeforeRBR(idx,:) predLickStimRBR(idx,:) predLickAfterRBR(idx,:)];

idx=idxGFPCS;
RBRCSGFPadj=[predLickBeforeRBR(idx,:) predLickStimRBR(idx,:) predLickAfterRBR(idx,:)];

idx=idxGFPRS;
RBRRSGFPadj=[predLickBeforeRBR(idx,:) predLickStimRBR(idx,:) predLickAfterRBR(idx,:)];

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
%
[~,p]=ttest2(reshape(A(:,[11:20]),[1 size(A,1)*10]),reshape(B(:,[11:20]),[1 size(B,1)*10]));

 % [pAnova,pMC,pLabel] = anovaRM2W_full_BH(A(:,11:20),B(:,11:20),1);
 % [pAnova2,pMC,pLabel] = anovaRM2W_full_BH(A(:,1:10),B(:,1:10),1);
 % [pAnova3,pMC,pLabel] = anovaRM2W_full_BH(A(:,21:30),B(:,21:30),1);



% title(['chr2 CS RS,P=',num2str([pAnova(1) pAnova2(1) pAnova3(1)])])
title(['chr2 CS RS',num2str(p)])

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

[~,p]=ttest2(reshape(A(:,[11:20]),[1 size(A,1)*10]),reshape(B(:,[11:20]),[1 size(B,1)*10]));

% ylim([0.6 1.3])
% title(['gfp CS RS, P=',num2str([pAnova(1) pAnova2(1) pAnova3(1)])])
title(['gfp CS RS',num2str(p)])

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

[~,p]=ttest2(reshape(A(:,[11:20]),[1 size(A,1)*10]),reshape(B(:,[11:20]),[1 size(B,1)*10]));

% ylim([0.6 1.3])
% title(['CS CHR2 GFP, P=',num2str([pAnova(1) pAnova2(1) pAnova3(1)])])
title(['CS CHR2 GFP',num2str(p)])

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

[~,p]=ttest2(reshape(A(:,[11:20]),[1 size(A,1)*10]),reshape(B(:,[11:20]),[1 size(B,1)*10]));

% ylim([0.6 1.3])
% title(['RS CHR2 GFP, P=',num2str([pAnova(1) pAnova2(1) pAnova3(1)])])
title(['RS CHR2 GFP',num2str(p)])
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
saveas(gcf,'RBRTogether.fig')

%% RBR shade: normalized to zero

idx=idxCHR2CS;
RBRCSadj=[predLickBeforeRBR(idx,:) predLickStimRBR(idx,:) predLickAfterRBR(idx,:)];

idx=idxCHR2RS;
RBRadj=[predLickBeforeRBR(idx,:) predLickStimRBR(idx,:) predLickAfterRBR(idx,:)];

idx=idxGFPCS;
RBRCSGFPadj=[predLickBeforeRBR(idx,:) predLickStimRBR(idx,:) predLickAfterRBR(idx,:)];

idx=idxGFPRS;
RBRRSGFPadj=[predLickBeforeRBR(idx,:) predLickStimRBR(idx,:) predLickAfterRBR(idx,:)];

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
%
[~,p1]=ttest2(reshape(A(:,[11:20]),[1 size(A,1)*10]),reshape(B(:,[11:20]),[1 size(B,1)*10]));
[~,p2]=ttest2(reshape(A(:,[1:10]),[1 size(A,1)*10]),reshape(B(:,[1:10]),[1 size(B,1)*10]));
[~,p3]=ttest2(reshape(A(:,[21:30]),[1 size(A,1)*10]),reshape(B(:,[21:30]),[1 size(B,1)*10]));
A2=reshape(A(:,[11:20]),[1 size(A,1)*10]);
A1=reshape(A(:,[1:10]),[1 size(A,1)*10]);

B2=reshape(B(:,[11:20]),[1 size(B,1)*10]);
B1=reshape(B(:,[1:10]),[1 size(B,1)*10]);

[~,pPreDurRS]=ttest2(A2,A1)
[~,pPreDurCS]=ttest2(B2,B1)

 % [pAnova,pMC,pLabel] = anovaRM2W_full_BH(A(:,11:20),B(:,11:20),1);
 % [pAnova2,pMC,pLabel] = anovaRM2W_full_BH(A(:,1:10),B(:,1:10),1);
 % [pAnova3,pMC,pLabel] = anovaRM2W_full_BH(A(:,21:30),B(:,21:30),1);



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
ylim([-0.8 0.6])

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
saveas(gcf,'RBRTogetherZero.fig')



%% below: pick a subset of RS runs so that the pre-stim levels are simimlar to that of CS 


%% bar graph per session
idxAll={};
idxAll{1}=idxCHR2CS;
idxAll{2}=idxCHR2RS;
idxAll{3}=idxGFPCS;
idxAll{4}=idxGFPRS;

figure
for t=1:length(idxAll)
allRBR=[];
idx=idxAll{t};

allRBR=[predLickBeforeRBR(idx,:) predLickStimRBR(idx,:) predLickAfterRBR(idx,:)];

allSession=[];
allSession(:,1)=nanmean(allRBR(:,[1:10]),2);
allSession(:,2)=nanmean(allRBR(:,[11:20]),2);
allSession(:,3)=nanmean(allRBR(:,[21:30]),2);
if t==2;
i=find(allSession(:,1)~=0);
% i=[4 5 6 7 8 9 10];
% i=[4 6 7 8 9 12];
allSession=allSession(i,:);
idxCHR2RS=idxCHR2RS(i);
idxCHR2RSHigherLick=idxCHR2RS;
save('idxCHR2RSHigherLick.mat','idxCHR2RSHigherLick')
end
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
saveas(gcf,'predictiveLickingPercentileRollingWindowPerSession_adjRSLevel.fig')
%% RBR shade

idx=idxCHR2CS;
RBRCSadj=[predLickBeforeRBR(idx,:) predLickStimRBR(idx,:) predLickAfterRBR(idx,:)];

idx=idxCHR2RS;
RBRadj=[predLickBeforeRBR(idx,:) predLickStimRBR(idx,:) predLickAfterRBR(idx,:)];

idx=idxGFPCS;
RBRCSGFPadj=[predLickBeforeRBR(idx,:) predLickStimRBR(idx,:) predLickAfterRBR(idx,:)];

idx=idxGFPRS;
RBRRSGFPadj=[predLickBeforeRBR(idx,:) predLickStimRBR(idx,:) predLickAfterRBR(idx,:)];

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
%
[~,p]=ttest2(reshape(A(:,[11:20]),[1 size(A,1)*10]),reshape(B(:,[11:20]),[1 size(B,1)*10]));

 % [pAnova,pMC,pLabel] = anovaRM2W_full_BH(A(:,11:20),B(:,11:20),1);
 % [pAnova2,pMC,pLabel] = anovaRM2W_full_BH(A(:,1:10),B(:,1:10),1);
 % [pAnova3,pMC,pLabel] = anovaRM2W_full_BH(A(:,21:30),B(:,21:30),1);



% title(['chr2 CS RS,P=',num2str([pAnova(1) pAnova2(1) pAnova3(1)])])
title(['chr2 CS RS',num2str(p)])

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

[~,p]=ttest2(reshape(A(:,[11:20]),[1 size(A,1)*10]),reshape(B(:,[11:20]),[1 size(B,1)*10]));

% ylim([0.6 1.3])
% title(['gfp CS RS, P=',num2str([pAnova(1) pAnova2(1) pAnova3(1)])])
title(['gfp CS RS',num2str(p)])

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

[~,p]=ttest2(reshape(A(:,[11:20]),[1 size(A,1)*10]),reshape(B(:,[11:20]),[1 size(B,1)*10]));

% ylim([0.6 1.3])
% title(['CS CHR2 GFP, P=',num2str([pAnova(1) pAnova2(1) pAnova3(1)])])
title(['CS CHR2 GFP',num2str(p)])

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

[~,p]=ttest2(reshape(A(:,[11:20]),[1 size(A,1)*10]),reshape(B(:,[11:20]),[1 size(B,1)*10]));

% ylim([0.6 1.3])
% title(['RS CHR2 GFP, P=',num2str([pAnova(1) pAnova2(1) pAnova3(1)])])
title(['RS CHR2 GFP',num2str(p)])
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
saveas(gcf,'RBRTogetherAdjRS.fig')

%% RBR shade: normalized to zero

idx=idxCHR2CS;
RBRCSadj=[predLickBeforeRBR(idx,:) predLickStimRBR(idx,:) predLickAfterRBR(idx,:)];

idx=idxCHR2RS;
RBRadj=[predLickBeforeRBR(idx,:) predLickStimRBR(idx,:) predLickAfterRBR(idx,:)];

idx=idxGFPCS;
RBRCSGFPadj=[predLickBeforeRBR(idx,:) predLickStimRBR(idx,:) predLickAfterRBR(idx,:)];

idx=idxGFPRS;
RBRRSGFPadj=[predLickBeforeRBR(idx,:) predLickStimRBR(idx,:) predLickAfterRBR(idx,:)];

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
%
[~,p1]=ttest2(reshape(A(:,[11:20]),[1 size(A,1)*10]),reshape(B(:,[11:20]),[1 size(B,1)*10]));
[~,p2]=ttest2(reshape(A(:,[1:10]),[1 size(A,1)*10]),reshape(B(:,[1:10]),[1 size(B,1)*10]));
[~,p3]=ttest2(reshape(A(:,[21:30]),[1 size(A,1)*10]),reshape(B(:,[21:30]),[1 size(B,1)*10]));

A2=reshape(A(:,[11:20]),[1 size(A,1)*10]);
A1=reshape(A(:,[1:10]),[1 size(A,1)*10]);

B2=reshape(B(:,[11:20]),[1 size(B,1)*10]);
B1=reshape(B(:,[1:10]),[1 size(B,1)*10]);

[~,pPreDurRS]=ttest2(A2,A1)
[~,pPreDurCS]=ttest2(B2,B1)
 % [pAnova,pMC,pLabel] = anovaRM2W_full_BH(A(:,11:20),B(:,11:20),1);
 % [pAnova2,pMC,pLabel] = anovaRM2W_full_BH(A(:,1:10),B(:,1:10),1);
 % [pAnova3,pMC,pLabel] = anovaRM2W_full_BH(A(:,21:30),B(:,21:30),1);



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
        plot(n,0.2,'k*')
    end
end
% ylim([-0.8 0.6])

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
        plot(n,0.4,'k*')
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
        plot(n,0.4,'k*')
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
        plot(n,0.4,'k*')
    end
end
saveas(gcf,'RBRTogetherZeroAdjRS.fig')