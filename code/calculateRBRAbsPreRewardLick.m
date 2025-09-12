beforeDistance=15;
afterDistance=30;
% load('Z:\papers\2025RewardOptogenetics\optRBR\dual\allIdxAll.mat')
%% stim
load('lickSegAllData.mat');
load('rewardSegAllData.mat');
load('positionSegAllData.mat');
load('threshold.mat');

[lickUsePositionStimAll,predLickPositionStimAll,absPreRewardLickStimRBR,absPreRewardLickStimAllRuns] = absPreRewardLicking_RBR(lickSegAllData,rewardSegAllData,positionSegAllData,threshold,beforeDistance,afterDistance);
save('absPreRewardLickStimRBR.mat','absPreRewardLickStimRBR')
save('absPreRewardLickStimAllRuns.mat','absPreRewardLickStimAllRuns')
%% 10 runs before
load('lickSegAllDataBefore.mat');
load('rewardSegAllDataBefore.mat');
load('positionSegAllDataBefore.mat');
load('threshold.mat');

[lickUsePositionBeforeAll,predLickPositionBeforeAll,absPreRewardLickBeforeRBR,absPreRewardLickBeforeAllRuns] = absPreRewardLicking_RBR(lickSegAllDataBefore,rewardSegAllDataBefore,positionSegAllDataBefore,threshold,beforeDistance,afterDistance);
save('absPreRewardLickBeforeRBR.mat','absPreRewardLickBeforeRBR')
save('absPreRewardLickBeforeAllRuns.mat','absPreRewardLickBeforeAllRuns')
%% 10 runs after
load('lickSegAllDataAfter.mat');
load('rewardSegAllDataAfter.mat');
load('positionSegAllDataAfter.mat');
load('threshold.mat');

[lickUsePositionAfterAll,predLickPositionAfterAll,absPreRewardLickAfterRBR,absPreRewardLickAfterAllRuns] = absPreRewardLicking_RBR(lickSegAllDataAfter,rewardSegAllDataAfter,positionSegAllDataAfter,threshold,beforeDistance,afterDistance);
save('absPreRewardLickAfterRBR.mat','absPreRewardLickAfterRBR')
save('absPreRewardLickAfterAllRuns.mat','absPreRewardLickAfterAllRuns')
%% 
idxCHR2RS=[1:15];
idxCHR2CS=[27:44];
idxGFPRS=[16:26];
idxGFPCS=[45:56];
allIdxAll={};%REMOVE CELLS WITH EXTREAMLY HIGH LICKS MORE THAN 600 PER 5 CM BIN BASED ON THE CALCULATION BELOW
%LINE 1195: do this:
% imagesc(countStimAllCHR2CS)
% imagesc(countStimAllCHR2RS)
%YOU will see that the cells excluded here have very high licks
% allIdxAll{1}=setdiff([1:1:length(idxCHR2RS)],[1 2 3 10 8 13]);%CHR2RS: remove poor performers because otherwise the licking was too low compared to CS: so people may think that this is not fair control

%criteria to remove: per bin more than 800: in the parameter below:
%countStimAllCHR2CS,countStimAllCHR2RS,countStimAllGFPCS, countStimAllGFPRS
%%%%%%%%%%%%%%%%%%%%
% allIdxAll{2}=setdiff([1:1:length(idxCHR2CS)],[12 14]); %chr2 CS: 
% allIdxAll{1}=setdiff([1:1:length(idxCHR2RS)],[8 9 13]);% chr2RS THESE TWO CELLS HAVE VERY HIGH LICK
% % allIdxAll{2}=setdiff([1:1:length(idxCHR2CS)],[]); %chr2 CS: 
% allIdxAll{3}=setdiff([1:1:length(idxGFPRS)],[]);%GFPRS
% allIdxAll{4}=setdiff([1:1:length(idxGFPCS)],[]); %GFP CS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idxCHR2CS=idxCHR2CS(allIdxAll{2});
idxCHR2RS=idxCHR2RS(allIdxAll{1});
idxGFPCS=idxGFPCS(allIdxAll{4});
idxGFPRS=idxGFPRS(allIdxAll{3});

save('idxCHR2CS.mat','idxCHR2CS')
save('idxCHR2RS.mat','idxCHR2RS')
save('idxGFPCS.mat','idxGFPCS')
save('idxGFPRS.mat','idxGFPRS')
% 
% idxCHR2CS=idxCHR2CS(idxUseCHR2CS);
% idxCHR2RS=idxCHR2RS(idxUseCHR2RS);
% idxGFPCS=idxGFPCS(idxUseGFPCS);
% idxGFPRS=idxGFPRS(idxUseGFPRS);
% save('idxCHR2RS.mat','idxCHR2RS')
% save('idxCHR2CS.mat','idxCHR2CS')
% save('idxGFPRS.mat','idxGFPRS')
% save('idxGFPCS.mat','idxGFPCS')
%% CHR2 C
allRBR=[];
idx=idxCHR2CS;
allRBR=[absPreRewardLickBeforeRBR(idx,:) absPreRewardLickStimRBR(idx,:) absPreRewardLickAfterRBR(idx,:)];
allAll=[absPreRewardLickBeforeAllRuns(idx,:) absPreRewardLickStimAllRuns(idx,:) absPreRewardLickAfterAllRuns(idx,:)];

figure
subplot(241)
errorbar([1:1:size(allRBR,2)],nanmean(allRBR,1),nansem(allRBR,1));
hold on
plot([11:1:20],nanmean(allRBR(:,[11:1:20]),1),'ro')
beforeAfter=allRBR(:,1:10);
beforeAfter=reshape(beforeAfter,[size(beforeAfter,1)*size(beforeAfter,2),1]);
stim=allRBR(:,11:20);
stim=reshape(stim,[size(stim,1)*size(stim,2),1]);
ylim([0 100])
[~,p]=ttest2(beforeAfter,stim);
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
allRBR=[absPreRewardLickBeforeRBR(idx,:) absPreRewardLickStimRBR(idx,:) absPreRewardLickAfterRBR(idx,:)];
allAll=[absPreRewardLickBeforeAllRuns(idx,:) absPreRewardLickStimAllRuns(idx,:) absPreRewardLickAfterAllRuns(idx,:)];


subplot(242)
errorbar([1:1:size(allRBR,2)],nanmean(allRBR,1),nansem(allRBR,1));
hold on
plot([11:1:20],nanmean(allRBR(:,[11:1:20]),1),'ro')
beforeAfter=allRBR(:,[1:10 21:30]);
beforeAfter=reshape(beforeAfter,[size(beforeAfter,1)*size(beforeAfter,2),1]);
stim=allRBR(:,11:20);
stim=reshape(stim,[size(stim,1)*size(stim,2),1]);
ylim([0 100])
[~,p]=ttest2(beforeAfter,stim);
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
allRBR=[absPreRewardLickBeforeRBR(idx,:) absPreRewardLickStimRBR(idx,:) absPreRewardLickAfterRBR(idx,:)];
allAll=[absPreRewardLickBeforeAllRuns(idx,:) absPreRewardLickStimAllRuns(idx,:) absPreRewardLickAfterAllRuns(idx,:)];

subplot(243)
errorbar([1:1:size(allRBR,2)],nanmean(allRBR,1),nansem(allRBR,1));
hold on
plot([11:1:20],nanmean(allRBR(:,[11:1:20]),1),'ro')
beforeAfter=allRBR(:,[1:10 21:30]);
beforeAfter=reshape(beforeAfter,[size(beforeAfter,1)*size(beforeAfter,2),1]);
stim=allRBR(:,11:20);
stim=reshape(stim,[size(stim,1)*size(stim,2),1]);
ylim([0 100])
[~,p]=ttest2(beforeAfter,stim);
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
allRBR=[absPreRewardLickBeforeRBR(idx,:) absPreRewardLickStimRBR(idx,:) absPreRewardLickAfterRBR(idx,:)];
allAll=[absPreRewardLickBeforeAllRuns(idx,:) absPreRewardLickStimAllRuns(idx,:) absPreRewardLickAfterAllRuns(idx,:)];


subplot(244)
errorbar([1:1:size(allRBR,2)],nanmean(allRBR,1),nansem(allRBR,1));
hold on
plot([11:1:20],nanmean(allRBR(:,[11:1:20]),1),'ro')
beforeAfter=allRBR(:,[1:10 21:30]);
beforeAfter=reshape(beforeAfter,[size(beforeAfter,1)*size(beforeAfter,2),1]);
stim=allRBR(:,11:20);
stim=reshape(stim,[size(stim,1)*size(stim,2),1]);
ylim([0 100])
[~,p]=ttest2(beforeAfter,stim);
title(['gfprs RBR',num2str(p)])


subplot(248)
errorbar([1:1:size(allAll,2)],nanmean(allAll,1),nansem(allAll,1));
hold on
plot(2,nanmean(allAll(:,2),1),'ro')
[~,p]=ttest(allAll(:,1),allAll(:,2));
title(['gfprs all runs',num2str(p)])

saveas(gcf,'absPreRewardictiveLicking.fig')

%% plot the shade version of this

%% shading

idx=idxCHR2RS;
RBR=[absPreRewardLickBeforeRBR(idx,:) absPreRewardLickStimRBR(idx,:) absPreRewardLickAfterRBR(idx,:)];

idx=idxCHR2CS;
RBRCS=[absPreRewardLickBeforeRBR(idx,:) absPreRewardLickStimRBR(idx,:) absPreRewardLickAfterRBR(idx,:)];

idx=idxGFPCS;
RBRCSGFP=[absPreRewardLickBeforeRBR(idx,:) absPreRewardLickStimRBR(idx,:) absPreRewardLickAfterRBR(idx,:)];

idx=idxGFPRS;
RBRRSGFP=[absPreRewardLickBeforeRBR(idx,:) absPreRewardLickStimRBR(idx,:) absPreRewardLickAfterRBR(idx,:)];

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
semshade(A,0.3,'k',[1:N])
hold on
semshade(B,0.3,'m',[1:N])
%

 % [pAnova,pMC,pLabel] = anovaRM2W_full_BH(A(:,11:20),B(:,11:20),1);
 % [pAnova2,pMC,pLabel] = anovaRM2W_full_BH(A(:,1:10),B(:,1:10),1);
 % [pAnova3,pMC,pLabel] = anovaRM2W_full_BH(A(:,21:30),B(:,21:30),1);


% ylim([0.6 1.3])
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
        plot(n,1.2,'k*')
    end
end
ylim([0 9])

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
semshade(A,0.1,'k',[1:N])
hold on
semshade(B,0.1,'m',[1:N])
 % [pAnova,pMC,pLabel] = anovaRM2W_full_BH(A(:,11:20),B(:,11:20),1);
 % [pAnova2,pMC,pLabel] = anovaRM2W_full_BH(A(:,1:10),B(:,1:10),1);
 % [pAnova3,pMC,pLabel] = anovaRM2W_full_BH(A(:,21:30),B(:,21:30),1);


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
        plot(n,1.2,'k*')
    end
end
ylim([0 9])
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
semshade(A,0.3,'m',[1:N])
hold on
semshade(B,0.1,'m',[1:N])
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
        plot(n,1.2,'k*')
    end
end
ylim([0 9])


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
semshade(A,0.3,'k',[1:N])
hold on
semshade(B,0.1,'k',[1:N])
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
        plot(n,1.2,'k*')
    end
end
ylim([0 9])
saveas(gcf,'RBRTogetherAbsBeforeReward.fig')

%% shading: each session is normalized by its total numbers of licks (before stim and after)

%this method amplifies the difference if there is only one lick.
%given how variable this data is, it is better to use the orignial number
%because if there is one lick in before session, it masses up the whole
%stats. but having just one lick is not a realible case.
idx=idxCHR2RS;
RBR=[absPreRewardLickBeforeRBR(idx,:) absPreRewardLickStimRBR(idx,:) absPreRewardLickAfterRBR(idx,:)];

idx=idxCHR2CS;
RBRCS=[absPreRewardLickBeforeRBR(idx,:) absPreRewardLickStimRBR(idx,:) absPreRewardLickAfterRBR(idx,:)];

idx=idxGFPCS;
RBRCSGFP=[absPreRewardLickBeforeRBR(idx,:) absPreRewardLickStimRBR(idx,:) absPreRewardLickAfterRBR(idx,:)];

idx=idxGFPRS;
RBRRSGFP=[absPreRewardLickBeforeRBR(idx,:) absPreRewardLickStimRBR(idx,:) absPreRewardLickAfterRBR(idx,:)];

N=30;

figure
subplot(221)
A=RBR(:,1:N);
A=A./nansum(A,2);
% for n=1:size(A,1);
    ma=nanmean(A(:,[1:10]),'all');
%     A(n,:)=A(n,:)/ma;
% end
% A(isinf(A))=nan;
A=A/ma;

B=RBRCS(:,1:N);
B=B./nansum(B,2);
% for n=1:size(B,1);
    ma=nanmean(B(:,[1:10]),'all');
%     B(n,:)=B(n,:)/ma;
% end
% B(isinf(B))=nan;
B=B/ma;
semshade(A,0.3,'k',[1:N])
hold on
semshade(B,0.3,'m',[1:N])
%

 % [pAnova,pMC,pLabel] = anovaRM2W_full_BH(A(:,11:20),B(:,11:20),1);
 % [pAnova2,pMC,pLabel] = anovaRM2W_full_BH(A(:,1:10),B(:,1:10),1);
 % [pAnova3,pMC,pLabel] = anovaRM2W_full_BH(A(:,21:30),B(:,21:30),1);


% ylim([0.6 1.3])
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
        plot(n,1.2,'k*')
    end
end
ylim([0 17])

subplot(222)
A=RBRRSGFP(:,1:N);
A=A./nansum(A,2);
% for n=1:size(A,1);
    ma=nanmean(A(:,[1:10]),'all');
%     A(n,:)=A(n,:)/ma;
% end
A=A/ma;

B=RBRCSGFP(:,1:N);
B=B./nansum(B,2);
% for n=1:size(B,1);
    ma=nanmean(B(:,[1:10]),'all');
%     B(n,:)=B(n,:)/ma;
% end
B=B/ma;
semshade(A,0.1,'k',[1:N])
hold on
semshade(B,0.1,'m',[1:N])
 % [pAnova,pMC,pLabel] = anovaRM2W_full_BH(A(:,11:20),B(:,11:20),1);
 % [pAnova2,pMC,pLabel] = anovaRM2W_full_BH(A(:,1:10),B(:,1:10),1);
 % [pAnova3,pMC,pLabel] = anovaRM2W_full_BH(A(:,21:30),B(:,21:30),1);


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
        plot(n,1.2,'k*')
    end
end
ylim([0 17])
%

subplot(224)
A=RBRCS(:,1:N);
A=A./nansum(A,2);
% for n=1:size(A,1);
    ma=nanmean(A(:,[1:10]),'all');
%     A(n,:)=A(n,:)/ma;
% end

A=A/ma;

B=RBRCSGFP(:,1:N);
B=B./nansum(B,2);
% for n=1:size(B,1);
    ma=nanmean(B(:,[1:10]),'all');
%     B(n,:)=B(n,:)/ma;
% end
B=B/ma;
semshade(A,0.3,'m',[1:N])
hold on
semshade(B,0.1,'m',[1:N])
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
        plot(n,1.2,'k*')
    end
end
ylim([0 17])


subplot(223)
A=RBR(:,1:N);
A=A./nansum(A,2);
% for n=1:size(A,1);
    ma=nanmean(A(:,[1:10]),'all');
%     A(n,:)=A(n,:)/ma;
% end
A=A/ma;
B=RBRRSGFP(:,1:N);
B=B./nansum(B,2);
% for n=1:size(B,1);
    ma=nanmean(B(:,[1:10]),'all');
%     B(n,:)=B(n,:)/ma;
% end
B=B/ma;
semshade(A,0.3,'k',[1:N])
hold on
semshade(B,0.1,'k',[1:N])
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
        plot(n,1.2,'k*')
    end
end
ylim([0 17])
saveas(gcf,'RBRTogetherAbsBeforeRewardNorm.fig')

%% plot the distribution of licks: JUST CHR2CS, CHR2RS, GFPCS, GFP RS

load('Z:\papers\2025RewardOptogenetics\4mData\250128_data\cueTemp_4m.mat')


idx=idxCHR2CS;
[lickUsePositionStimAll,~,~,~] = absPreRewardLicking_RBR(lickSegAllData(idx),rewardSegAllData(idx),positionSegAllData(idx),threshold(idx),beforeDistance,afterDistance);
[lickUsePositionBeforeAll,~,~,~] = absPreRewardLicking_RBR(lickSegAllDataBefore(idx),rewardSegAllDataBefore(idx),positionSegAllDataBefore(idx),threshold(idx),beforeDistance,afterDistance);
lickUsePositionStimAll=lickUsePositionStimAll(lickUsePositionStimAll>50&lickUsePositionStimAll<366);
lickUsePositionBeforeAll=lickUsePositionBeforeAll(lickUsePositionBeforeAll>50&lickUsePositionBeforeAll<366);

h=histogram(lickUsePositionBeforeAll,50);
x1=h.BinEdges(1:end-1)+h.BinWidth/2;
y1=h.Values/length(h.Data);

h=histogram(lickUsePositionStimAll,50);
x2=h.BinEdges(1:end-1)+h.BinWidth/2;
y2=h.Values/length(h.Data);



idx=idxCHR2RS;
[lickUsePositionStimAll,~,~,~] = absPreRewardLicking_RBR(lickSegAllData(idx),rewardSegAllData(idx),positionSegAllData(idx),threshold(idx),beforeDistance,afterDistance);
[lickUsePositionBeforeAll,~,~,~] = absPreRewardLicking_RBR(lickSegAllDataBefore(idx),rewardSegAllDataBefore(idx),positionSegAllDataBefore(idx),threshold(idx),beforeDistance,afterDistance);
lickUsePositionStimAll=lickUsePositionStimAll(lickUsePositionStimAll>50&lickUsePositionStimAll<366);
lickUsePositionBeforeAll=lickUsePositionBeforeAll(lickUsePositionBeforeAll>50&lickUsePositionBeforeAll<366);
h=histogram(lickUsePositionBeforeAll,50);
x3=h.BinEdges(1:end-1)+h.BinWidth/2;
y3=h.Values/length(h.Data);

h=histogram(lickUsePositionStimAll,50);
x4=h.BinEdges(1:end-1)+h.BinWidth/2;
y4=h.Values/length(h.Data);


idx=idxGFPCS;
[lickUsePositionStimAll,~,~,~] = absPreRewardLicking_RBR(lickSegAllData(idx),rewardSegAllData(idx),positionSegAllData(idx),threshold(idx),beforeDistance,afterDistance);
[lickUsePositionBeforeAll,~,~,~] = absPreRewardLicking_RBR(lickSegAllDataBefore(idx),rewardSegAllDataBefore(idx),positionSegAllDataBefore(idx),threshold(idx),beforeDistance,afterDistance);
lickUsePositionStimAll=lickUsePositionStimAll(lickUsePositionStimAll>50&lickUsePositionStimAll<366);
lickUsePositionBeforeAll=lickUsePositionBeforeAll(lickUsePositionBeforeAll>50&lickUsePositionBeforeAll<366);

h=histogram(lickUsePositionBeforeAll,50);
x5=h.BinEdges(1:end-1)+h.BinWidth/2;
y5=h.Values/length(h.Data);

h=histogram(lickUsePositionStimAll,50);
x6=h.BinEdges(1:end-1)+h.BinWidth/2;
y6=h.Values/length(h.Data);


idx=idxGFPRS;
[lickUsePositionStimAll,~,~,~] = absPreRewardLicking_RBR(lickSegAllData(idx),rewardSegAllData(idx),positionSegAllData(idx),threshold(idx),beforeDistance,afterDistance);
[lickUsePositionBeforeAll,~,~,~] = absPreRewardLicking_RBR(lickSegAllDataBefore(idx),rewardSegAllDataBefore(idx),positionSegAllDataBefore(idx),threshold(idx),beforeDistance,afterDistance);
lickUsePositionStimAll=lickUsePositionStimAll(lickUsePositionStimAll>50&lickUsePositionStimAll<366);
lickUsePositionBeforeAll=lickUsePositionBeforeAll(lickUsePositionBeforeAll>50&lickUsePositionBeforeAll<366);


h=histogram(lickUsePositionBeforeAll,50);
x7=h.BinEdges(1:end-1)+h.BinWidth/2;
y7=h.Values/length(h.Data);

h=histogram(lickUsePositionStimAll,50);
x8=h.BinEdges(1:end-1)+h.BinWidth/2;
y8=h.Values/length(h.Data);

close all

figure
subplot(221)
bar(x1,y1)
hold on
bar(x2,y2)
hold on
line([366 366],[0 0.2])
hold on
plot(resample(cueTemp,400,80)*0.2)
title('CHR2CS')

subplot(222)
bar(x3,y3)
hold on
bar(x4,y4)
hold on
line([366 366],[0 0.2])
hold on
plot(resample(cueTemp,400,80)*0.2)
title('CHR2RS')

subplot(223)
bar(x5,y5)
hold on
bar(x6,y6)
hold on
line([366 366],[0 0.2])
hold on
plot(resample(cueTemp,400,80)*0.2)
title('GFPCS')

subplot(224)
bar(x7,y7)
hold on
bar(x8,y8)
hold on
line([366 366],[0 0.2])
hold on
plot(resample(cueTemp,400,80)*0.2)

title('GFPRS')

saveas(gcf,'lickDistributionALL.fig')


%% plot the distribution of licks: JUST CHR2CS, CHR2RS, GFPCS, GFP RS_NO NORMALIZE

load('Z:\papers\2025RewardOptogenetics\4mData\250128_data\cueTemp_4m.mat')
figure

idx=idxCHR2CS;
[lickUsePositionStimAll,~,~,~] = absPreRewardLicking_RBR(lickSegAllData(idx),rewardSegAllData(idx),positionSegAllData(idx),threshold(idx),beforeDistance,afterDistance);
[lickUsePositionBeforeAll,~,~,~] = absPreRewardLicking_RBR(lickSegAllDataBefore(idx),rewardSegAllDataBefore(idx),positionSegAllDataBefore(idx),threshold(idx),beforeDistance,afterDistance);
lickUsePositionStimAll=lickUsePositionStimAll(lickUsePositionStimAll>50&lickUsePositionStimAll<366);
lickUsePositionBeforeAll=lickUsePositionBeforeAll(lickUsePositionBeforeAll>50&lickUsePositionBeforeAll<366);

subplot(221)
histogram(lickUsePositionBeforeAll,50);
% x1=h.BinEdges(1:end-1)+h.BinWidth/2;
% y1=h.Values/length(h.Data);
hold on
histogram(lickUsePositionStimAll,50);
% x2=h.BinEdges(1:end-1)+h.BinWidth/2;
% y2=h.Values/length(h.Data);

title('CHR2CS')

idx=idxCHR2RS;
[lickUsePositionStimAll,~,~,~] = absPreRewardLicking_RBR(lickSegAllData(idx),rewardSegAllData(idx),positionSegAllData(idx),threshold(idx),beforeDistance,afterDistance);
[lickUsePositionBeforeAll,~,~,~] = absPreRewardLicking_RBR(lickSegAllDataBefore(idx),rewardSegAllDataBefore(idx),positionSegAllDataBefore(idx),threshold(idx),beforeDistance,afterDistance);
lickUsePositionStimAll=lickUsePositionStimAll(lickUsePositionStimAll>50&lickUsePositionStimAll<366);
lickUsePositionBeforeAll=lickUsePositionBeforeAll(lickUsePositionBeforeAll>50&lickUsePositionBeforeAll<366);
subplot(222)
histogram(lickUsePositionBeforeAll,50);
% x1=h.BinEdges(1:end-1)+h.BinWidth/2;
% y1=h.Values/length(h.Data);
hold on
histogram(lickUsePositionStimAll,50);
% x2=h.BinEdges(1:end-1)+h.BinWidth/2;
% y2=h.Values/length(h.Data);

title('CHR2RS')


idx=idxGFPCS;
[lickUsePositionStimAll,~,~,~] = absPreRewardLicking_RBR(lickSegAllData(idx),rewardSegAllData(idx),positionSegAllData(idx),threshold(idx),beforeDistance,afterDistance);
[lickUsePositionBeforeAll,~,~,~] = absPreRewardLicking_RBR(lickSegAllDataBefore(idx),rewardSegAllDataBefore(idx),positionSegAllDataBefore(idx),threshold(idx),beforeDistance,afterDistance);
lickUsePositionStimAll=lickUsePositionStimAll(lickUsePositionStimAll>50&lickUsePositionStimAll<366);
lickUsePositionBeforeAll=lickUsePositionBeforeAll(lickUsePositionBeforeAll>50&lickUsePositionBeforeAll<366);

subplot(223)
histogram(lickUsePositionBeforeAll,50);
% x1=h.BinEdges(1:end-1)+h.BinWidth/2;
% y1=h.Values/length(h.Data);
hold on
histogram(lickUsePositionStimAll,50);
% x2=h.BinEdges(1:end-1)+h.BinWidth/2;
% y2=h.Values/length(h.Data);

title('GFPCS')


idx=idxGFPRS;
[lickUsePositionStimAll,~,~,~] = absPreRewardLicking_RBR(lickSegAllData(idx),rewardSegAllData(idx),positionSegAllData(idx),threshold(idx),beforeDistance,afterDistance);
[lickUsePositionBeforeAll,~,~,~] = absPreRewardLicking_RBR(lickSegAllDataBefore(idx),rewardSegAllDataBefore(idx),positionSegAllDataBefore(idx),threshold(idx),beforeDistance,afterDistance);
lickUsePositionStimAll=lickUsePositionStimAll(lickUsePositionStimAll>50&lickUsePositionStimAll<366);
lickUsePositionBeforeAll=lickUsePositionBeforeAll(lickUsePositionBeforeAll>50&lickUsePositionBeforeAll<366);

subplot(224)
histogram(lickUsePositionBeforeAll,50);
% x1=h.BinEdges(1:end-1)+h.BinWidth/2;
% y1=h.Values/length(h.Data);
hold on
histogram(lickUsePositionStimAll,50);
% x2=h.BinEdges(1:end-1)+h.BinWidth/2;
% y2=h.Values/length(h.Data);

title('GFPRS')



saveas(gcf,'lickDistributionALL_noNormalize.fig')

%% calculate distance: must use the same bin edges
edges=linspace(0,366,101);
width=nanmean(diff(edges));
binWindow=7;

load('Z:\papers\2025RewardOptogenetics\4mData\250128_data\cueTemp_4m.mat')

figure
idx=idxCHR2CS;
[lickUsePositionStimAll,~,~,~] = absPreRewardLicking_RBR(lickSegAllData(idx),rewardSegAllData(idx),positionSegAllData(idx),threshold(idx),beforeDistance,afterDistance);
[lickUsePositionBeforeAll,~,~,~] = absPreRewardLicking_RBR(lickSegAllDataBefore(idx),rewardSegAllDataBefore(idx),positionSegAllDataBefore(idx),threshold(idx),beforeDistance,afterDistance);
% lickUsePositionStimAll=lickUsePositionStimAll(lickUsePositionStimAll>50&lickUsePositionStimAll<366);
% lickUsePositionBeforeAll=lickUsePositionBeforeAll(lickUsePositionBeforeAll>50&lickUsePositionBeforeAll<366);
lickUsePositionStimAll=lickUsePositionStimAll(lickUsePositionStimAll<366);
lickUsePositionBeforeAll=lickUsePositionBeforeAll(lickUsePositionBeforeAll<366);

subplot(421)
h=histogram(lickUsePositionBeforeAll,edges);
countBefore=h.Values;
hold on
h=histogram(lickUsePositionStimAll,edges);
countStim=h.Values;

title('CHR2CS')

subplot(422)
diffCount=countStim-countBefore;

bar(edges(1:end-1)+width/2,diffCount)

avgDiffCount=[];
for n=1:length(diffCount)-binWindow+1;
    avgDiffCount(n)=nanmean(diffCount(n:n+binWindow-1));
end

p=length(find(avgDiffCount(1:end-1)<avgDiffCount(end)))/length(avgDiffCount);

title(['CHR2CS diff',num2str(p)])

idx=idxCHR2RS;
[lickUsePositionStimAll,~,~,~] = absPreRewardLicking_RBR(lickSegAllData(idx),rewardSegAllData(idx),positionSegAllData(idx),threshold(idx),beforeDistance,afterDistance);
[lickUsePositionBeforeAll,~,~,~] = absPreRewardLicking_RBR(lickSegAllDataBefore(idx),rewardSegAllDataBefore(idx),positionSegAllDataBefore(idx),threshold(idx),beforeDistance,afterDistance);
% lickUsePositionStimAll=lickUsePositionStimAll(lickUsePositionStimAll>50&lickUsePositionStimAll<366);
% lickUsePositionBeforeAll=lickUsePositionBeforeAll(lickUsePositionBeforeAll>50&lickUsePositionBeforeAll<366);
lickUsePositionStimAll=lickUsePositionStimAll(lickUsePositionStimAll<366);
lickUsePositionBeforeAll=lickUsePositionBeforeAll(lickUsePositionBeforeAll<366);

subplot(423)
h=histogram(lickUsePositionBeforeAll,edges);
countBefore=h.Values;
hold on
h=histogram(lickUsePositionStimAll,edges);
countStim=h.Values;

title('CHR2RS')

subplot(424)
diffCount=countStim-countBefore;
bar(edges(1:end-1)+width/2,diffCount)

avgDiffCount=[];
for n=1:length(diffCount)-binWindow+1;
    avgDiffCount(n)=nanmean(diffCount(n:n+binWindow-1));
end

p=length(find(avgDiffCount(1:end-1)<avgDiffCount(end)))/length(avgDiffCount);

title(['CHR2RS diff',num2str(p)])

idx=idxGFPCS;
[lickUsePositionStimAll,~,~,~] = absPreRewardLicking_RBR(lickSegAllData(idx),rewardSegAllData(idx),positionSegAllData(idx),threshold(idx),beforeDistance,afterDistance);
[lickUsePositionBeforeAll,~,~,~] = absPreRewardLicking_RBR(lickSegAllDataBefore(idx),rewardSegAllDataBefore(idx),positionSegAllDataBefore(idx),threshold(idx),beforeDistance,afterDistance);
% lickUsePositionStimAll=lickUsePositionStimAll(lickUsePositionStimAll>50&lickUsePositionStimAll<366);
% lickUsePositionBeforeAll=lickUsePositionBeforeAll(lickUsePositionBeforeAll>50&lickUsePositionBeforeAll<366);
lickUsePositionStimAll=lickUsePositionStimAll(lickUsePositionStimAll<366);
lickUsePositionBeforeAll=lickUsePositionBeforeAll(lickUsePositionBeforeAll<366);

subplot(425)
h=histogram(lickUsePositionBeforeAll,edges);
countBefore=h.Values;
hold on
h=histogram(lickUsePositionStimAll,edges);
countStim=h.Values;

title('GFPCS')

subplot(426)
diffCount=countStim-countBefore;
bar(edges(1:end-1)+width/2,diffCount)

avgDiffCount=[];
for n=1:length(diffCount)-binWindow+1;
    avgDiffCount(n)=nanmean(diffCount(n:n+binWindow-1));
end

p=length(find(avgDiffCount(1:end-1)<avgDiffCount(end)))/length(avgDiffCount);

title(['GFPCS diff',num2str(p)])

idx=idxGFPRS;
[lickUsePositionStimAll,~,~,~] = absPreRewardLicking_RBR(lickSegAllData(idx),rewardSegAllData(idx),positionSegAllData(idx),threshold(idx),beforeDistance,afterDistance);
[lickUsePositionBeforeAll,~,~,~] = absPreRewardLicking_RBR(lickSegAllDataBefore(idx),rewardSegAllDataBefore(idx),positionSegAllDataBefore(idx),threshold(idx),beforeDistance,afterDistance);
% lickUsePositionStimAll=lickUsePositionStimAll(lickUsePositionStimAll>50&lickUsePositionStimAll<366);
% lickUsePositionBeforeAll=lickUsePositionBeforeAll(lickUsePositionBeforeAll>50&lickUsePositionBeforeAll<366);
lickUsePositionStimAll=lickUsePositionStimAll(lickUsePositionStimAll<366);
lickUsePositionBeforeAll=lickUsePositionBeforeAll(lickUsePositionBeforeAll<366);

subplot(427)
h=histogram(lickUsePositionBeforeAll,edges);
countBefore=h.Values;
hold on
h=histogram(lickUsePositionStimAll,edges);
countStim=h.Values;

title('GFPRS')

subplot(428)
diffCount=countStim-countBefore;
bar(edges(1:end-1)+width/2,diffCount)

avgDiffCount=[];
for n=1:length(diffCount)-binWindow+1;
    avgDiffCount(n)=nanmean(diffCount(n:n+binWindow-1));
end

p=length(find(avgDiffCount(1:end-1)<avgDiffCount(end)))/length(avgDiffCount);

title(['GFPRS diff',num2str(p)])


saveas(gcf,'lickDistributionALL_noNormalizeChange.fig')


%% %per session: calculate distance: must use the same bin edges
edges=linspace(0,366,81);
width=nanmean(diff(edges));
binWindow=7;
cues=[63 85;163 185;279 301;374 396];
stimLoc=[];
stimLoc(1,:)=[63+3.5 63+3.5+15];
stimLoc(2,:)=[163+3.5 163+3.5+15];
stimLoc(3,:)=[279+3.5 279+3.5+15];
% stimLoc(4,:)=[366-15 366];
binWidth=1;
rewardLoc=366/binWidth;
aroundCueBin=ceil(cues(1:3,:)/width);
adjBin=3;
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


figure
subplot(431)
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[0 50],'Color','r')

bar(edges(1:end-1)+width/2,nanmean(countBeforeAllCHR2CS,1),'FaceColor',[0.5 0.5 0.5])
hold on
bar(edges(1:end-1)+width/2,nanmean(countStimAllCHR2CS,1),'FaceColor','b')
xlim([0 400])


subplot(432)

for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-30 50],'Color','r')

diffCount=countStimAllCHR2CS-countBeforeAllCHR2CS;
bar(edges(1:end-1)+width/2,nanmean(diffCount,1),'FaceColor','m')
xlim([0 400])
ylim([-30 100])
avgDiffCount=[];
for n=1:size(diffCount,2)-binWindow+1;
    avgDiffCount(n)=nanmean(diffCount(:,n:n+binWindow-1),'all');
end

p=length(find(avgDiffCount(1:end-1)<avgDiffCount(end)))/length(avgDiffCount);

title(['CHR2CS diff',num2str(p)])
%
subplot(433)
diffCue=[];
for n=1:size(aroundCueBin,1);
    c=diffCount(:,[aroundCueBin(n,1):aroundCueBin(n,2)]);
    diffCue=[diffCue;c];
end

bar([1:1:size(diffCue,2)],nanmean(diffCue,1))
hold on
line([adjBin adjBin]+0.5,[0 10])

%

subplot(434)
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
%  for i=1:size(stimLoc,1)
% hold on
%         patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
%  end
hold on
line([rewardLoc rewardLoc],[0 50],'Color','r')

bar(edges(1:end-1)+width/2,nanmean(countBeforeAllCHR2RS,1),'FaceColor',[0.5 0.5 0.5])
hold on
bar(edges(1:end-1)+width/2,nanmean(countStimAllCHR2RS,1),'FaceColor','b')
xlim([0 400])
%
subplot(435)

for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
%  for i=1:size(stimLoc,1)
% hold on
%         patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
%  end
hold on
line([rewardLoc rewardLoc],[-50 50],'Color','r')
diffCount=countStimAllCHR2RS-countBeforeAllCHR2RS;
bar(edges(1:end-1)+width/2,nanmean(diffCount,1),'FaceColor',[0.5 0.5 0.5])
xlim([0 400])
avgDiffCount=[];
for n=1:size(diffCount,2)-binWindow+1;
    avgDiffCount(n)=nanmean(diffCount(:,n:n+binWindow-1),'all');
end

p=length(find(avgDiffCount(1:end-1)<avgDiffCount(end)))/length(avgDiffCount);

title(['CHR2RS diff',num2str(p)])

subplot(436)
diffCue=[];
for n=1:size(aroundCueBin,1);
    c=diffCount(:,[aroundCueBin(n,1):aroundCueBin(n,2)]);
    diffCue=[diffCue;c];
end

bar([1:1:size(diffCue,2)],nanmean(diffCue,1))
hold on
line([adjBin adjBin]+0.5,[0 10])
%

subplot(437)
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[0 50],'Color','r')
bar(edges(1:end-1)+width/2,nanmean(countBeforeAllGFPCS,1),'FaceColor',[0.5 0.5 0.5])
hold on
bar(edges(1:end-1)+width/2,nanmean(countStimAllGFPCS,1),'FaceColor','b')
xlim([0 400])
subplot(438)
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 40 40 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 40 40 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-50 40],'Color','r')
diffCount=countStimAllGFPCS-countBeforeAllGFPCS;
bar(edges(1:end-1)+width/2,nanmean(diffCount,1),'FaceColor','m')
xlim([0 400])
avgDiffCount=[];
for n=1:size(diffCount,2)-binWindow+1;
    avgDiffCount(n)=nanmean(diffCount(:,n:n+binWindow-1),'all');
end

p=length(find(avgDiffCount(1:end-1)<avgDiffCount(end)))/length(avgDiffCount);

title(['GFPCS diff',num2str(p)])

subplot(439)
diffCue=[];
for n=1:size(aroundCueBin,1);
    c=diffCount(:,[aroundCueBin(n,1):aroundCueBin(n,2)]);
    diffCue=[diffCue;c];
end

bar([1:1:size(diffCue,2)],nanmean(diffCue,1))
hold on
line([adjBin adjBin]+0.5,[0 10])

subplot(4,3,10)
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
%  for i=1:size(stimLoc,1)
% hold on
%         patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
%  end
hold on
line([rewardLoc rewardLoc],[0 50],'Color','r')
bar(edges(1:end-1)+width/2,nanmean(countBeforeAllGFPRS,1),'FaceColor',[0.5 0.5 0.5])
hold on
bar(edges(1:end-1)+width/2,nanmean(countStimAllGFPRS,1),'FaceColor','b')
xlim([0 400])




subplot(4,3,11)
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
%  for i=1:size(stimLoc,1)
% hold on
%         patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
%  end
hold on
line([rewardLoc rewardLoc],[-50 50],'Color','r')
diffCount=countStimAllGFPRS-countBeforeAllGFPRS;
bar(edges(1:end-1)+width/2,nanmean(diffCount,1),'FaceColor',[0.5 0.5 0.5])
xlim([0 400])
avgDiffCount=[];
for n=1:size(diffCount,2)-binWindow+1;
    avgDiffCount(n)=nanmean(diffCount(:,n:n+binWindow-1),'all');
end
ylim([-20 50])
p=length(find(avgDiffCount(1:end-1)<avgDiffCount(end)))/length(avgDiffCount);

title(['GFPRS diff',num2str(p)])

subplot(4,3,12)
diffCue=[];
for n=1:size(aroundCueBin,1);
    c=diffCount(:,[aroundCueBin(n,1):aroundCueBin(n,2)]);
    diffCue=[diffCue;c];
end

bar([1:1:size(diffCue,2)],nanmean(diffCue,1))
hold on
line([adjBin adjBin]+0.5,[0 10])

saveas(gcf,'lickDistributionALL_noNormalizeChangeAvgSession.fig')

%% %per session: Normalize: each session number was normzed by the total number of licks before and stim.calculate distance: must use the same bin edges
edges=linspace(0,366,81);
width=nanmean(diff(edges));
binWindow=7;
cues=[63 85;163 185;279 301;374 396];
stimLoc=[];
stimLoc(1,:)=[63+3.5 63+3.5+15];
stimLoc(2,:)=[163+3.5 163+3.5+15];
stimLoc(3,:)=[279+3.5 279+3.5+15];
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


figure,
subplot(421)
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 0.01 0.01 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 0.01 0.01 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-0.01 0.01],'Color','r')
M=nansum([countBeforeAllCHR2CS countStimAllCHR2CS],2);
bar(edges(1:end-1)+width/2,nanmean(countBeforeAllCHR2CS./M,1),'FaceColor',[0.5 0.5 0.5])
hold on
bar(edges(1:end-1)+width/2,nanmean(countStimAllCHR2CS./M,1),'FaceColor','b')
xlim([0 400])


subplot(422)

for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 0.01 0.01 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 0.01 0.01 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-0.01 0.01],'Color','r')

diffCount=countStimAllCHR2CS./M-countBeforeAllCHR2CS./M;
bar(edges(1:end-1)+width/2,nanmean(diffCount,1),'FaceColor','m')
xlim([0 400])
% ylim([-30 100])
avgDiffCount=[];
for n=1:size(diffCount,2)-binWindow+1;
    avgDiffCount(n)=nanmean(diffCount(:,n:n+binWindow-1),'all');
end

p=length(find(avgDiffCount(1:end-1)<avgDiffCount(end)))/length(avgDiffCount);

title(['CHR2CS diff',num2str(p)])


subplot(423)
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 0.01 0.01 0],'y')
end
%  for i=1:size(stimLoc,1)
% hold on
%         patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
%  end
hold on
line([rewardLoc rewardLoc],[-0.01 0.01],'Color','r')
M=nansum([countBeforeAllCHR2RS countStimAllCHR2RS],2);
bar(edges(1:end-1)+width/2,nanmean(countBeforeAllCHR2RS./M,1),'FaceColor',[0.5 0.5 0.5])
hold on
bar(edges(1:end-1)+width/2,nanmean(countStimAllCHR2RS./M,1),'FaceColor','b')
xlim([0 400])
subplot(424)

for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 0.01 0.01 0],'y')
end
%  for i=1:size(stimLoc,1)
% hold on
%         patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
%  end
hold on
line([rewardLoc rewardLoc],[-0.01 0.01],'Color','r')
diffCount=countStimAllCHR2RS./M-countBeforeAllCHR2RS./M;
bar(edges(1:end-1)+width/2,nanmean(diffCount,1),'FaceColor',[0.5 0.5 0.5])
xlim([0 400])
avgDiffCount=[];
for n=1:size(diffCount,2)-binWindow+1;
    avgDiffCount(n)=nanmean(diffCount(:,n:n+binWindow-1),'all');
end

p=length(find(avgDiffCount(1:end-1)<avgDiffCount(end)))/length(avgDiffCount);

title(['CHR2RS diff',num2str(p)])




subplot(425)
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 0.01 0.01 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 0.01 0.01 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-0.01 0.01],'Color','r')
M=nansum([countBeforeAllGFPCS countStimAllGFPCS],2);
bar(edges(1:end-1)+width/2,nanmean(countBeforeAllGFPCS./M,1),'FaceColor',[0.5 0.5 0.5])
hold on
bar(edges(1:end-1)+width/2,nanmean(countStimAllGFPCS./M,1),'FaceColor','b')
xlim([0 400])
subplot(426)
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 0.01 0.01 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 0.01 0.01 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-0.01 0.01],'Color','r')
diffCount=countStimAllGFPCS./M-countBeforeAllGFPCS./M;
bar(edges(1:end-1)+width/2,nanmean(diffCount,1),'FaceColor','m')
xlim([0 400])
avgDiffCount=[];
for n=1:size(diffCount,2)-binWindow+1;
    avgDiffCount(n)=nanmean(diffCount(:,n:n+binWindow-1),'all');
end

p=length(find(avgDiffCount(1:end-1)<avgDiffCount(end)))/length(avgDiffCount);

title(['GFPCS diff',num2str(p)])



subplot(427)
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 0.01 0.01 0],'y')
end
%  for i=1:size(stimLoc,1)
% hold on
%         patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
%  end
hold on
line([rewardLoc rewardLoc],[-0.01 0.01],'Color','r')
M=nansum([countBeforeAllGFPRS countStimAllGFPRS],2);
bar(edges(1:end-1)+width/2,nanmean(countBeforeAllGFPRS./M,1),'FaceColor',[0.5 0.5 0.5])
hold on
bar(edges(1:end-1)+width/2,nanmean(countStimAllGFPRS./M,1),'FaceColor','b')
xlim([0 400])

subplot(428)
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 0.01 0.01 0],'y')
end
%  for i=1:size(stimLoc,1)
% hold on
%         patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
%  end
hold on
line([rewardLoc rewardLoc],[-0.01 0.01],'Color','r')
diffCount=countStimAllGFPRS./M-countBeforeAllGFPRS./M;
bar(edges(1:end-1)+width/2,nanmean(diffCount,1),'FaceColor',[0.5 0.5 0.5])
xlim([0 400])
avgDiffCount=[];
for n=1:size(diffCount,2)-binWindow+1;
    avgDiffCount(n)=nanmean(diffCount(:,n:n+binWindow-1),'all');
end
% ylim([-20 50])
p=length(find(avgDiffCount(1:end-1)<avgDiffCount(end)))/length(avgDiffCount);

title(['GFPRS diff',num2str(p)])

saveas(gcf,'lickDistributionALL_NormalizeChangeAvgSession.fig')

%% RBR licks during stimulation
idx=idxCHR2CS;


[lickUsePositionStimRBR,~,~,~] = absPreRewardLicking_RBRRBRLick(lickSegAllData(idx),rewardSegAllData(idx),positionSegAllData(idx),threshold(idx),beforeDistance,afterDistance);
for session=1:length(lickUsePositionStimRBR);
    for lap=1:length(lickUsePositionStimRBR{session});
lickUsePositionStimRBR{session}{lap}=lickUsePositionStimRBR{session}{lap}(lickUsePositionStimRBR{session}{lap}<366);
    end
end

% group licks across runs
countRBR=[];
allStimRBR={};
for lap=1:10;
    allStimRBR{lap}=[];
    for session=1:length(lickUsePositionStimRBR);
        allStimRBR{lap}=[allStimRBR{lap};lickUsePositionStimRBR{session}{lap}];
    end
    h=histogram(allStimRBR{lap},edges);
    countRBR(lap,:)=h.Values;
end

figure,
for n=1:size(countRBR,1)
subplot(10,2,1+2*(n-1))
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-30 50],'Color','r')
hold on
bar(edges(1:end-1)+width/2,countRBR(n,:));
xlim([0 400])

end

% comparing the first 5 and last 5
dayThresh=3;
countRBREarly=nanmean(countRBR(1:dayThresh,:));
countRBRLate=nanmean(countRBR(dayThresh+1:end,:));

subplot(10,2,[2 4 6 8 10])
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-30 50],'Color','r')
hold on
bar(edges(1:end-1)+width/2,countRBREarly);
xlim([0 400])

title(['early runs',num2str(1),'to',num2str(dayThresh)])


subplot(10,2,[12 14 16 18 20])
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-30 50],'Color','r')
hold on
bar(edges(1:end-1)+width/2,countRBRLate);
xlim([0 400])

title(['late runs',num2str(dayThresh+1),'to',num2str(size(countRBR,1))])
tightfig
saveas(gcf,'lickDistriRBRCHR2CS.fig')

% RBR licks during stimulation
idx=idxCHR2RS;


[lickUsePositionStimRBR,~,~,~] = absPreRewardLicking_RBRRBRLick(lickSegAllData(idx),rewardSegAllData(idx),positionSegAllData(idx),threshold(idx),beforeDistance,afterDistance);
for session=1:length(lickUsePositionStimRBR);
    for lap=1:length(lickUsePositionStimRBR{session});
lickUsePositionStimRBR{session}{lap}=lickUsePositionStimRBR{session}{lap}(lickUsePositionStimRBR{session}{lap}<366);
    end
end

% group licks across runs
countRBR=[];
allStimRBR={};
for lap=1:10;
    allStimRBR{lap}=[];
    for session=1:length(lickUsePositionStimRBR);
        allStimRBR{lap}=[allStimRBR{lap};lickUsePositionStimRBR{session}{lap}];
    end
    h=histogram(allStimRBR{lap},edges);
    countRBR(lap,:)=h.Values;
end

figure,
for n=1:size(countRBR,1)
subplot(10,2,1+2*(n-1))
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-30 50],'Color','r')
hold on
bar(edges(1:end-1)+width/2,countRBR(n,:));
xlim([0 400])

end

% comparing the first 5 and last 5
dayThresh=3;
countRBREarly=nanmean(countRBR(1:dayThresh,:));
countRBRLate=nanmean(countRBR(dayThresh+1:end,:));

subplot(10,2,[2 4 6 8 10])
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-30 50],'Color','r')
hold on
bar(edges(1:end-1)+width/2,countRBREarly);
xlim([0 400])

title(['early runs',num2str(1),'to',num2str(dayThresh)])


subplot(10,2,[12 14 16 18 20])
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-30 50],'Color','r')
hold on
bar(edges(1:end-1)+width/2,countRBRLate);
xlim([0 400])

title(['late runs',num2str(dayThresh+1),'to',num2str(size(countRBR,1))])
saveas(gcf,'lickDistriRBRCHR2RS.fig')

%% RBR GFP
idx=idxGFPCS;


[lickUsePositionStimRBR,~,~,~] = absPreRewardLicking_RBRRBRLick(lickSegAllData(idx),rewardSegAllData(idx),positionSegAllData(idx),threshold(idx),beforeDistance,afterDistance);
for session=1:length(lickUsePositionStimRBR);
    for lap=1:length(lickUsePositionStimRBR{session});
lickUsePositionStimRBR{session}{lap}=lickUsePositionStimRBR{session}{lap}(lickUsePositionStimRBR{session}{lap}<366);
    end
end

% group licks across runs
countRBR=[];
allStimRBR={};
for lap=1:10;
    allStimRBR{lap}=[];
    for session=1:length(lickUsePositionStimRBR);
        allStimRBR{lap}=[allStimRBR{lap};lickUsePositionStimRBR{session}{lap}];
    end
    h=histogram(allStimRBR{lap},edges);
    countRBR(lap,:)=h.Values;
end

figure,
for n=1:size(countRBR,1)
subplot(10,2,1+2*(n-1))
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-30 50],'Color','r')
hold on
bar(edges(1:end-1)+width/2,countRBR(n,:));
xlim([0 400])

end

% comparing the first 5 and last 5
dayThresh=3;
countRBREarly=nanmean(countRBR(1:dayThresh,:));
countRBRLate=nanmean(countRBR(dayThresh+1:end,:));

subplot(10,2,[2 4 6 8 10])
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-30 50],'Color','r')
hold on
bar(edges(1:end-1)+width/2,countRBREarly);
xlim([0 400])

title(['early runs',num2str(1),'to',num2str(dayThresh)])


subplot(10,2,[12 14 16 18 20])
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-30 50],'Color','r')
hold on
bar(edges(1:end-1)+width/2,countRBRLate);
xlim([0 400])

title(['late runs',num2str(dayThresh+1),'to',num2str(size(countRBR,1))])
tightfig
saveas(gcf,'lickDistriRBRGFPCS.fig')

%
idx=idxGFPRS;


[lickUsePositionStimRBR,~,~,~] = absPreRewardLicking_RBRRBRLick(lickSegAllData(idx),rewardSegAllData(idx),positionSegAllData(idx),threshold(idx),beforeDistance,afterDistance);
for session=1:length(lickUsePositionStimRBR);
    for lap=1:length(lickUsePositionStimRBR{session});
lickUsePositionStimRBR{session}{lap}=lickUsePositionStimRBR{session}{lap}(lickUsePositionStimRBR{session}{lap}<366);
    end
end

% group licks across runs
countRBR=[];
allStimRBR={};
for lap=1:10;
    allStimRBR{lap}=[];
    for session=1:length(lickUsePositionStimRBR);
        allStimRBR{lap}=[allStimRBR{lap};lickUsePositionStimRBR{session}{lap}];
    end
    h=histogram(allStimRBR{lap},edges);
    countRBR(lap,:)=h.Values;
end

figure,
for n=1:size(countRBR,1)
subplot(10,2,1+2*(n-1))
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-30 50],'Color','r')
hold on
bar(edges(1:end-1)+width/2,countRBR(n,:));
xlim([0 400])

end

% comparing the first 5 and last 5
dayThresh=3;
countRBREarly=nanmean(countRBR(1:dayThresh,:));
countRBRLate=nanmean(countRBR(dayThresh+1:end,:));

subplot(10,2,[2 4 6 8 10])
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-30 50],'Color','r')
hold on
bar(edges(1:end-1)+width/2,countRBREarly);
xlim([0 400])

title(['early runs',num2str(1),'to',num2str(dayThresh)])


subplot(10,2,[12 14 16 18 20])
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-30 50],'Color','r')
hold on
bar(edges(1:end-1)+width/2,countRBRLate);
xlim([0 400])

title(['late runs',num2str(dayThresh+1),'to',num2str(size(countRBR,1))])
saveas(gcf,'lickDistriRBRGFPRS.fig')

%% before stimulation sessions
%% RBR licks during stimulation
idx=idxCHR2CS;


[lickUsePositionBeforeRBR,~,~,~] = absPreRewardLicking_RBRRBRLick(lickSegAllDataBefore(idx),rewardSegAllDataBefore(idx),positionSegAllDataBefore(idx),threshold(idx),beforeDistance,afterDistance);
for session=1:length(lickUsePositionBeforeRBR);
    for lap=1:length(lickUsePositionBeforeRBR{session});
lickUsePositionBeforeRBR{session}{lap}=lickUsePositionBeforeRBR{session}{lap}(lickUsePositionBeforeRBR{session}{lap}<366);
    end
end

% group licks across runs
countRBR=[];
allBeforeRBR={};
for lap=1:10;
    allBeforeRBR{lap}=[];
    for session=1:length(lickUsePositionBeforeRBR);
        allBeforeRBR{lap}=[allBeforeRBR{lap};lickUsePositionBeforeRBR{session}{lap}];
    end
    h=histogram(allBeforeRBR{lap},edges);
    countRBR(lap,:)=h.Values;
end

figure,
for n=1:size(countRBR,1)
subplot(10,2,1+2*(n-1))
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-30 50],'Color','r')
hold on
bar(edges(1:end-1)+width/2,countRBR(n,:));
xlim([0 400])

end

% comparing the first 5 and last 5
dayThresh=3;
countRBREarly=nanmean(countRBR(1:dayThresh,:));
countRBRLate=nanmean(countRBR(dayThresh+1:end,:));

subplot(10,2,[2 4 6 8 10])
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-30 50],'Color','r')
hold on
bar(edges(1:end-1)+width/2,countRBREarly);
xlim([0 400])

title(['early runs',num2str(1),'to',num2str(dayThresh)])


subplot(10,2,[12 14 16 18 20])
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-30 50],'Color','r')
hold on
bar(edges(1:end-1)+width/2,countRBRLate);
xlim([0 400])

title(['late runs',num2str(dayThresh+1),'to',num2str(size(countRBR,1))])
tightfig
saveas(gcf,'lickDistriRBRCHR2CSBefore.fig')

% RBR licks during stimulation
idx=idxCHR2RS;


[lickUsePositionBeforeRBR,~,~,~] = absPreRewardLicking_RBRRBRLick(lickSegAllDataBefore(idx),rewardSegAllDataBefore(idx),positionSegAllDataBefore(idx),threshold(idx),beforeDistance,afterDistance);
for session=1:length(lickUsePositionBeforeRBR);
    for lap=1:length(lickUsePositionBeforeRBR{session});
lickUsePositionBeforeRBR{session}{lap}=lickUsePositionBeforeRBR{session}{lap}(lickUsePositionBeforeRBR{session}{lap}<366);
    end
end

% group licks across runs
countRBR=[];
allBeforeRBR={};
for lap=1:10;
    allBeforeRBR{lap}=[];
    for session=1:length(lickUsePositionBeforeRBR);
        allBeforeRBR{lap}=[allBeforeRBR{lap};lickUsePositionBeforeRBR{session}{lap}];
    end
    h=histogram(allBeforeRBR{lap},edges);
    countRBR(lap,:)=h.Values;
end

figure,
for n=1:size(countRBR,1)
subplot(10,2,1+2*(n-1))
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-30 50],'Color','r')
hold on
bar(edges(1:end-1)+width/2,countRBR(n,:));
xlim([0 400])

end

% comparing the first 5 and last 5
dayThresh=3;
countRBREarly=nanmean(countRBR(1:dayThresh,:));
countRBRLate=nanmean(countRBR(dayThresh+1:end,:));

subplot(10,2,[2 4 6 8 10])
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-30 50],'Color','r')
hold on
bar(edges(1:end-1)+width/2,countRBREarly);
xlim([0 400])

title(['early runs',num2str(1),'to',num2str(dayThresh)])


subplot(10,2,[12 14 16 18 20])
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-30 50],'Color','r')
hold on
bar(edges(1:end-1)+width/2,countRBRLate);
xlim([0 400])

title(['late runs',num2str(dayThresh+1),'to',num2str(size(countRBR,1))])
saveas(gcf,'lickDistriRBRCHR2RSBefore.fig')

%% RBR GFP
idx=idxGFPCS;


[lickUsePositionBeforeRBR,~,~,~] = absPreRewardLicking_RBRRBRLick(lickSegAllDataBefore(idx),rewardSegAllDataBefore(idx),positionSegAllDataBefore(idx),threshold(idx),beforeDistance,afterDistance);
for session=1:length(lickUsePositionBeforeRBR);
    for lap=1:length(lickUsePositionBeforeRBR{session});
lickUsePositionBeforeRBR{session}{lap}=lickUsePositionBeforeRBR{session}{lap}(lickUsePositionBeforeRBR{session}{lap}<366);
    end
end

% group licks across runs
countRBR=[];
allBeforeRBR={};
for lap=1:10;
    allBeforeRBR{lap}=[];
    for session=1:length(lickUsePositionBeforeRBR);
        allBeforeRBR{lap}=[allBeforeRBR{lap};lickUsePositionBeforeRBR{session}{lap}];
    end
    h=histogram(allBeforeRBR{lap},edges);
    countRBR(lap,:)=h.Values;
end

figure,
for n=1:size(countRBR,1)
subplot(10,2,1+2*(n-1))
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-30 50],'Color','r')
hold on
bar(edges(1:end-1)+width/2,countRBR(n,:));
xlim([0 400])

end

% comparing the first 5 and last 5
dayThresh=3;
countRBREarly=nanmean(countRBR(1:dayThresh,:));
countRBRLate=nanmean(countRBR(dayThresh+1:end,:));

subplot(10,2,[2 4 6 8 10])
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-30 50],'Color','r')
hold on
bar(edges(1:end-1)+width/2,countRBREarly);
xlim([0 400])

title(['early runs',num2str(1),'to',num2str(dayThresh)])


subplot(10,2,[12 14 16 18 20])
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-30 50],'Color','r')
hold on
bar(edges(1:end-1)+width/2,countRBRLate);
xlim([0 400])

title(['late runs',num2str(dayThresh+1),'to',num2str(size(countRBR,1))])
tightfig
saveas(gcf,'lickDistriRBRGFPCSBefore.fig')

%
idx=idxGFPRS;


[lickUsePositionBeforeRBR,~,~,~] = absPreRewardLicking_RBRRBRLick(lickSegAllDataBefore(idx),rewardSegAllDataBefore(idx),positionSegAllDataBefore(idx),threshold(idx),beforeDistance,afterDistance);
for session=1:length(lickUsePositionBeforeRBR);
    for lap=1:length(lickUsePositionBeforeRBR{session});
lickUsePositionBeforeRBR{session}{lap}=lickUsePositionBeforeRBR{session}{lap}(lickUsePositionBeforeRBR{session}{lap}<366);
    end
end

% group licks across runs
countRBR=[];
allBeforeRBR={};
for lap=1:10;
    allBeforeRBR{lap}=[];
    for session=1:length(lickUsePositionBeforeRBR);
        allBeforeRBR{lap}=[allBeforeRBR{lap};lickUsePositionBeforeRBR{session}{lap}];
    end
    h=histogram(allBeforeRBR{lap},edges);
    countRBR(lap,:)=h.Values;
end

figure,
for n=1:size(countRBR,1)
subplot(10,2,1+2*(n-1))
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-30 50],'Color','r')
hold on
bar(edges(1:end-1)+width/2,countRBR(n,:));
xlim([0 400])

end

% comparing the first 5 and last 5
dayThresh=3;
countRBREarly=nanmean(countRBR(1:dayThresh,:));
countRBRLate=nanmean(countRBR(dayThresh+1:end,:));

subplot(10,2,[2 4 6 8 10])
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-30 50],'Color','r')
hold on
bar(edges(1:end-1)+width/2,countRBREarly);
xlim([0 400])

title(['early runs',num2str(1),'to',num2str(dayThresh)])


subplot(10,2,[12 14 16 18 20])
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-30 50],'Color','r')
hold on
bar(edges(1:end-1)+width/2,countRBRLate);
xlim([0 400])

title(['late runs',num2str(dayThresh+1),'to',num2str(size(countRBR,1))])
saveas(gcf,'lickDistriRBRGFPRSBefore.fig')

%% after stimulation
%% RBR licks during stimulation
idx=idxCHR2CS;


[lickUsePositionAfterRBR,~,~,~] = absPreRewardLicking_RBRRBRLick(lickSegAllDataAfter(idx),rewardSegAllDataAfter(idx),positionSegAllDataAfter(idx),threshold(idx),beforeDistance,afterDistance);
for session=1:length(lickUsePositionAfterRBR);
    for lap=1:length(lickUsePositionAfterRBR{session});
lickUsePositionAfterRBR{session}{lap}=lickUsePositionAfterRBR{session}{lap}(lickUsePositionAfterRBR{session}{lap}<366);
    end
end

% group licks across runs
countRBR=[];
allAfterRBR={};
for lap=1:10;
    allAfterRBR{lap}=[];
    for session=1:length(lickUsePositionAfterRBR);
        allAfterRBR{lap}=[allAfterRBR{lap};lickUsePositionAfterRBR{session}{lap}];
    end
    h=histogram(allAfterRBR{lap},edges);
    countRBR(lap,:)=h.Values;
end

figure,
for n=1:size(countRBR,1)
subplot(10,2,1+2*(n-1))
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-30 50],'Color','r')
hold on
bar(edges(1:end-1)+width/2,countRBR(n,:));
xlim([0 400])

end

% comparing the first 5 and last 5
dayThresh=3;
countRBREarly=nanmean(countRBR(1:dayThresh,:));
countRBRLate=nanmean(countRBR(dayThresh+1:end,:));

subplot(10,2,[2 4 6 8 10])
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-30 50],'Color','r')
hold on
bar(edges(1:end-1)+width/2,countRBREarly);
xlim([0 400])

title(['early runs',num2str(1),'to',num2str(dayThresh)])


subplot(10,2,[12 14 16 18 20])
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-30 50],'Color','r')
hold on
bar(edges(1:end-1)+width/2,countRBRLate);
xlim([0 400])

title(['late runs',num2str(dayThresh+1),'to',num2str(size(countRBR,1))])
tightfig
saveas(gcf,'lickDistriRBRCHR2CSAfter.fig')

% RBR licks during stimulation
idx=idxCHR2RS;


[lickUsePositionAfterRBR,~,~,~] = absPreRewardLicking_RBRRBRLick(lickSegAllDataAfter(idx),rewardSegAllDataAfter(idx),positionSegAllDataAfter(idx),threshold(idx),beforeDistance,afterDistance);
for session=1:length(lickUsePositionAfterRBR);
    for lap=1:length(lickUsePositionAfterRBR{session});
lickUsePositionAfterRBR{session}{lap}=lickUsePositionAfterRBR{session}{lap}(lickUsePositionAfterRBR{session}{lap}<366);
    end
end

% group licks across runs
countRBR=[];
allAfterRBR={};
for lap=1:10;
    allAfterRBR{lap}=[];
    for session=1:length(lickUsePositionAfterRBR);
        allAfterRBR{lap}=[allAfterRBR{lap};lickUsePositionAfterRBR{session}{lap}];
    end
    h=histogram(allAfterRBR{lap},edges);
    countRBR(lap,:)=h.Values;
end

figure,
for n=1:size(countRBR,1)
subplot(10,2,1+2*(n-1))
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-30 50],'Color','r')
hold on
bar(edges(1:end-1)+width/2,countRBR(n,:));
xlim([0 400])

end

% comparing the first 5 and last 5
dayThresh=3;
countRBREarly=nanmean(countRBR(1:dayThresh,:));
countRBRLate=nanmean(countRBR(dayThresh+1:end,:));

subplot(10,2,[2 4 6 8 10])
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-30 50],'Color','r')
hold on
bar(edges(1:end-1)+width/2,countRBREarly);
xlim([0 400])

title(['early runs',num2str(1),'to',num2str(dayThresh)])


subplot(10,2,[12 14 16 18 20])
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-30 50],'Color','r')
hold on
bar(edges(1:end-1)+width/2,countRBRLate);
xlim([0 400])

title(['late runs',num2str(dayThresh+1),'to',num2str(size(countRBR,1))])
saveas(gcf,'lickDistriRBRCHR2RSAfter.fig')

%% RBR GFP
idx=idxGFPCS;


[lickUsePositionAfterRBR,~,~,~] = absPreRewardLicking_RBRRBRLick(lickSegAllDataAfter(idx),rewardSegAllDataAfter(idx),positionSegAllDataAfter(idx),threshold(idx),beforeDistance,afterDistance);
for session=1:length(lickUsePositionAfterRBR);
    for lap=1:length(lickUsePositionAfterRBR{session});
lickUsePositionAfterRBR{session}{lap}=lickUsePositionAfterRBR{session}{lap}(lickUsePositionAfterRBR{session}{lap}<366);
    end
end

% group licks across runs
countRBR=[];
allAfterRBR={};
for lap=1:10;
    allAfterRBR{lap}=[];
    for session=1:length(lickUsePositionAfterRBR);
        allAfterRBR{lap}=[allAfterRBR{lap};lickUsePositionAfterRBR{session}{lap}];
    end
    h=histogram(allAfterRBR{lap},edges);
    countRBR(lap,:)=h.Values;
end

figure,
for n=1:size(countRBR,1)
subplot(10,2,1+2*(n-1))
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-30 50],'Color','r')
hold on
bar(edges(1:end-1)+width/2,countRBR(n,:));
xlim([0 400])

end

% comparing the first 5 and last 5
dayThresh=3;
countRBREarly=nanmean(countRBR(1:dayThresh,:));
countRBRLate=nanmean(countRBR(dayThresh+1:end,:));

subplot(10,2,[2 4 6 8 10])
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-30 50],'Color','r')
hold on
bar(edges(1:end-1)+width/2,countRBREarly);
xlim([0 400])

title(['early runs',num2str(1),'to',num2str(dayThresh)])


subplot(10,2,[12 14 16 18 20])
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-30 50],'Color','r')
hold on
bar(edges(1:end-1)+width/2,countRBRLate);
xlim([0 400])

title(['late runs',num2str(dayThresh+1),'to',num2str(size(countRBR,1))])
tightfig
saveas(gcf,'lickDistriRBRGFPCSAfter.fig')

%
idx=idxGFPRS;


[lickUsePositionAfterRBR,~,~,~] = absPreRewardLicking_RBRRBRLick(lickSegAllDataAfter(idx),rewardSegAllDataAfter(idx),positionSegAllDataAfter(idx),threshold(idx),beforeDistance,afterDistance);
for session=1:length(lickUsePositionAfterRBR);
    for lap=1:length(lickUsePositionAfterRBR{session});
lickUsePositionAfterRBR{session}{lap}=lickUsePositionAfterRBR{session}{lap}(lickUsePositionAfterRBR{session}{lap}<366);
    end
end

% group licks across runs
countRBR=[];
allAfterRBR={};
for lap=1:10;
    allAfterRBR{lap}=[];
    for session=1:length(lickUsePositionAfterRBR);
        allAfterRBR{lap}=[allAfterRBR{lap};lickUsePositionAfterRBR{session}{lap}];
    end
    h=histogram(allAfterRBR{lap},edges);
    countRBR(lap,:)=h.Values;
end

figure,
for n=1:size(countRBR,1)
subplot(10,2,1+2*(n-1))
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-30 50],'Color','r')
hold on
bar(edges(1:end-1)+width/2,countRBR(n,:));
xlim([0 400])

end

% comparing the first 5 and last 5
dayThresh=3;
countRBREarly=nanmean(countRBR(1:dayThresh,:));
countRBRLate=nanmean(countRBR(dayThresh+1:end,:));

subplot(10,2,[2 4 6 8 10])
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-30 50],'Color','r')
hold on
bar(edges(1:end-1)+width/2,countRBREarly);
xlim([0 400])

title(['early runs',num2str(1),'to',num2str(dayThresh)])


subplot(10,2,[12 14 16 18 20])
for i=1:size(cues,1)
    hold on
    patch([cues(i,1) cues(i,1) cues(i,2) cues(i,2)],[0 50 50 0],'y')
end
 for i=1:size(stimLoc,1)
hold on
        patch([stimLoc(i,1) stimLoc(i,1) stimLoc(i,2) stimLoc(i,2)],[0 50 50 0],'b')
 end
hold on
line([rewardLoc rewardLoc],[-30 50],'Color','r')
hold on
bar(edges(1:end-1)+width/2,countRBRLate);
xlim([0 400])

title(['late runs',num2str(dayThresh+1),'to',num2str(size(countRBR,1))])
saveas(gcf,'lickDistriRBRGFPRSAfter.fig')