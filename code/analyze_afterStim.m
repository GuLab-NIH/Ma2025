
%this is the best. removed data reward stim GFP when mice changed rigs,
%which caused a huge change in behavior.
%organize the excel data here (20250124-After CS or RS.xlsx) baed on all the data used here: 
% Z:\papers\2025RewardOptogenetics\optRBR\reward

%all stim was used,even though some stil was removed from the bar graph
%(comparing pre, during, and post) due tl baseline issue. But this analysis
%is independent of the baseline issue. So we included all data here.
%Included all sessions after a poor performer CS or RS (the after session
%could be anything, like good performer session) as long at they are on the
%next day or the same day.



%15cmBR-ChR2-summary_inYGData: get chr2RewardCS, chr2RewardRS
%15cmBR-GFP-summary_inYGData: get gfpRewardCS, gfpRewardRS
%3x5cm_15cm-ChR2-summaryinYGData: get chr2DualCS, chr2DualRS
%3x5cm_15cm-GFP-summaryinYGData: get gfpDualCS, gfpDualRS
load('chr2RewardCS.mat')
load('chr2RewardRS.mat')
load('gfpRewardCS.mat')
load('gfpRewardRS.mat')

load('chr2DualCS.mat')
load('chr2DualRS.mat')
load('gfpDualCS.mat')
load('gfpDualRS.mat')


p=[];
[~,p(1)]=ttest(chr2RewardCS,0);
[~,p(2)]=ttest(chr2RewardRS,0);
[~,p(3)]=ttest(gfpRewardCS,0);
[~,p(4)]=ttest(gfpRewardRS,0);

[~,p(5)]=ttest(chr2DualCS,0);
[~,p(6)]=ttest(chr2DualRS,0);
[~,p(7)]=ttest(gfpDualCS,0);
[~,p(8)]=ttest(gfpDualRS,0);

stimRewardAll={};
stimRewardAll{1}=chr2RewardCS;
stimRewardAll{2}=chr2RewardRS;
stimRewardAll{3}=gfpRewardCS;
stimRewardAll{4}=gfpRewardRS;

stimDualAll={};
stimDualAll{1}=chr2DualCS;
stimDualAll{2}=chr2DualRS;
stimDualAll{3}=gfpDualCS;
stimDualAll{4}=gfpDualRS;

save('stimRewardAll.mat','stimRewardAll')
save('stimDualAll.mat','stimDualAll')
%% 
r=0.3;
figure
load('stimRewardAll.mat');
a=stimRewardAll;
a=a([1:1:4]);
M=[];
E=[];
for n=1:length(a)
    M(n)=nanmean(a{n},1);
    E(n)=nansem(a{n},1);
end

x=[1 2 4 5];
% bar(x,M);

subplot(311)
bar(x,M,'FaceColor','w')
hold on
e=errorbar(x,M,E,'k.');


for n=1:length(a);
    x1=x(n);
    x1=x1+(rand(length(a{n}),1)-0.5)*r;
    hold on
    if n==1 || n==3;
plot(x1,a{n},'r.')
    else
        plot(x1,a{n},'k.')
    end
end

xlim([0 6])    
ylim([-35 25])

for n=1:length(a)
[~,sig(2,n)]=ttest(a{n},0);
end
pp=[];
[~,pp(1)]=ttest2(a{1},a{2});
[~,pp(2)]=ttest2(a{3},a{4});
[~,pp(3)]=ttest2(a{1},a{3});
[~,pp(4)]=ttest2(a{2},a{4});
title(['15 p=',num2str([p(1:4) pp])])

%%


load('stimDualAll.mat');
a=stimDualAll;
a=a([1:1:4]);
M=[];
E=[];
for n=1:length(a)
    M(n)=nanmean(a{n},1);
    E(n)=nansem(a{n},1);
end

x=[1 2 4 5];
% bar(x,M);

subplot(312)
bar(x,M,'FaceColor','w')
hold on
e=errorbar(x,M,E,'k.');


for n=1:length(a);
    x1=x(n);
    x1=x1+(rand(length(a{n}),1)-0.5)*r;
    hold on
    if n==1 || n==3;
plot(x1,a{n},'r.')
    else
        plot(x1,a{n},'k.')
    end
end

xlim([0 6])    
ylim([-35 25])

for n=1:length(a)
[~,sig(2,n)]=ttest(a{n},0);
end
pp=[];
[~,pp(1)]=ttest2(a{1},a{2});
[~,pp(2)]=ttest2(a{3},a{4});
[~,pp(3)]=ttest2(a{1},a{3});
[~,pp(4)]=ttest2(a{2},a{4});
title(['3x515 p=',num2str([p(5:8) pp])])

%%
load('stimDualAll.mat');
load('stimRewardAll.mat');

a={};
for n=1:4;
    a{n}=[stimDualAll{n};stimRewardAll{n}];
end
a=a([1:1:4]);
M=[];
E=[];
for n=1:length(a)
    M(n)=nanmean(a{n},1);
    E(n)=nansem(a{n},1);
end

x=[1 2 4 5];
% bar(x,M);

subplot(313)
bar(x,M,'FaceColor','w')
hold on
e=errorbar(x,M,E,'k.');


for n=1:length(a);
    x1=x(n);
    x1=x1+(rand(length(a{n}),1)-0.5)*r;
    hold on
    if n==1 || n==3;
plot(x1,a{n},'r.')
    else
        plot(x1,a{n},'k.')
    end
end

xlim([0 6])    
ylim([-35 25])

for n=1:length(a)
[~,sig(2,n)]=ttest(a{n},0);
end

pCombine=[];
[~,pCombine(1)]=ttest(a{1},0);
[~,pCombine(2)]=ttest(a{2},0);
[~,pCombine(3)]=ttest(a{3},0);
[~,pCombine(4)]=ttest(a{4},0);
pp=[];
[~,pp(1)]=ttest2(a{1},a{2});
[~,pp(2)]=ttest2(a{3},a{4});
[~,pp(3)]=ttest2(a{1},a{3});
[~,pp(4)]=ttest2(a{2},a{4});
title(['combine p=',num2str([pCombine pp])])

saveas(gcf,'postStim.fig')
 print -painters -depsc postStim.eps


