load('activity4m.mat')
load('behavior4m.mat')

figure
subplot(1,4,1)
a=behavior4m.slow.good;
m=nanmean(a,1);
e=nansem(a,1);

errorbar([1:1:10],m,e,'g')
[r,p]=corr([1:1:10]',m');
title(['behavior',num2str(p)])

subplot(1,4,2)
b=activity4m.cell.good.preRew;
m=[];
e=[];
for n=1:length(b);
    m(n)=nanmean(b{n});
    e(n)=nansem(b{n},1);
end

errorbar([1:1:10],m,e,'k')
ylim([0.11 0.34])
[r,p]=corr([1:1:10]',m');
title(['pre reward',num2str(p)])

subplot(1,4,3)
b=activity4m.cell.good.inCue;
m=[];
e=[];
for n=1:length(b);
    m(n)=nanmean(b{n});
    e(n)=nansem(b{n},1);
end

errorbar([1:1:10],m,e,'k')
ylim([0.11 0.34])
[r,p]=corr([1:1:10]',m');
title(['in cue',num2str(p)])

subplot(1,4,4)
b=activity4m.cell.good.outCue;
m=[];
e=[];
for n=1:length(b);
    m(n)=nanmean(b{n});
    e(n)=nansem(b{n},1);
end

errorbar([1:1:10],m,e,'k')
ylim([0.11 0.34])
[r,p]=corr([1:1:10]',m');
title(['out cue',num2str(p)])

saveas(gcf,'behaviorActivity.fig')
print -painters -depsc behaviorActivity.eps

%% plotting consistency and PS together
load('activity4m.mat')
load('behavior4m.mat')
figure
a=behavior4m.slow.good;
m=nanmean(a,1);
e=nansem(a,1);
subplot(131)
b=activity4m.cell.good.preRew;
m1=[];
e1=[];
for n=1:length(b);
    m1(n)=nanmean(b{n});
    e1(n)=nansem(b{n},1);
end
yyaxis left
errorbar([1:1:10],m1,e1,'k')
ylim([0.11 0.37])
hold on
yyaxis right
errorbar([1:1:10],m,e,'g')
ylim([50 85])
xlim([0 11])
title('reward')

subplot(132)
b=activity4m.cell.good.inCue;
m1=[];
e1=[];
for n=1:length(b);
    m1(n)=nanmean(b{n});
    e1(n)=nansem(b{n},1);
end
yyaxis left
errorbar([1:1:10],m1,e1,'k')
ylim([0.11 0.37])
hold on
yyaxis right
errorbar([1:1:10],m,e,'g')
ylim([53 105])
xlim([0 11])
title('cue')


subplot(133)
b=activity4m.cell.good.outCue;
m1=[];
e1=[];
for n=1:length(b);
    m1(n)=nanmean(b{n});
    e1(n)=nansem(b{n},1);
end
yyaxis left
errorbar([1:1:10],m1,e1,'k')
ylim([0.11 0.37])
hold on
yyaxis right
errorbar([1:1:10],m,e,'g')
ylim([55 115])
xlim([0 11])
title('out cue')

saveas(gcf,'behaviorActivityPlotTogether.fig')
print -painters -depsc behaviorActivityPlotTogether.eps

%% cdf before and after learning
beforeDay=[1 2];
afterDay=[7 8 9 10];

figure
subplot(1,4,1)
a=behavior4m.slow.good;
before=a(:,beforeDay);
before=reshape(before,[size(before,1)*size(before,2),1]);
after=a(:,afterDay);
after=reshape(after,[size(after,1)*size(after,2),1]);

cdfplot(before)
hold on
cdfplot(after)
[~,p]=kstest2(before, after);
title(['behavior',num2str(p)])

subplot(1,4,2)
a=activity4m.cell.good.preRew;
for n=1:length(a)
    a{n}=a{n}';
end
before=cell2mat(a(beforeDay));

after=cell2mat(a(afterDay));

cdfplot(before)
hold on
cdfplot(after)
[~,p]=kstest2(before, after);
title(['pre reward',num2str(p)])

subplot(1,4,3)
a=activity4m.cell.good.inCue;
for n=1:length(a)
    a{n}=a{n}';
end
before=cell2mat(a(beforeDay));

after=cell2mat(a(afterDay));

cdfplot(before)
hold on
cdfplot(after)
[~,p]=kstest2(before, after);
title(['in cue',num2str(p)])

subplot(1,4,4)
a=activity4m.cell.good.outCue;
for n=1:length(a)
    a{n}=a{n}';
end
before=cell2mat(a(beforeDay));

after=cell2mat(a(afterDay));

cdfplot(before)
hold on
cdfplot(after)
[~,p]=kstest2(before, after);
title(['out cue',num2str(p)])
saveas(gcf,'behaviorActivityCDF.fig')
print -painters -depsc behaviorActivityCDF.eps
