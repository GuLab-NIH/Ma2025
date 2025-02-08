%this is based on this code: 
% Z:\papers\2025RewardOptogenetics\optRBR\reward\testChR2CS\removingSessionsBasedOnPreCondition_fixed30
load('RBRCS.mat')
load('RBR.mat');
load('RBRCSGFP.mat')
load('RBRRSGFP.mat');
all={};
all{1}=RBRCS;
all{2}=RBR;
all{3}=RBRCSGFP;
all{4}=RBRRSGFP;
allIdxAll={};

figure
for t=1:length(all);
    this=all{t};
    prePoor=this(:,1:10);
    durPoor=this(:,11:20);
    postPoor=this(:,21:30);

    baseline=nanmean(prePoor,1);
lastDays=nanmean(baseline(8:10));
earlyDays=nanmean(baseline(1:7));
d=(lastDays-earlyDays)/earlyDays;
factor=round(d,1);
% factor=d;

pCompare=[];
for n=1:size(prePoor,1);
     a=nanmean(prePoor(n,8:10));
    b=nanmean(prePoor(n,1:7));
    pCompare(n,1)=(a-b)/b;
end


if factor>0.0
allIdx=find(pCompare(:,1)<=0.5);%need to remove the ones that the later part was larger than the before part
elseif factor<-0
allIdx=find(pCompare(:,1)>=-0.5);%need to remove the ones that the later part was smaller than the before part
else %factor==0: no corrections
    allIdx=[1:1:size(prePoor)];

end
allIdxAll{t}=allIdx;

subplot(2,4,t)
a=[prePoor durPoor postPoor];
% allIdx=size(all,1);
errorbar([1:1:30],nanmean(a(allIdx,:),1),nansem(a(allIdx,:),1),'k')
hold on
plot([11:1:20],nanmean(durPoor(allIdx,:),1),'b.','MarkerSize',20)
title(['N session',num2str(length(allIdx))])


subplot(2,4,t+4);
allSession=[];
allSession(:,1)=nanmean(prePoor(allIdx,:),2);
allSession(:,2)=nanmean(durPoor(allIdx,:),2);
allSession(:,3)=nanmean(postPoor(allIdx,:),2);
p=[];
[~,p(1)]=ttest(allSession(:,1),allSession(:,2));
[~,p(2)]=ttest(allSession(:,1),allSession(:,3));
[~,p(3)]=ttest(allSession(:,2),allSession(:,3));
[isSignificant,adjusted_pvals,~]= bonferroni_holm(p,0.05); % adjusted_pvals are the corrected p values
% adjusted_pvals=p;

for i=1:size(allSession,1);
    hold on
    plot([1:1:3],allSession(i,:),'m-')
end

bar([1:1:3],nanmean(allSession,1),'FaceColor','w');
hold on
errorbar([1:1:3],nanmean(allSession,1),nansem(allSession,1),'k.');
ylim([0 90])
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

saveas(gcf,'barGraphs.fig')
RBRCSadj=RBRCS(allIdxAll{1},:);
RBRadj=RBR(allIdxAll{2},:);
RBRCSGFPadj=RBRCSGFP(allIdxAll{3},:);
RBRRSGFPadj=RBRRSGFP(allIdxAll{4},:);
save('RBRCSadj.mat','RBRCSadj')
save('RBRadj.mat','RBRadj')
save('RBRCSGFPadj.mat','RBRCSGFPadj')
save('RBRRSGFPadj.mat','RBRRSGFPadj')
save('allIdxAll.mat','allIdxAll');
%% shading
load('RBRCSadj.mat');
load('RBRadj.mat');
load('RBRCSGFPadj.mat')
load('RBRRSGFPadj.mat')


RBR=RBRadj;
RBRCS=RBRCSadj;
RBRCSGFP=RBRCSGFPadj;
RBRRSGFP=RBRRSGFPadj;
N=32;

figure
subplot(221)
A=RBR(:,1:N);
for n=1:size(A,1);
    ma=mean(A(n,[1:10]));
    A(n,:)=A(n,:)/ma;
end

B=RBRCS(:,1:N);
for n=1:size(B,1);
    ma=mean(B(n,[1:10]));
    B(n,:)=B(n,:)/ma;
end
semshade(A,0.3,'k',[1:N])
hold on
semshade(B,0.3,'r',[1:N])

 [pAnova,pMC,pLabel] = anovaRM2W_full_BH(A(:,11:20),B(:,11:20),1);
 [pAnova2,pMC,pLabel] = anovaRM2W_full_BH(A(:,1:10),B(:,1:10),1);
 [pAnova3,pMC,pLabel] = anovaRM2W_full_BH(A(:,21:30),B(:,21:30),1);


% ylim([0.6 1.3])
title(['chr2 CS RS,P=',num2str([pAnova(1) pAnova2(1) pAnova3(1)])])
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
    ma=mean(A(n,[1:10]));
    A(n,:)=A(n,:)/ma;
end

B=RBRCSGFP(:,1:N);
for n=1:size(B,1);
    ma=mean(B(n,[1:10]));
    B(n,:)=B(n,:)/ma;
end
semshade(A,0.1,'k',[1:N])
hold on
semshade(B,0.1,'r',[1:N])
 [pAnova,pMC,pLabel] = anovaRM2W_full_BH(A(:,11:20),B(:,11:20),1);
 [pAnova2,pMC,pLabel] = anovaRM2W_full_BH(A(:,1:10),B(:,1:10),1);
 [pAnova3,pMC,pLabel] = anovaRM2W_full_BH(A(:,21:30),B(:,21:30),1);


% ylim([0.6 1.3])
title(['gfp CS RS, P=',num2str([pAnova(1) pAnova2(1) pAnova3(1)])])
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



subplot(224)
A=RBRCS(:,1:N);
for n=1:size(A,1);
    ma=mean(A(n,[1:10]));
    A(n,:)=A(n,:)/ma;
end

B=RBRCSGFP(:,1:N);
for n=1:size(B,1);
    ma=mean(B(n,[1:10]));
    B(n,:)=B(n,:)/ma;
end
semshade(A,0.3,'r',[1:N])
hold on
semshade(B,0.1,'r',[1:N])
 [pAnova,pMC,pLabel] = anovaRM2W_full_BH(A(:,11:20),B(:,11:20),1);
 [pAnova2,pMC,pLabel] = anovaRM2W_full_BH(A(:,1:10),B(:,1:10),1);
 [pAnova3,pMC,pLabel] = anovaRM2W_full_BH(A(:,21:30),B(:,21:30),1);


% ylim([0.6 1.3])
title(['CS CHR2 GFP, P=',num2str([pAnova(1) pAnova2(1) pAnova3(1)])])
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
    ma=mean(A(n,[1:10]));
    A(n,:)=A(n,:)/ma;
end

B=RBRRSGFP(:,1:N);
for n=1:size(B,1);
    ma=mean(B(n,[1:10]));
    B(n,:)=B(n,:)/ma;
end
semshade(A,0.3,'k',[1:N])
hold on
semshade(B,0.1,'k',[1:N])
 [pAnova,pMC,pLabel] = anovaRM2W_full_BH(A(:,11:20),B(:,11:20),1);
  [pAnova2,pMC,pLabel] = anovaRM2W_full_BH(A(:,1:10),B(:,1:10),1);
 [pAnova3,pMC,pLabel] = anovaRM2W_full_BH(A(:,21:30),B(:,21:30),1);


% ylim([0.6 1.3])
title(['RS CHR2 GFP, P=',num2str([pAnova(1) pAnova2(1) pAnova3(1)])])
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
saveas(gcf,'RBRTogether.fig')