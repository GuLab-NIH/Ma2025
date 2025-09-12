
%%reduce CS licking
load('idxCHR2RS.mat')
load('idxCHR2CS.mat')
load('idxGFPRS.mat')
load('idxGFPCS.mat')
load('predLickAfterRBR.mat')
load('predLickBeforeRBR.mat')
load('predLickStimRBR.mat')
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
if t==1;
% i=find(allSession(:,1)~=0);
% i=[4 5 6 7 8 9 10];
i=[1 2 3 4 5 6 7 8 13 14]; %can have 5
% 
% i=[1 2 3 4 5 7 8 13 14];%could be an option
allSession=allSession(i,:);
idxCHR2CS=idxCHR2CS(i);
idxCHR2CSLowerLick=idxCHR2CS;
save('idxCHR2CSLowerLick.mat','idxCHR2CSLowerLick')
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
saveas(gcf,'predictiveLickingPercentileRollingWindowPerSession_adjCSLevel.fig')
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
saveas(gcf,'RBRTogetherAdjCS.fig')

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
saveas(gcf,'RBRTogetherZeroAdjCS.fig')