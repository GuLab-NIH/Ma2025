%% RBR example plots

clear; close all; clc

p1 = 'Z:\papers\2025RewardOptogenetics\4mData\250128_data';
cd(p1)

load('TM221026_day8_loc2_cell11_highRBR.mat')
load('ID20220628A_day1_loc1_cell46_lowRBR.mat')

load('ID20220628A_day1_loc1_cell46_lowRBR_dfof.mat')
load('TM221026_day8_loc2_cell11_highRBR_dfof.mat')


dfofM = {rbrLow_dfof,rbrHigh_dfof};
% dfofMAll = {dfofMAllLow,dfofMAllHigh};
% dfofSmooth = [dfofSmoothLow, dfofSmoothHigh];

load('D:\4m_optoPaper\RBRLoc.mat','RBRLoc')

% corrInfoLocation = dfofM_RBRbyLoc_Subset(dfofM,dfofSmooth);



%%

% define plot parameters
plotX = 12.5:5:387.5;
plotX = 0:5:375;

% define cue/reward info
colorsCue = [0.5 0.5 0.5];
colorsRew = [0 1 1];
cueX = 2.5:5:397.5;

% load cue/reward
load('cueTemp_4m.mat')
rewLoc = 366;

% calculate spatial RBR consistency
corrInfoLocation = dfofM_RBRbyLoc_Subset(dfofM,ones(80,2));
toOthers = corrInfoLocation.toOthersMean;

useDay = [1 8];
maxRuns = 15;

figure

for ii = 1:2
    %% Plot run-by-run activity
    subplot(4,2,ii)
    imagesc(dfofM{ii}(1:min(end,maxRuns),:))
    axis('off')
    ylabel('Run')

    % plot spatial RBR consistency
    subplot(4,2,ii+2)
    plot(plotX,toOthers(ii,:))
    ylabel('RBR correlation')
    xlim([0 400])

    
    subplot(4,2,ii+4)
    plot(mean(dfofM{ii},1))
    ylabel('dfof')



    % %% Plot single cell spatial rbr conistency
    % subplot(3,2,ii+2); hold on
    % ylim([-1 1])
    % ylabel('RBR correlation')
    % 
    % % plot data with shading
    % curData = corrInfoLocation.toOthersMean(ii,:);
    % plot(plotX,curData,'Color','k');
    % 
    % % plot cues
    % for jj = 1:size(cueTemp,2)
    %     plotCues(cueX,cueTemp,max(ylim),colorsCue,min(ylim));
    % end
    % 
    % % plot reward
    % xline(rewLoc,'Color',colorsRew)
    % 

    %% Plot population rbr conistency
    subplot(4,2,ii+6); hold on
    ylim([0 0.5])
    ylabel('RBR correlation')

    % plot data with shading
    curData = RBRLoc{useDay(ii),12};
    semshade(curData,0.3,'g',plotX);

    title(['Day ' num2str(useDay(ii))])

    % plot cues
    for jj = 1:size(cueTemp,2)
        plotCues(cueX,cueTemp,max(ylim),colorsCue,min(ylim));
    end

    set(gca,'FontSize',16)

    % plot reward
    xline(rewLoc,'Color',colorsRew)

end

savefig('rbrExample.fig')



