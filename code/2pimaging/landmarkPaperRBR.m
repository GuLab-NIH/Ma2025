%% landmarkPaperRBR
% Calculate and plot RBR consistency as a function of track position for
% the landmark paper data. Used to compare pre-versus post-reward
% consistency.
%

clear; close all; clc

% move to base folder
p1 = 'D:\4m_optoPaper\RevisionAnalysis\LandmarkPaperData';
cd(p1)

% load data folders
load('Folders_For_Taylor.mat','Folders_For_Taylor')
nMouse = size(Folders_For_Taylor,1);
nDays = size(Folders_For_Taylor,2);
nFOV = size(Folders_For_Taylor,3);

% load cue template
load('cueTemp.mat')


%% Calculate spatial consistency

% % define data file
% dFile = 'RunByRun_dfof\corrInfoLocation.mat';
% 
% dataRBR = cell(nMouse,nDays,nFOV);
% 
% tic
% for mm = 1:nMouse
%     for dd = 1:nDays
%         for ff = 1:nFOV
%             % define current folder
%             curFolder = erase(Folders_For_Taylor(mm,dd,ff),'SNM\');
%             cd(curFolder)
% 
%             if isfile(dFile)
%                 load(dFile,'corrInfoLocation')
%                 dataRBR{mm,dd,ff} = corrInfoLocation.toOthersMean;
%             else
%                 dfofM_LocationConsistency
%             end
%         end
%     end
% end
% 
% cd(p1)
% save('dataRBR_lm.mat','dataRBR')
% toc


%%

% load data
load('dataRBR_lm.mat','dataRBR')

% calculate FOV means
meanRBR = cellfun(@(x) mean(x,'omitnan'),dataRBR,'UniformOutput',false);

% concatenate RBR by day
catCell = cell(1,nDays);
catMean = cell(1,nDays);
for dd = 1:nDays
    catCell{dd} = cat(1,dataRBR{:,dd,:});
    catMean{dd} = cat(1,meanRBR{:,dd,:});
end


%% Plot RBR distributions

useData = {catCell,catMean};
useLabs = {'by Cell','by FOV'};

colors2 = {[1 0 0];[0 0 1]};
labs2 = {'pre-reward','post-reward'};


% define cue info
nBins = length(temp);
cueX = (1:nBins)*5-2.5;

figure; tiledlayout(2,round(nDays))
X = (1:171)*5-2.5;
rewLoc = 63*5-2.5;
yLims = [0 0.08];

for ii = 1:2
    for dd = 1:nDays
        % current data
        plotData = useData{ii}{dd};

        % plot distribution
        nexttile(); hold on
        semshade(plotData,0.3,'g',X);

        plotCues(cueX,temp,yLims(2),'k',yLims(1));

        xline(rewLoc,'b')
        title(['Day ' num2str(dd) ': ' useLabs{ii}])
        xlabel('Track position')
        ylabel('RBR consistency')
        xlim([0 875])
        ylim(yLims)

        % add rectangles
        x1 = [275 305 305 275];
        x2 = [320 350 350 320];
        y = [yLims(1) yLims(1) yLims(2) yLims(2)];
        patch(x1, y, colors2{1}, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        patch(x2, y, colors2{2}, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    end
end


%% Plot RBR time course

figure; tiledlayout(1,2)
finalData = {};

for ii = 1:2
    % quantify distributions
    quantStruct = distrQuant(useData{ii},temp,62:64,2,[],4:6);

    % concatenate struct
    quantData = [quantStruct.RewPre;quantStruct.RewPost];

    % plot RBR time course
    nexttile(); hold on
    [~,pAll] = errorSig(1:nDays,quantData,colors2,labs2);

    % calculate correlation p values
    plotMean = cellfun(@(x) mean(x,'omitnan'),quantData);
    [r,p] = corr((1:nDays)',plotMean');
    disp(r)
    disp(p)

    xlabel('Day')
    ylabel('RBR consistency')
    xlim([0.5 nDays+0.5])
    title(['Spatial RBR Conistency (' useLabs{ii} ')'])

    if ii==2
        for jj = 1:2
            finalData{jj} = cat(2,quantData{jj,:});
        end
    end
end














