%% calculateLocationConsistency
% combines spatial consistency information for all FOV and calculates
% change over learning for spatial subsets


%% Idenitfy pcaica folders for all locs/days

clear; clc; close all

% set current folder
p1 = 'D:\4m_optoPaper';
cd(p1)

% save folder
curFolder = 'SpatialActivity\forFigsRev';
if ~isfolder(curFolder)
    mkdir(curFolder)
end

% load data
load([curFolder '\locDataRev.mat'],'rawData')
load([curFolder '\activity4m.mat'],'activity4m')
subLocRBR = activity4m;

% define categories
cellTypes = {'all','grid','cue','other'};
distTypes = 'LocRBR';
subLocLabs = {'preRew','inCue','outCue'};
nCT = length(cellTypes);
nDT = length(distTypes);
nSub = length(subLocLabs);

% define groups
groupMice = {[1 2],[3 4],5,[6 7]};
groupIDs = 'good';
colors = 'g';
nFOV = length(cat(2,groupMice{:}));
nMice = length(groupMice);

% load cue info
load('cueTemp_4m.mat','cueTemp')
cueX = 2.5:5:397.5;
xBase = 2.5:5:377.5;

% load behavior
load('behavior4m.mat','behavior4m')
behData = struct();
behData.slowing = behavior4m.slow.good;
load('Z:\papers\2025RewardOptogenetics\4mData\VirmenLogs\allLickEachDay.mat','allLickEachDay')
behData.licking = allLickEachDay;
% behData.licking = behavior4m.slow.good;
behTypes = {'slowing','licking'};
nBeh = length(behTypes);

% define averaging methods
typeA = {'byCell','byFOV','byMouse'};
nA = length(typeA);


%% Extract sublocation shuffles

binOffSet = 0;
rewCutSize = 0;

% define zones
bRewBin = 73;
beforeReward = 67:bRewBin-rewCutSize; %before reward
temp = cueTemp;
temp(beforeReward) = 2;
tempUse = temp(1+binOffSet:bRewBin-rewCutSize);

% define sublocations
outCue = find(tempUse==0);  % out cue bins
inCue = find(tempUse==1);   % in cue bins (excludes reward)
preRew = find(tempUse==2);  % pre-reward bins
subLocBins = {preRew,inCue,outCue};

% initialize shuffle settings
rng(42)
nShuff = 200;
nBins = 76;
rndBins = zeros(nShuff,length(preRew));
for rr = 1:nShuff
    curRnd = randperm(nBins);
    rndBins(rr,:) = curRnd(preRew);
end

% initialize data struct
splitRBR = struct();

for ct = 1:nCT
    % load current data
    curFull = subLocRBR.full.(cellTypes{ct}).(distTypes).(groupIDs);

    % calculate sub regions per FOV
    for ss = 1:nSub
        splitRBR.(cellTypes{ct}).(subLocLabs{ss}) = cell(size(curFull));
        for ff = 1:numel(curFull)
            curFOV = mean(curFull{ff}(:,subLocBins{ss}),2,'omitnan');
            splitRBR.(cellTypes{ct}).(subLocLabs{ss}){ff} = curFOV;
        end
    end

    % initialize shuffle output
    splitRBR.(cellTypes{ct}).shuffle = cell(size(curFull));
    for ff = 1:numel(curFull)
        curFOV = curFull{ff};
        nCells = size(curFOV,1);

        % shuffle current FOV
        shuffCur = zeros(nCells,nShuff);
        for rr = 1:nShuff
           shuffCur(:,rr) = mean(curFOV(:,rndBins(rr,:)),2,'omitnan');
        end

        % take average of shuffles
        splitRBR.(cellTypes{ct}).shuffle{ff} = mean(shuffCur,2,'omitnan');
    end
end


%% Plot all cell spatial conistency across learning

% define plot info
useDays = 1:10;
nDays = length(useDays);

finalData = struct();

for ct = 1:nCT
    % initialize figure
    figure; tiledlayout (nA,nDays)
    sgtitle(['Spatial Consistency Distribution ('  cellTypes{ct} ' cells)'])

    curCellsN = [];

    for aa = 1:nA
        % combine all
        distData = subLocRBR.full.(cellTypes{ct}).(distTypes).(groupIDs);

        for dd = 1:nDays
            curDay = useDays(dd);

            % initialize panel
            nexttile(); hold on
            title([typeA{aa} ': Day ' num2str(curDay)])

            if strcmp(typeA{aa},'byCell')
                % use all cells
                distCat = cat(1,distData{curDay,:});
            elseif strcmp(typeA{aa},'byFOV')
                % average by FOV and concatenate
                distMean = cellfun(@(x) mean(x,1,'omitnan'),distData(curDay,:),...
                    'UniformOutput',false);
                distCat = cat(1,distMean{:});
            elseif strcmp(typeA{aa},'byMouse')
                % concatenate and average by mouse and concatenate
                distMouse = cell(1,nMice);
                for mm = 1:nMice
                    distMouse{mm} = cat(1,distData{curDay,groupMice{mm}});
                end
                distMean = cellfun(@(x) mean(x,1,'omitnan'),distMouse,'UniformOutput',false);
                distCat = cat(1,distMean{:});

                % adjust title
                curDM = cellfun(@(x) size(x,1),distMouse);
                curDMM = round(mean(curDM,'all'));
                curDMS = round(nansem(curDM,'all'));
                title({[typeA{aa} ': Day ' num2str(curDay)],...
                    ['n = ' num2str(curDMM) ' +/- ' num2str(curDMS) ' cells']})
            else
                error('Invalid averaging type')
            end

            % plot distribution
            semshade(distCat,0.3,colors,xBase);
            plotCues(cueX,cueTemp,max(ylim),[0.5 0.5 0.5],0);
            xline(366)
            ylim([0 0.5])

            % add cell counting
            if aa==3
                curCellsN = [curCellsN,curDM];
            end

            finalData.(cellTypes{ct}){aa,dd} = distCat;
        end
    end

    % display current cell numbers
    globDMM = round(mean(curCellsN,'all'));
    globDMS = round(nansem(curCellsN,'all'));
    disp([cellTypes{ct} ' cells: n = ' num2str(globDMM) ' +/- ' num2str(globDMS)])

    savefig([curFolder '\activitySpatial_learning_' cellTypes{ct} '.fig'])
end


%% Plot single-day spatial conistency by cell type

% define plot info
useDay = 9;

% initialize figure
figure; tiledlayout (nA,nCT)
sgtitle(['Spatial Consistency Distribution: Day ' num2str(useDay)])

for aa = 1:nA
    for ct = 1:nCT
        nexttile(); hold on
        title([typeA{aa} ' ' cellTypes{ct}])

        % combine all
        distData = subLocRBR.full.(cellTypes{ct}).(distTypes).(groupIDs);

        if strcmp(typeA{aa},'byCell')
            % use all cells
            distCat = cat(1,distData{useDay,:});
        elseif strcmp(typeA{aa},'byFOV')
            % average by FOV and concatenate
            distMean = cellfun(@(x) mean(x,1,'omitnan'),distData(useDay,:),...
                'UniformOutput',false);
            distCat = cat(1,distMean{:});
        elseif strcmp(typeA{aa},'byMouse')
            % concatenate and average by mouse and concatenate
            distMouse = cell(1,nMice);
            for mm = 1:nMice
                distMouse{mm} = cat(1,distData{useDay,groupMice{mm}});
            end
            distMean = cellfun(@(x) mean(x,1,'omitnan'),distMouse,'UniformOutput',false);
            distCat = cat(1,distMean{:});
        else
            error('Invalid averaging type')
        end

        % plot distribution
        semshade(distCat,0.3,colors,xBase);
        plotCues(cueX,cueTemp,max(ylim),[0.5 0.5 0.5],0);
        xline(366)
        ylim([0 0.4])

    end
end

savefig([curFolder '\activitySpatial_D' num2str(useDay) '.fig'])


%% Plot learning curves

% define curves
subLabs = fieldnames(splitRBR.(cellTypes{1}));
colors4 = {[0 0 1],[1 1 0],[1 1 1],[0.5 0.5 0.5]};
nSubs = length(subLabs);
nDays = size(splitRBR.(cellTypes{1}).(subLabs{1}),1);

% initialize bar data
barData = struct();
finalData = struct();

% initialize figure
figure; tiledlayout(nA,nCT)
sgtitle('Learning Curves')

for aa = 1:nA
    yLims = [NaN NaN];

    for ct = 1:nCT
        % set current axis
        nexttile((nCT*(aa-1)+ct)); hold on
        title([typeA{aa} ' ' cellTypes{ct}])

        for ss = 1:nSubs
            % get current data
            plotData = splitRBR.(cellTypes{ct}).(subLabs{ss});

            if strcmp(typeA{aa},'byCell')
                % average by cell (i.e. unchanged)
                plotMean = plotData;
            elseif strcmp(typeA{aa},'byFOV')
                % average by FOV
                plotMean = cellfun(@(x) mean(x,1,'omitnan'),plotData,...
                    'UniformOutput',false);
            elseif strcmp(typeA{aa},'byMouse')
                % concatenate and average by mouse
                plotMouse = cell(size(plotData,1),nMice);
                for mm = 1:nMice
                    for dd = 1:nDays
                        plotMouse{dd,mm} = cat(1,plotData{dd,groupMice{mm}});
                    end
                end
                plotMean = cellfun(@(x) mean(x,1,'omitnan'),plotMouse,'UniformOutput',false);

                % pring number of cells
                if ss==1
                    nCellPerMouse = cellfun(@length,plotMouse);
                    curDMM = round(mean(nCellPerMouse,'all'));
                    curDMS = round(nansem(nCellPerMouse,'all'));
                    disp([typeA{aa} '-' cellTypes{ct} ': n = ' num2str(curDMM) ' +/- ' num2str(curDMS) ' cells'])
                end
            else
                error('Invalid averaging type')
            end

            plotCat = cell(nDays,1);
            for dd = 1:nDays
                plotCat{dd} = cat(1,plotMean{dd,:});
            end

            % plot current data
            plotM = cellfun(@(x) mean(x,1,'omitnan'),plotCat);
            plotSEM = cellfun(@(x) nansem(x,1),plotCat);
            errorbar(plotM,plotSEM,'Color',colors4{ss})

            % set plot labels
            xlabel('Days in NE')
            ylabel('RBR consistency')
            fontsize(12,'points')
            
            % calculate and print correlation
            [rCur,pCur] = corr((1:length(plotM))',plotM);
            fprintf('%s %s %s: r = %.2f, p = %.2g\n',typeA{aa},...
                cellTypes{ct},subLabs{ss},rCur,pCur)

            % store data for bar plot
            barData.(typeA{aa}).(cellTypes{ct}).(subLabs{ss}) = cat(1,plotCat{:});

            if aa==1
                finalData.(typeA{aa}).(cellTypes{ct}).(subLabs{ss}) = plotCat;
            else
                finalData.(typeA{aa}).(cellTypes{ct}).(subLabs{ss}) = cat(2,plotCat{:});
            end
        end

        % update y limits
        yLims(1) = min(yLims(1),min(ylim),'omitnan');
        yLims(2) = max(yLims(2),max(ylim),'omitnan');
    end

    % set global y limits
    for ct = 1:nCT
        nexttile((nCT*(aa-1)+ct))
        xlim([0.5 11.5])
        % ylim(yLims)
        ylim([0 0.4])
        set(gca,'Color','k')
    end

end

savefig([curFolder '\activityTimecourse.fig'])


%% Compare learning curves to behavior

% define behavior colors
colors2 = {'g','m'};

for bb = 1:nBeh
    curBeh = behData.(behTypes{bb});

    % initialize figure
    figure; tiledlayout(nA,nSubs)
    sgtitle(['RBR Consistency vs. ' behTypes{bb}])

    for aa = 1:nA
        for ss = 1:nSubs
            % set current axis
            nexttile((nSubs*(aa-1)+ss)); hold on

            % get current data
            plotData = splitRBR.all.(subLabs{ss});

            if strcmp(typeA{aa},'byCell')
                % average by cell (i.e. unchanged)
                plotMean = plotData;
            elseif strcmp(typeA{aa},'byFOV')
                % average by FOV
                plotMean = cellfun(@(x) mean(x,1,'omitnan'),plotData,...
                    'UniformOutput',false);
            elseif strcmp(typeA{aa},'byMouse')
                % concatenate and average by mouse
                plotMouse = cell(size(plotData,1),nMice);
                for mm = 1:nMice
                    for dd = 1:nDays
                        plotMouse{dd,mm} = cat(1,plotData{dd,groupMice{mm}});
                    end
                end
                plotMean = cellfun(@(x) mean(x,1,'omitnan'),plotMouse,'UniformOutput',false);
            else
                error('Invalid averaging type')
            end

            plotCat = cell(nDays,1);
            for dd = 1:nDays
                plotCat{dd} = cat(1,plotMean{dd,:});
            end

            % plot RBR data
            yyaxis left
            AplotM = cellfun(@(x) mean(x,1,'omitnan'),plotCat);
            AplotSEM = cellfun(@(x) nansem(x,1),plotCat);
            errorbar(AplotM,AplotSEM,'Color','k')
            ylabel('RBR consistency')
            
            % plot RBR data
            yyaxis right
            BplotM = mean(curBeh,1,'omitnan')';
            BplotSEM = nansem(curBeh,1)';
            errorbar(BplotM,BplotSEM,'Color',colors2{bb})
            ylabel(['Predictive ' behTypes{bb}])

            % set plot labels
            xlabel('Days in NE')
            fontsize(12,'points')
            xlim([0 11])

            % calculate correlation
            [rCur,pCur] = corr(AplotM,BplotM);

            % set title
            title([typeA{aa} ' ' subLabs{ss} ' (r = ' num2str(rCur,'%.2f')...
                ', p = ' num2str(pCur,'%.2g') ')']);
        end
    end

    savefig([curFolder '\activityTimecourse_' behTypes{bb} '.fig'])
end


%% Compare learning curves to behavior by mouse

% define behavior colors
colors2 = {'g','m'};

for bb = 1:nBeh
    curBeh = behData.(behTypes{bb});

    % initialize figure
    figure; tiledlayout(nMice,nSubs)
    sgtitle(['RBR Consistency vs. ' behTypes{bb}])

    for ss = 1:nSubs
        % get current data
        plotData = splitRBR.all.(subLabs{ss});

        % concatenate and average by mouse
        plotMouse = cell(size(plotData,1),nMice);
        for mm = 1:nMice
            for dd = 1:nDays
                plotMouse{dd,mm} = cat(1,plotData{dd,groupMice{mm}});
            end
        end
        plotMean = cellfun(@(x) mean(x,1,'omitnan'),plotMouse);

        for mm = 1:nMice
            % set current axis
            nexttile((nSubs*(mm-1)+ss)); hold on

            plotMean(:,mm);

            % plot RBR data
            yyaxis left
            plot(AplotM,'Color','k')
            ylabel('RBR consistency')
            
            % plot RBR data
            yyaxis right
            BplotM = curBeh(mm,:)';
            plot(BplotM,'Color',colors2{bb})
            ylabel(['Predictive ' behTypes{bb}])

            % set plot labels
            xlabel('Days in NE')
            fontsize(12,'points')
            xlim([0 11])

            % calculate correlation
            [rCur,pCur] = corr(AplotM,BplotM);

            % set title
            title(['Mouse ' num2str(mm) ' ' subLabs{ss} ' (r = ' num2str(rCur,'%.2f')...
                ', p = ' num2str(pCur,'%.2g') ')']);
        end
    end

    savefig([curFolder '\activityTimecourse_' behTypes{bb} '_byMouse.fig'])
end


%% Compare learning curves to behavior by mouse correlation

% define behavior colors
colors2 = {'g','m'};

% initialize figure
figure; tiledlayout(nBeh,nSubs)
sgtitle('RBR Consistency vs. Behavior')

for bb = 1:nBeh
    curBeh = behData.(behTypes{bb});

    for ss = 1:nSubs
        % get current data
        plotData = splitRBR.all.(subLabs{ss});

        % concatenate and average by mouse
        plotMouse = cell(size(plotData,1),nMice);
        for mm = 1:nMice
            for dd = 1:nDays
                plotMouse{dd,mm} = cat(1,plotData{dd,groupMice{mm}});
            end
        end
        plotMean = cellfun(@(x) mean(x,1,'omitnan'),plotMouse);

        % set current axis
        nexttile((nSubs*(bb-1)+ss)); hold on

        % plot scatter
        scatter(plotMean,curBeh')

        % set plot labels
        xlabel('RBR consistency')
        ylabel(behTypes{bb})
        fontsize(12,'points')

        % calculate correlation
        [rCur,pCur] = corr(plotMean(:),curBeh(:));

        % set title
        title([behTypes{bb} ' ' subLabs{ss} ' (r = ' num2str(rCur,'%.2f')...
            ', p = ' num2str(pCur,'%.2g') ')']);
    end
end

savefig([curFolder '\activityBehaviorCorrelation.fig'])


%% Plot pooled bar graph

% define plot parameters
pShow = [1 2;1 3;1 4;2 3];
gOrder = [1 2 3 4];
yLimsUse = {[0 0.3],[0 0.4],[0 0.4],[0 0.4]};

% initialize figure
figure; tiledlayout(nA,nCT)
sgtitle('Pooled spatial activity')

for aa = 1:nA
    yLims = [NaN NaN];

    for ct = 1:nCT
        % set current axis
        nexttile((nCT*(aa-1)+ct)); hold on
        title([typeA{aa} ' ' cellTypes{ct}])

        plotData = cell(1,nSubs);
        for ss = 1:nSubs
            % get current data
            plotData{ss} = barData.(typeA{aa}).(cellTypes{ct}).(subLabs{ss});
        end

        % plot pooled bar
        [~,pBar] = barGroup(subLabs(gOrder),plotData(gOrder),'bar',colors4,pShow,'pair');
        fprintf('%s %s\n',typeA{aa},cellTypes{ct})
        disp(num2str(pBar,3))
        disp(' ')

        % set plot labels
        ylabel('RBR consistency')
        fontsize(12,'points')

        % update y limits
        yLims(1) = min(yLims(1),min(ylim),'omitnan');
        yLims(2) = max(yLims(2),max(ylim),'omitnan');
    end

    % set global y limits
    for ct = 1:nCT
        nexttile((nCT*(aa-1)+ct))
        % ylim(yLims)
        ylim(yLimsUse{ct})
    end

end

savefig([curFolder '\activityBar.fig'])


%% Plot prereward percentile with respect to other locations

% define averaging methods
typeAb = {'bySessionDay','byFOV','byMouse'};
nAb = length(typeAb);

figure; tiledlayout(nA,nCT)
sgtitle('Pre-reward percentile')

finalData = {};

for ct = 1:nCT
    % get expanded activity by bin
    perData = subLocRBR.expand.(cellTypes{ct}).(distTypes).good;

    % get pre-reward mean and non-reward concatenation
    perDataRew = cellfun(@(x) mean(x,2,'omitnan'),perData.preRew,'UniformOutput',false);
    perDataNon = cellfun(@(x,y) cat(2,x,y),perData.outCue,perData.inCue,'UniformOutput',false);

    % calculate percentile by cell
    perRew = cellfun(@(x,y) invPercentile(x',y')',perDataRew,perDataNon,'UniformOutput',false);

    % calculate cdf per fov
    [cdfY,cdfX] = cellfun(@(x) ecdf(x),perRew,'UniformOutput',false);

    % interpolate cdf
    intX = 0:0.1:100;
    nX = length(intX);
    intY = cellfun(@(x,y) interp1(x(2:end),y(2:end),intX),cdfX,cdfY,'UniformOUtput',false);

    % correct NaNs at extreme values
    for ii = 1:numel(intY)
        YY = intY{ii};

        % correct start NaNs
        idxStart = find(~isnan(YY),1);
        if idxStart>1
            YY(1:idxStart-1) = 0;
        end

        % correct end NaNs
        idxEnd = find(~isnan(YY),1,'last');
        if idxEnd<nX
            YY(idxEnd+1:nX) = 1;
        end

        % check for data errors
        if any(YY(1:end-1)>YY(2:end)) || any(isnan(YY))
            disp(ii)
            return
        end

        intY{ii} = YY;
    end

    for aa = 1:nAb
        if strcmp(typeAb{aa},'bySessionDay')
            % average by session-day (i.e. unchanged)
            intPlot = intY;
        elseif strcmp(typeAb{aa},'byFOV')
            % average by FOV
            intPlot = cell(nFOV,1);
            for ff = 1:nFOV
                intPlot{ff} = mean(cat(1,intY{ff,:}),1,'omitnan');
            end
        elseif strcmp(typeAb{aa},'byMouse')
            % concatenate and average by mouse
            intPlot = cell(nMice,1);
            for mm = 1:nMice
                intPlot{mm} = mean(cat(1,intY{groupMice{mm},:}),1,'omitnan');
            end
        else
            error('Invalid averaging type')
        end

        % concatenate interpolated cdf
        cdfPlot = cat(1,intPlot{:});

        % plot groups data with shading
        nexttile((aa-1)*nCT+ct); hold on
        title(['CDF: ' typeAb{aa} ' ' cellTypes{ct}])
        semshade(cdfPlot,0.3,colors,intX);

        finalData{ct,aa} = cdfPlot';
    end
end

savefig([curFolder '\activityPercentile.fig'])

