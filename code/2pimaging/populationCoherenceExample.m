%% calculatePopulationCoherence
% calculates population coherence across days and analyzes trends

clear; clc
close all

p1 = 'D:\4m_optoPaper';
cd(p1)

% load folders
load('folders4m.mat')
nDays = size(allLocs,1);
nFOV = size(allLocs,2);

% define settings
nBins = 80;
suff = 'dfof';
fNames = 'PVC';
curCMap = customcolormap([0 0.25 0.5 0.75 1],[0 0 1; 0.5 0.5 1; 1 1 1; 1 0.5 0.5; 1 0 0]);
cLims = [-1 1];
plotLabs = {'Single Lap','Mean Data'};

% define exmaple
useRuns = [1 2];

% initialize data array
data = struct();


for ff = nFOV
    for dd = nDays
        % move to current folder
        cd(allLocs{dd,ff})

        % load use rois
        load('roiMasks.mat','masks')
        useMask = masks(:,1);


        % load dfof
        load('RunByRun_dfof\dfofMInterpM.mat')
        dfof = cat(3,dfofMInterpM{:});

        % get lap-to-lap PVC
        [meanPVC,~,PVC] = popCoherence(dfof(:,:,useMask));

        % extract single PVC
        snPVC = PVC(:,:,useRuns(1),useRuns(2));

        %% Plot example figure

        figure; tiledlayout(1,2)
        plotData = {snPVC,meanPVC};

        for ii = 1:2
            % calculate 2D mean
            plotMap = plotData{ii};

            % plot population coherence heatmap
            nexttile(ii); hold on

            imAlpha = ones(size(plotMap));
            imAlpha(isnan(plotMap))=0;
            imagesc(plotMap,'AlphaData',imAlpha);
            set(gca,'color',0*[1 1 1]);
            set(gca,'FontSize',15);

            % set axes
            xlim([0.5 nBins+0.5])
            ylim([0.5 nBins+0.5])
            clim(cLims)
            colormap(curCMap)
            axis('square')
            set(gca,'XTick',[],'YTick',[],'XColor','none','YColor','none','Box','off')
            title([fNames ': ' plotLabs{ii}])
        end
    end
end
