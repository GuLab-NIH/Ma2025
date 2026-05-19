%% plotConsistencyExample_forRev
% Plots examples of spatial consistency information for all FOV and cell
% types


%% Load data and initialize settings

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
load([curFolder '\activity4m.mat'],'activity4m')
subLocRBR = activity4m;

% define categories
cellTypes = {'all','grid','cue','other'};
nCT = length(cellTypes);

% load cue info
load('cueTemp_4m.mat','cueTemp')
cueX = 2.5:5:397.5;
xBase = 2.5:5:377.5;

% define plot settings
colors = 'g';
nRuns = 15;

% define FOV, 1 example per mouse
% [1 2],[3 4],5,[6 7]
mMatch = [1 1 2 2 3 4 4];


%% Plot example spatial distribution

% define examples
useDay = [4 9 4 6; 4 10 7 9; 9 8 10 9];
useFOV = [1 3 5 6; 2 3 5 6; 2 4 5 7];
useCell = [3 80 8 21;8 17 2 6;11 66 2 27];

finalData = {};

figure; hold on
tiledlayout(3,4)
sgtitle('Spatial distribution examples')
for ct = 2:nCT
    for ex = 1:4
        curD = useDay(ct-1,ex);
        curF = useFOV(ct-1,ex);
        curC = useCell(ct-1,ex);
        curM = mMatch(curF);

        % combine all
        distData = subLocRBR.full.(cellTypes{ct}).LocRBR.good;
        distCat = cat(1,distData{curD,curF}(curC,:));

        % plot distribution
        nexttile(); hold on
        semshade(distCat,0.3,colors,xBase);

        plotCues(cueX,cueTemp,max(ylim),[0.5 0.5 0.5],0);

        xline(366)
        title([cellTypes{ct} ': D' num2str(curD) ', M' num2str(curM)...
            ', F' num2str(curF) ', C' num2str(curC)])
        ylim([-.2 1])

        finalData{ct-1,ex} = distCat;
    end
end

savefig([curFolder '\activitySpatialExamples_forRevision.fig'])


%% Plot example dfofM matrices

figure; hold on
tiledlayout(3,4)
sgtitle('Spatial distribution examples')

finalData = {};

for ct = 2:nCT
    for ex = 1:4
        curD = useDay(ct-1,ex);
        curF = useFOV(ct-1,ex);
        curC = useCell(ct-1,ex);
        curM = mMatch(curF);

        % combine all
        dfofMData = subLocRBR.dfofM.(cellTypes{ct}).LocRBR.good;
        dfofMEx = dfofMData{curD,curF}{curC}(1:min(nRuns,end),:);

        % plot distribution
        nexttile(); hold on
        imagesc(flipud(dfofMEx))
        axis('off')

        title([cellTypes{ct} ': D' num2str(curD) ', M' num2str(curM)...
            ', F' num2str(curF) ', C' num2str(curC)])

        finalData{ct-1,ex} = dfofMEx;
    end
end

savefig([curFolder '\activitySpatialExamplesMatrix_forRevision.fig'])

