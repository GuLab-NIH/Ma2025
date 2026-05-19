%% calculateCrossLag
% calculate cross-lag correaltion between rolling RBR consitency and
% behavior

clear; close all; clc

p1 = 'D:\4m_optoPaper';
cd(p1)

% load folders
load('folders4m.mat')
nDays = size(allLocs,1);
nFolds = size(allLocs,2);

% define use mice
useMice = {[1 2], [8 9], 10, [11 12]};
useFOV = cat(2,useMice{:})';
behCode = [1 1 0 0 0 0 0 2 2 3 4 4];
nFOV = length(useFOV);

% define settings
rollRuns = 5;
fileNameRBR = ['rollRBRSpatial-' num2str(rollRuns) '.mat'];
fileNameLaps = 'lapIdx.mat';

% load cue template
load('cueTemp.mat');
bRewBin = 73;
beforeReward = 67:bRewBin; %before reward
temp = cueTemp;
temp(beforeReward) = 2;
tempUse = temp(1:bRewBin);

% define sublocations
idxs = struct();
idxs.all = 1:73;
idxs.outCue = find(tempUse==0);  % out cue bins
idxs.inCue = find(tempUse==1);   % in cue bins (excludes reward)
idxs.preRew = find(tempUse==2);  % pre-reward bins

% define location categories
locCats = {'all','outCue','inCue','preRew'};
nCats = length(locCats);

% initialize data array
data = struct();
data.RBR = cell(nDays,nFOV);
for cc = 1:nCats
    data.(['RBR_' locCats{cc}]) = cell(nDays,nFOV);
end
data.Slow = cell(nDays,nFOV);
data.Lick = cell(nDays,nFOV);
data.LickB = cell(nDays,nFOV);
data.LickC = cell(nDays,nFOV);
data.Merge = cell(nDays,nFOV);

% load behavior data
load(['RevisionAnalysis\data\behRolling_' num2str(rollRuns) '.mat'],'dataBeh')

for ff = 1:nFolds
    ffIdx = find(useFOV==ff);
    if isempty(ffIdx)
        continue
    end

    for dd = 1:nDays
        % move to current folder
        cd(allLocs{dd,ff})


        %% Calculate rolling RBR

        % get rolling RBR
        if isfile(fileNameRBR)
            load(fileNameRBR,'rollRBRSpatial')
        else
            % load dfof
            load('RunByRun_dfof\corrInfoLocation.mat','corrInfoLocation')

            rollRBRSpatial = temporalRBR_Spatial(corrInfoLocation,rollRuns);
            save(fileNameRBR,'rollRBRSpatial')
            copyfile('D:\4m_optoPaper\RevisionAnalysis\temporalRBR_Spatial.m','postAnalysis\temporalRBR_Spatial.m')
        end

        % store RBR data
        data.RBR{dd,ffIdx} = rollRBRSpatial;

        
        %% Split data by spatial location

        for cc = 1:nCats
            % get current sub-location data
            dataTemp = rollRBRSpatial(:,idxs.(locCats{cc}));
            dataMean = mean(dataTemp,2,'omitnan');

            data.(['RBR_' locCats{cc}]){dd,ffIdx} = dataMean;
        end


        %% Calculate rolling behavior

        % get true lap indices
        if isfile(fileNameLaps)
            load(fileNameLaps,'imageLapIdx')
        else
            disp(pwd)
            try
                imageLapIdx = getImagingLaps();
                save(fileNameLaps,'imageLapIdx')
                copyfile('D:\4m_optoPaper\RevisionAnalysis\getImagingLaps.m','postAnalysis\getImagingLaps.m')
            catch
                % try
                %     cpFold = erase(replace(pwd,p1,'G:\Mice\Raw'),'\pcaica');
                %     d = dir([cpFold '\T*']);
                %     folderPath = [d(1).folder '\' d(1).name];
                %     f2 = '*VoltageRecording*.csv';
                %     d2=dir([folderPath '\' f2]);
                %     m = csvread([folderPath '\' d2(1).name],2,0);
                %     save([folderPath '\m.mat'],'m','-v7.3');
                %
                %     copyfile([folderPath '\m.mat'],pwd)
                % catch
                % end
                error('Image laps not synced')
            end
        end

        % get mouse behavior
        mSlow = dataBeh.Slow{dd,behCode(ff)};
        mLick = dataBeh.Lick{dd,behCode(ff)};
        mLickB = dataBeh.LickB{dd,behCode(ff)};
        mLickC = dataBeh.LickC{dd,behCode(ff)};

        % sync behavior and imaging parameters
        % try
        data.Slow{dd,ffIdx} = mSlow(imageLapIdx(1:end-rollRuns+1))/100;
        data.Lick{dd,ffIdx} = mLick(imageLapIdx(1:end-rollRuns+1));
        data.LickB{dd,ffIdx} = mLickB(imageLapIdx(1:end-rollRuns+1));
        data.LickC{dd,ffIdx} = mLickC(imageLapIdx(1:end-rollRuns+1));
        data.Merge{dd,ffIdx} = mean([data.Slow{dd,ffIdx},data.LickB{dd,ffIdx}],2);

        % return to original folder
        cd(p1)
    end
end


%% Calculate cross lag correlation

% define variables to run
varBeh = {'Slow','LickC','LickB','Merge'};
nBeh = length(varBeh);

% set cross lag parameters
maxLag = 0;
lagFact = 2;
useDiff = 0;
useNorm = 1;
nanLim = 2;

% set plot parameters
colsAct = [0 0 0];
colsBeh = {[0 1 0],[1 0 1],[0 1 1],[1 0 0]};

% intialize data output
dataLags = struct();

% nDays = 5;
% nFOV = 4;

for cc = 1:nCats
    for bb = 1:nBeh
        % initialize current output
        dataLags.(varBeh{bb}) = zeros(nDays,nFOV);

        % initialize visualization
        figure; tiledlayout(nDays,nFOV)
        sgtitle([locCats{cc} ' RBR vs. ' varBeh{bb}])

        for dd = 1:nDays
            for ff = 1:nFOV
                nexttile(); hold on

                % get current rolling averages
                curAct = data.(['RBR_' locCats{cc}]){dd,ff};
                curBeh = data.(varBeh{bb}){dd,ff}(1:length(curAct));

                % skip comparisons with too many NaNs
                if sum(isnan(curBeh))>length(curBeh)/nanLim
                    dataLags.(varBeh{bb})(dd,ff) = NaN;
                    continue
                end

                % take difference
                if useDiff==1
                    curAct = diff(curAct);
                    curBeh = diff(curBeh);
                end

                % normalize vectors
                if useNorm==1
                    if length(unique(curAct(~isnan(curAct))))==1
                        curActN = curAct;
                        curActN(~isnan(curActN)) = 0;
                    else
                        curActN = normalize(curAct);
                    end

                    if length(unique(curBeh(~isnan(curBeh))))==1
                        curBehN = curBeh;
                        curBehN(~isnan(curBehN)) = 0;
                    else
                        curBehN = normalize(curBeh);
                    end
                else
                    curActN = curAct;
                    curBehN = curBeh;
                end

                % determine use max lag
                if maxLag==0
                    useMaxLag = round(length(curActN)/lagFact);
                else
                    useMaxLag = maxLag;
                end

                % perform cross correlation
                [corrAB,lags] = xcov(curActN,curBehN,useMaxLag,'normalized');
                if all(isnan(corrAB)) || length(unique(corrAB))==1
                    curLag = NaN;
                else
                    [mxCorr,mxIdx] = max(corrAB);
                    curLag = lags(mxIdx);
                end

                % store cross lag
                dataLags.locCats{cc}.(varBeh{bb})(dd,ff) = curLag;

                % plot curves
                title(['Lag = ' num2str(curLag)])
                plot(curActN,'Color',colsAct)
                plot(curBehN,'Color',colsBeh{bb})
                % legend(varAct,varBeh{bb},'Location','best')
            end
        end

        close

    end
end


%% Plot lag distributions

f1 = figure; tiledlayout(nCats,nBeh)
f2 = figure; tiledlayout(nCats,nBeh)
f3 = figure; tiledlayout(nCats,nBeh)

sgtitle(['Rolling Lap Cross-Correlation (' num2str(rollRuns) ' runs)'])
rng(42)

for cc = 1:nCats
    for bb = 1:nBeh
        % get current distribution
        curDist = dataLags.locCats{cc}.(varBeh{bb})(:);

        % calculate proportion negative
        testDists = curDist(~isnan(curDist));
        obs = mean(testDists<0);

        % calculate expected proportion negative
        nperm = 10000;
        null = zeros(nperm,1);
        for i = 1:nperm
            flips = sign(randn(size(testDists)));
            null(i) = mean(testDists.*flips<0);
        end

        % calculate sign permuation p-value
        pCur = mean(null>obs)+mean(null==obs);

        % perform t-test
        % [~,pCur] = ttest(testDists);

        % plot current distribution
        figure(f1)
        nexttile((cc-1)*nBeh+bb); hold on
        [ks,Xi] = ksdensity(curDist,'Function','pdf');
        plot(Xi,ks,'Color',[0 0 1])

        % plot histogram
        % histogram(curDist)

        % plot mean and zero
        % curM = median(curDist,'omitnan');
        % xline(curM,'Color',[0 0 1])
        xlim([-15 15])
        xline(0,':k')
        % yline(0.5,':k')
        ylim([0 0.2])

        % set figure labels
        ylabel('Probability density function')
        xlabel('Lag')
        title([locCats{cc} ' vs. ' varBeh{bb} ': % neg = ' num2str(obs,2)])


        %% Plot shuffles

        % plot shuffle
        figure(f2)
        nexttile((cc-1)*nBeh+bb); hold on

        % plot shuffles
        [yi,xi] = ksdensity(null,0:0.01:1);
        plot(xi,yi,'Color',[0.5 0.5 0.5])

        % plot true data
        xline(obs,'Color',[0 0 1])

        % set figure labels
        ylabel('Probability density function')
        xlabel('% negative lags')
        title([varBeh{bb} ': p = ' num2str(pCur,2)])


        %% Plot pie chart

        figure(f3)
        nexttile((cc-1)*nBeh+bb); hold on

        % generate pie data
        pieData = [sum(testDists<0) sum(testDists==0) sum(testDists>0)];

        % plot pie chart
        pie(pieData)

        % define plot settings
        axis('square','off')
        legend('neagtive','zero','positive','Location','best')
        title([varBeh{bb} ': p = ' num2str(pCur,2)])
    end
end

