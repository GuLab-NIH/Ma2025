%% GLM_4m
% Use GLMs to predict mouse behavior as a function of activity.
%

clear; close all; clc

p1 = 'D:\4m_optoPaper';
cd(p1)

% load data
load('RevisionAnalysis\dataRevision_vRPR.mat','data')
fNames = fieldnames(data.spatial);
nFields = length(fNames);
nDays = size(data.spatial.(fNames{1}),1);

% define use bin range
useBins = 1:73;
nBins = length(useBins);

% set permuation number
rng(42)
nPerm = 100;
nBX = 6;
nBY = 6;

% select activity mice
actMice = {[1 2], [8 9], 10, [11 12]};
nMice = length(actMice);
behStart = 5;

% define track positions
binX = abs((2.5:5:397.5)-366);

% generate data structure
dataSp = struct();
for ff = 1:nFields
    curData = data.spatial.(fNames{ff})(:,:,useBins);

    % process data
    if ff<behStart
        curSp = zeros(nDays,nMice,nBins);

        % take mouse averages
        for mm = 1:nMice
            curSp(:,mm,:) = mean(curData(:,actMice{mm},:),2);
        end
    else
        % keep mouse data unchanged
        curSp = curData;
    end

    % remove infinities
    curSp(abs(curSp)==Inf) = NaN;

    % store data
    dataSp.(fNames{ff}) = curSp;
end

% generate grouping variables
dataSp.session = categorical(repmat((1:nDays)',[1 nMice nBins]));
dataSp.mouse = categorical(repmat(1:nMice,[nDays 1 nBins]));
dataSp.rewDist = repmat(permute(binX(useBins),[1 3 2]),[nDays nMice 1]);
fNames = fieldnames(dataSp);
nFields = length(fNames);

% generate table
normIdx = [1:6 11];
tableSp = table();
for ff = 1:nFields
    if ismember(ff,normIdx)
        tableSp.(fNames{ff}) = normalize(dataSp.(fNames{ff})(:));
    else
        tableSp.(fNames{ff}) = dataSp.(fNames{ff})(:);
    end
end

% correct lick numbers
tableSp.LickingCount = round(tableSp.Licking.*tableSp.laps);
fNamesTbl = tableSp.Properties.VariableNames;


%% Perform GLMM analysis to predict behavior from activity

rng(42)

dataMI = struct();

% define activity parameters
colsX = 4;

% define behavior parameters
labsBeh = {'velocity','licking'};
colsPred = {[6 7 11],[5 6 11]};
colsY = [5 12];
colsYMI = [5 7];
nBeh = length(labsBeh);

powX = 1;

% define random effects
colsRand = [9 10];


for xx = 1:length(nBX)
    for yy = 1:length(nBY)
        % loop through behavior types
        for bb = 1:nBeh
            % remove NaNs
            useIdx = ~any(ismissing(tableSp(:,[colsX colsPred{bb} colsRand colsY(bb)])),2);
            useTbl = tableSp(useIdx,:);

            % scale X variable
            useTbl{:,fNamesTbl{colsX}} = real(useTbl{:,fNamesTbl{colsX}}.^powX);


            %% Calculate mutual information

            % calculate mutual information
            varMice = useTbl.mouse;
            varSess = useTbl.session;
            varMiceIDs = unique(varMice);
            varSessIDs = unique(varSess);
            varX = useTbl{:,colsX};
            varY = useTbl{:,colsYMI(bb)};

            % calculate MI per session
            sMI = zeros(length(varMiceIDs),length(varSessIDs));
            permMI = zeros(length(varMiceIDs),length(varSessIDs),nPerm);
            for mm = 1:length(varMiceIDs)
                for ss = 1:length(varSessIDs)
                    % calculate MI
                    useIdx = find(varMice==varMiceIDs(mm) & varSess==varSessIDs(ss));
                    sMI(mm,ss) = mutInfoShannon(varX(useIdx),varY(useIdx),nBX(xx),nBY(yy));

                    % perform permutation
                    for ii = 1:nPerm
                        idxPerm = randperm(length(useIdx));
                        permMI(mm,ss,ii) = mutInfoShannon(varX(useIdx),varY(useIdx(idxPerm)),nBX(xx),nBY(yy));
                    end

                end
            end

            % calculate mean MI
            sessMI = mean(sMI,'all');

            % calculate MI significance
            permMIP = mean(sessMI <= squeeze(mean(permMI,[1,2])));

            % store results
            dataMI.(labsBeh{bb}).sessMI{xx,yy} = sMI;
            dataMI.(labsBeh{bb}).permMI{xx,yy} = permMI;
            dataMI.(labsBeh{bb}).permMIP(xx,yy) = permMIP;
        end
    end
end

% save('RevisionAnalysis\figures\MI_map.mat','dataMI')


%% Plot significance heatmap

figure; tiledlayout(1,nBeh)
for bb = 1:nBeh
    % load p-values
    curP = dataMI.(labsBeh{bb}).permMIP;

    % transform p-values
    plotP = 1./((1-curP)*nPerm);

    % plot heatmap
    nexttile(bb); hold on
    imagesc(plotP)

    % set labeles
    xlim([0.5 length(nBX)+0.5])
    ylim([0.5 length(nBY)+0.5])
    xticks(1:length(nBX))
    yticks(1:length(nBY))
    xticklabels(nBX)
    yticklabels(nBY)
    clim([0.0001 1])
    set(gca, 'ColorScale', 'log')
    colorbar
end



