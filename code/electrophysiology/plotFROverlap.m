close all; clear;

list=dir('/Users/kimj65/Documents/Data/MEC_optogenetics/AAV8-hsyn_GFP/20211111/5MaxSNR_Ch/Record1*_processed.mat');
noList = length(list);


% chList = [1,2,3,4,9,11,12,17,21,22,23,24,31]; %Rec1-2021110
% chList = [1,8,11,12,13,14,21,22,23,24,28,30,31,32];% %Rec1-2021116
chList = 1:32;
noCh = length(chList);

sampleRate = 30000;
stimBlockPeriodSec = 1; % time period of 1 stimulation sequence (sec)
stablizePeriodSec = 1; % time after stim for washing out the effect of stim 2s for GFP 1s for ChR2
stimBlockPeriod = stimBlockPeriodSec * sampleRate;
stablizePeriod = stablizePeriodSec * sampleRate;

exZoneIdx = 224;
exZoneTime2 = 80; % (ms) exclude baseline for caculation of firing rate 80ms means using first 20 ms
exZoneIdx2 = exZoneTime2 / 1000 * sampleRate;

preStim = [1:stimBlockPeriod];
durStim = stimBlockPeriod + preStim;
postStim = stimBlockPeriod + stablizePeriod + durStim;

preStim = preStim(1:end-150); %% for ChR2

avgData = [];
for i = 1:noList
    load(fullfile(list(i).folder,list(i).name));

    tmp=split(list(i).name,'.');
    tmp = split(tmp{1}, '_');
    dataFileName=[tmp{1},'_',tmp{2}]

    stimConditions = fieldnames(spikeData.PSTH.spike);
    stimIdx = find(spikeData.stimuli_info(:,1) == 0.002 & spikeData.stimuli_info(:,2) == 10);

    spikeData.Count.avgSpikeRate = struct();
    % spikeData.Count.maxSpikeRate = struct();
    mergedFR=[];
    mergedPreFR=[];
    mergedPostFR=[];
    figure;
    for stimNum = 1%:length(stimIdx)
        spikeData.Count.spike.(stimConditions{stimIdx(stimNum)}) = cell(noCh,1);%%
    
        % pulse timing
        x1 = double([spikeData.stimON{stimIdx(stimNum)}]-[spikeData.stimON{stimIdx(stimNum)}(1)])+stimBlockPeriod;%% % pulse on times
        x2 = double([spikeData.stimOFF{stimIdx(stimNum)}]-[spikeData.stimON{stimIdx(stimNum)}(1)])+stimBlockPeriod;%% % pulse off times
        pulInterval = min(unique(diff(x1)));%%
        
        noPulse = length(spikeData.stimON{stimIdx(stimNum)});
        tmpShift = - spikeData.stimON{stimIdx(stimNum)}(1) + min(durStim);
        
        for chIdx = 1:noCh
            spikeData.Count.avgSpikeRate.(['ch',num2str(chList(chIdx)),'_',stimConditions{stimIdx(stimNum)}]) = zeros(1,3);

            % spikeData.Count.maxSpikeRate.(['ch',num2str(chList(chIdx)),'_',stimConditions{stimIdx(stimNum)}]) = zeros(3,3);

            tmpSpkBin = zeros(size(spikeData.SpikeSeg{chList(chIdx)},2),5);
            tmpSpkBin(spikeData.SpkTimeAmp{chList(chIdx)}{stimIdx(stimNum)}(:,1),1) = 1; % spikes
            tmpSpkBin(preStim,2) = 1; % pre stim area
            tmpSpkBin(postStim,4) = 1; % post stim area
    
            tmpSpkBin(durStim, 3) = 1; % stim area
            tmpSpkBin(durStim, 5) = 1; % stim area with artivact exclusion
            
            spikeData.Count.spike.(stimConditions{stimIdx(stimNum)}){chIdx} = zeros(pulInterval+1,1);%%
            preCount = zeros(pulInterval+1,1);%%
            postCount = zeros(pulInterval+1,1);%%
            for pulseIdx = 1:noPulse%%
                tmpStim1 = [spikeData.stimON{stimIdx(stimNum)}(pulseIdx):spikeData.stimOFF{stimIdx(stimNum)}(pulseIdx)-1] + tmpShift;%%
                tmpSpkBin(tmpStim1, 3) = 0;%%
                tmpSpkBin(tmpStim1, 5) = 0;%%

                tmpStim2 = [spikeData.stimON{stimIdx(stimNum)}(pulseIdx)-exZoneIdx2:spikeData.stimON{stimIdx(stimNum)}(pulseIdx)+exZoneIdx] + tmpShift;%%
                tmpSpkBin(tmpStim2, 5) = 0;%%
    
                % tmpSpikeSeg = spikeData.SpikeSeg{goodChList(chIdx)}(stimIdx(stimNum),[x1(pulseIdx):x1(pulseIdx)+pulInterval])';%%
                tmpSpikeSeg = tmpSpkBin([x1(pulseIdx):x1(pulseIdx)+pulInterval],1);%%
                % tmprawSeg = spikeData.rawSeg{goodChList(chIdx)}(stimIdx(stimNum),[x1(pulseIdx):x1(pulseIdx)+pulInterval])';%%
    
                % plot([x1(pulseIdx):x1(pulseIdx)+pulInterval],tmprawSeg);%%
    
                spikeData.Count.spike.(stimConditions{stimIdx(stimNum)}){chIdx} = spikeData.Count.spike.(stimConditions{stimIdx(stimNum)}){chIdx} + tmpSpikeSeg;%%
            end
            spikeData.Count.spike.(stimConditions{stimIdx(stimNum)}){chIdx} = spikeData.Count.spike.(stimConditions{stimIdx(stimNum)}){chIdx} / noPulse;%%
            preCount = tmpSpkBin(preStim,1);%%
            postCount = tmpSpkBin(postStim,1);%%
            spkPreStim = tmpSpkBin(:,1).*tmpSpkBin(:,2);
            spkDurStim = tmpSpkBin(:,1).*tmpSpkBin(:,3);
            spkDurStim2 = tmpSpkBin(:,1).*tmpSpkBin(:,5);
            spkPostStim = tmpSpkBin(:,1).*tmpSpkBin(:,4);

            tmpDurFR = spkFR(spikeData.Count.spike.(stimConditions{stimIdx(stimNum)}){chIdx}, 0.005, sampleRate);
            tmpPreFR = spkFR(preCount, 0.005, sampleRate);
            tmpPostFR = spkFR(postCount,0.005, sampleRate);

            % spkDurStim = tmpSpkBin(414:end-15,1).*tmpSpkBin(414:end-15,3); % why 414:end-15 ????
            spkDurStim = tmpSpkBin(:,1).*tmpSpkBin(:,3);
    
            tmpPeriod = sum(tmpSpkBin);
            meanSpkPreStim = sum(spkPreStim) / (tmpPeriod(2) / sampleRate);
            meanSpkDurStim = sum(spkDurStim2) / (tmpPeriod(5) / sampleRate);
            meanSpkPostStim = sum(spkPostStim) / (tmpPeriod(4) / sampleRate);
            spikeData.Count.avgSpikeRate.(['ch',num2str(chList(chIdx)),'_',stimConditions{stimIdx(stimNum)}]) = [meanSpkPreStim, meanSpkDurStim, meanSpkPostStim];

            % maxSpkPreStim = max(tmpPreFR(:,2));
            % maxSpkDurStim = max(tmpDurFR(:,2));
            % maxSpkPostStim = max(tmpPostFR(:,2));
            % spikeData.Count.MaxSpikeRate.(['ch',num2str(chList(chIdx)),'_',stimConditions{stimIdx(stimNum)}])(stimNum,:) = [maxSpkPreStim, maxSpkDurStim, maxSpkPostStim];

            mergedFR = [mergedFR;tmpDurFR(:,2)'];
            mergedPreFR = [mergedPreFR;tmpPreFR(:,2)'];
            mergedPostFR = [mergedPostFR;tmpPostFR(:,2)'];
            plot(tmpDurFR(1:end-150,2),'-','Color','#AAAAAA');
            hold on;
        end
    end
    meanFR = mean(mergedFR);
    meanPreFR = mean(mergedPreFR);
    meanPostFR = mean(mergedPostFR);
    plot(meanPreFR(1:pulInterval-150),'g-','LineWidth',2);
    plot(meanPostFR(1+pulInterval*2:pulInterval+pulInterval*2-150),'b-','LineWidth',2);

    plot(meanFR(1:end-150),'k-','LineWidth',2);
    
    % finding light-artifact exclusive time zone from GFP
    % totMeanFR = mean(meanFR(500:end));
    % totStdFR = std(meanFR(500:end));
    % threshold = totMeanFR + totStdFR * 4.152;
    % noPts = numel(mergedFR(:,500:end));
    % noOverThrePts = nnz(mergedFR(:,500:end) > threshold);
    % ratioOverThrePts = noOverThrePts / noPts
    % [~, tmpIdx] = find (meanFR > threshold);
    % exZoneIdx = max(tmpIdx)
    % plot([0,3000],[threshold, threshold],'r-', 'LineWidth', 2);
    h = fill([1 1 exZoneIdx exZoneIdx],[0 1600 1600 0],'r','EdgeColor','none');
    h(1).FaceColor = '#FF0000'; 
    h(1).FaceAlpha = 0.1;
    ax = gca;
    ax.XTick = [0:600:3000];
    ax.XTickLabel = {'0', '20', '40', '60', '80', '100'};
    ax.XLabel.String = 'Time (msec)';
    ax.YLabel.String = 'Firing rate (Hz)';
    % pp = findobj('type', 'line');
    % lgd = legend(pp(1), ['avg + (std x 4.152) = ', num2str(exZoneIdx / 30), 'ms']);
    hold off;
    % fig=gcf;
    % outFileName = [dataFileName,'_SpkFROverlap.fig'];
    % saveas(fig,fullfile(list(i).folder,outFileName));
    % 
    % save(fullfile(list(i).folder,list(i).name),'spikeData', '-v7.3');

    data = fieldnames(spikeData.Count.avgSpikeRate);
    noData = length(data);
    for i=1:noData
        avgData = [avgData; spikeData.Count.avgSpikeRate.(data{i})];
    end
end

% plot avg FR
% data = fieldnames(spikeData.Count.avgSpikeRate);
% noData = length(data);
% avgData=nan(noData,3);

noData = length(avgData);
figure;
for i=1:noData
    % avgData(i,:) = spikeData.Count.avgSpikeRate.(data{i});
    % plot([1:3],spikeData.Count.avgSpikeRate.(data{i}),'.-');
    plot([1:3],avgData(i,:),'.-','Color',[0.5, 0.5, 0.5]);
    hold on;
end
totavgData = mean(avgData,'omitnan');
totsteData = std(avgData,'omitnan')/sqrt(length(avgData));
% plot([1:3],totavgData,'k.-','LineWidth',3,'MarkerSize',10);
err=errorbar([1:3],totavgData,totsteData,'Marker','.','Color','black','MarkerSize',20,'LineWidth',2,'CapSize',5);
ax=gca;
ax.XTick=[1:3];
ax.XLim = [0.7,3.3];
ax.XTickLabel = {'Pre','Stim','Post'};
ax.YLabel.String='mean FR (Hz)';

[p,tbl,stats] = friedman(avgData);
figure;
[c,m,h,gnames] = multcompare(stats);
disp(c);