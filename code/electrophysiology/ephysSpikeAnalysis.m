%%%%% Open ephys data analysis %%%%%%
% First run "load_open_ephys_binary.m" to extract data for each recording,
% then modify and run "Extract_and_save_20211021.m" to organize the data format 
% then load the data and stimulus timestmaps

%% ephys data analysis
% - spike/LFP filtering and segment trances by stimulation sessions
% - spike detection
% - spike rate
% - PSTH (spikes/LFP)
% by June Hoan Kim, updated on 12/8/2025
%
% 

close all; clear;

fpath='/Users/kimj65/Documents/Data/MEC_optogenetics/AAV8-hsyn_ChR2/20211110';
outpath = [fpath,'/selectedCh4'];
dataList = dir([fpath,'/Record1*.mat']);
noData = length(dataList);

for dataIdx = 1:noData
    close all;
    clearvars -except fpath outpath noData dataList dataIdx;
%% laod data
% fpath='/Users/kimj65/Documents/Data/MEC_optogenetics/AAV8-hsyn_ChR2/20211110';
% fname='Record1_20211110.mat';
fname = dataList(dataIdx).name
load(fullfile(fpath,fname));
variables = who('R*');

tmp = split(fname,'.');
prefix = tmp{1};
tmp = split(prefix,'_');
dataName = [tmp{1},'-',tmp{2}];

stim = eval(variables{2}).Timestamps;
voltage = eval(variables{1}).Data;
time = eval(variables{1}).Timestamps; 

sampleRate = 30000;
stimBlockPeriodSec = 1; % time period of 1 stimulation sequence (sec)
stablizePeriodSec = 1; % time after stim for washing out the effect of stim
stimBlockPeriod = stimBlockPeriodSec * sampleRate;
stablizePeriod = stablizePeriodSec * sampleRate;

totalNoCh = size(voltage,1);
chList = 1:totalNoCh;
% totalNoCh = length(chList);

noGoodCh = 5; % number of channel to analyze
goodChThreSNR = 4; % good channel: average peak amplitude / std > this value
adjShift = 0; % adjust mismatch between time and voltage sequences after bandpass filtering (unit: samples)

thresFold = 3; % threshold = avg + std * threshold Fold
threshold_noise = 9; % Peaks above this value are regarded as noise

tmp=split(fname,'.');
dataFileName=tmp{1};

spikeData=struct();
spikeData.Data=dataFileName;

%% segment data

segTimeRange = [stimBlockPeriodSec, stablizePeriodSec + stimBlockPeriodSec*5]; % [a, b]: -(t1-a) - (t2 + b) seconds @ t1 = stim on time and t2 = stim off time
spikeData=segEphysData(stim, voltage, time, segTimeRange, chList, sampleRate, adjShift);

%% peak detection for low freq activities (0.7 - 500 Hz)

spikeData.LFPThreshold = cell(totalNoCh, 1);
spikeData.LFPTimeAmp = cell(totalNoCh, 1);
% spikeData.SNR = zeros(totalNoCh, 36); % noStim 36
for chIdx = chList
    % tmpData = spikeData.LFPSeg{chList(chIdx)}(:,1:stimBlockPeriod)';
    tmpData = lowpass(spikeData.LFPSeg{chList(chIdx)}(:,1:stimBlockPeriod)', 50, sampleRate);
    tmpAvg = mean(tmpData)';
    tmpStd = std(tmpData)';
    spikeData.LFPThreshold{chList(chIdx)} = tmpAvg - tmpStd * thresFold;
    noiseThreshold = tmpAvg - tmpStd * threshold_noise;

    noStim = size(spikeData.stimuli_info,1);
    spikeData.LFPTimeAmp{chList(chIdx)} = cell(noStim,1);
    for stimIdx = 1:noStim
        % tmpData = spikeData.SpikeSeg{chList(chIdx)}(stimIdx,stimBlockPeriod:stimBlockPeriod+stimBlockPeriod);
        tmpData = spikeData.LFPSeg{chList(chIdx)}(stimIdx,:);

        % prominence to remove small peaks
        k_prom = 3;                   % pick between 1.0 - 3.0
        minProm = k_prom * tmpStd(stimIdx);
        
        [pks, locs] = findpeaks(-tmpData, 'MinPeakProminence', minProm,'MinPeakDistance',400);
        pks = - pks;
        peakIdx = find(pks < spikeData.LFPThreshold{chList(chIdx)}(stimIdx) & locs >= stimBlockPeriod  & locs <= stimBlockPeriod * 2);% & pks > noiseThreshold(stimIdx));
        
        if (numel(peakIdx) > 0)
            spikeData.LFPTimeAmp{chList(chIdx)}{stimIdx} = [locs(peakIdx)', pks(peakIdx)'];
            spikeData.LFPSNR(chList(chIdx), stimIdx) = -mean(spikeData.LFPTimeAmp{chList(chIdx)}{stimIdx}(:,2))/tmpStd(stimIdx); % SNR = average peak amplitudes / std
        end

        % plot(tmpData)
        % hold on;
        % plot(locs(peakIdx),pks(peakIdx),'r+');
        % hold off;
    end
end

%% peak detection for spikes

spikeData.spkThreshold = cell(totalNoCh, 1);
spikeData.SpkTimeAmp = cell(totalNoCh, 1);
spikeData.SNR = zeros(totalNoCh, 36); % noStim 36
for chIdx = chList
    tmpData = spikeData.SpikeSeg{chList(chIdx)}(:,1:stimBlockPeriod)';
    tmpAvg = mean(tmpData)';
    tmpStd = std(tmpData)';
    spikeData.spkThreshold{chList(chIdx)} = tmpAvg - tmpStd * thresFold;
    noiseThreshold = tmpAvg - tmpStd * threshold_noise;

    noStim = size(spikeData.stimuli_info,1);
    spikeData.SpkTimeAmp{chList(chIdx)} = cell(noStim,1);
    for stimIdx = 1:noStim
        % tmpData = spikeData.SpikeSeg{chList(chIdx)}(stimIdx,stimBlockPeriod:stimBlockPeriod+stimBlockPeriod);
        tmpData = spikeData.SpikeSeg{chList(chIdx)}(stimIdx,:);
        
        [pks, locs] = findpeaks(-tmpData);
        pks = - pks;
        spkIdx = find(pks < spikeData.spkThreshold{chList(chIdx)}(stimIdx) & pks > noiseThreshold(stimIdx));
        
        spikeData.SpkTimeAmp{chList(chIdx)}{stimIdx} = [locs(spkIdx)', pks(spkIdx)'];
        spikeData.SNR(chList(chIdx), stimIdx) = -mean(spikeData.SpkTimeAmp{chList(chIdx)}{stimIdx}(:,2))/tmpStd(stimIdx); % SNR = average peak amplitudes / std
    end
end

%% channel selection 
% avgSNR = mean(spikeData.SNR, 2);
% h=histogram(avgSNR,'NumBins',6);
% histoFig = gcf;
% ax=gca;
% ax.Title.String = dataFileName;
% ax.YLabel.String = 'Counts';
% ax.XLabel.String = 'mean SNR of channels';
% saveas(histoFig, [outpath,'/histo_avgSNR_',dataFileName,'.fig']);
% 
% [maxSNR, goodChList] = maxk(avgSNR, noGoodCh)
% spikeData.ChannelList = goodChList;
% % [goodChList, ~] = find(avgSNR > goodChThreSNR)
% % noGoodCh = numel(goodChList);
% % if(noGoodCh < 1)
% %     disp("there's no good channel!");
% %     return;
% % end

goodChList=[6,14,24,28];
noGoodCh = numel(goodChList);

%% spike count in period of light off in stim (stimNoLight) and pre/post stim (preStim/postStim) + PSTH (%%)

preStim = [1:stimBlockPeriod];
durStim = stimBlockPeriod + preStim;
postStim = stimBlockPeriod + stablizePeriod + durStim;

spikeData.SpikeRate = struct();
spikeData.SpikeRate.unit = 'Hz';
spikeData.SpikeRate.Data = 'preStim (1s) | durStim | postStim (1s) - 1s after last pulse pulse';

spikeData.PSTH = struct();%%
spikeData.PSTH.spike = struct();%%
spikeData.PSTH.LFP = struct();%%
for stimIdx = 1:noStim
    info = [num2str(spikeData.stimuli_info(stimIdx,1)*1000),'ms_' , num2str(spikeData.stimuli_info(stimIdx,2)),'Hz'];%%
    spikeData.PSTH.spike.(['stim',num2str(stimIdx),'_',info]) = cell(noGoodCh,1);%%
    spikeData.PSTH.LFP.(['stim',num2str(stimIdx),'_',info]) = cell(noGoodCh,1);%%

    % pulse timing
    x1 = double([spikeData.stimON{stimIdx}]-[spikeData.stimON{stimIdx}(1)])+stimBlockPeriod;%% % pulse on times
    x2 = double([spikeData.stimOFF{stimIdx}]-[spikeData.stimON{stimIdx}(1)])+stimBlockPeriod;%% % pulse off times
    pulInterval = min(unique(diff(x1)));%%
    
    noPulse = length(spikeData.stimON{stimIdx});
    tmpShift = - spikeData.stimON{stimIdx}(1) + min(durStim);

    for chIdx = 1:noGoodCh
        tmpSpkBin = zeros(size(spikeData.SpikeSeg{goodChList(chIdx)},2),4);
        tmpSpkBin(spikeData.SpkTimeAmp{goodChList(chIdx)}{stimIdx}(:,1),1) = 1;
        tmpSpkBin(preStim,2) = 1;
        tmpSpkBin(postStim,4) = 1;

        tmpSpkBin(durStim, 3) = 1;
        
        spikeData.PSTH.spike.(['stim',num2str(stimIdx),'_',info]){chIdx} = zeros(pulInterval+1,1);%%
        spikeData.PSTH.LFP.(['stim',num2str(stimIdx),'_',info]){chIdx} = zeros(pulInterval+1,1);%%
        for pulseIdx = 1:noPulse%%
            tmpStim = [spikeData.stimON{stimIdx}(pulseIdx):spikeData.stimOFF{stimIdx}(pulseIdx)-1] + tmpShift;%%
            tmpSpkBin(tmpStim, 3) = 0;%%

            % tmpSpikeSeg = spikeData.SpikeSeg{goodChList(chIdx)}(stimIdx,[x1(pulseIdx):x1(pulseIdx)+pulInterval])';%%
            tmpSpikeSeg = tmpSpkBin([x1(pulseIdx):x1(pulseIdx)+pulInterval],1);%%
            tmpLFPSeg = spikeData.LFPSeg{goodChList(chIdx)}(stimIdx,[x1(pulseIdx):x1(pulseIdx)+pulInterval])';%%
            % tmprawSeg = spikeData.rawSeg{goodChList(chIdx)}(stimIdx,[x1(pulseIdx):x1(pulseIdx)+pulInterval])';%%

            % plot([x1(pulseIdx):x1(pulseIdx)+pulInterval],tmprawSeg);%%

            spikeData.PSTH.spike.(['stim',num2str(stimIdx),'_',info]){chIdx} = spikeData.PSTH.spike.(['stim',num2str(stimIdx),'_',info]){chIdx} + tmpSpikeSeg;%%
            spikeData.PSTH.LFP.(['stim',num2str(stimIdx),'_',info]){chIdx} = spikeData.PSTH.LFP.(['stim',num2str(stimIdx),'_',info]){chIdx} + tmpLFPSeg;%%
        end
        spikeData.PSTH.spike.(['stim',num2str(stimIdx),'_',info]){chIdx} = spikeData.PSTH.spike.(['stim',num2str(stimIdx),'_',info]){chIdx} / noPulse;%%
        spikeData.PSTH.LFP.(['stim',num2str(stimIdx),'_',info]){chIdx} = spikeData.PSTH.LFP.(['stim',num2str(stimIdx),'_',info]){chIdx} / noPulse;%%

        tmpPeriod = sum(tmpSpkBin);
        spkPreStim = sum(tmpSpkBin(:,1).*tmpSpkBin(:,2)) / (tmpPeriod(2) / sampleRate);
        spkDurStim = sum(tmpSpkBin(:,1).*tmpSpkBin(:,3)) / (tmpPeriod(3) / sampleRate);
        spkPostStim = sum(tmpSpkBin(:,1).*tmpSpkBin(:,4)) / (tmpPeriod(4) / sampleRate);
        spikeData.SpikeRate.(['ch',num2str(goodChList(chIdx))])(stimIdx,:) = [spkPreStim, spkDurStim, spkPostStim];
    end
end

%% reorganize data by sitm condition

spikeData.SpikeRate = reorganDataStimCondition(spikeData.SpikeRate, spikeData.stimuli_info, goodChList);

%% save data
% outFileName = [dataFileName,'_processed.mat'];
% save(fullfile(outpath,outFileName),'spikeData', '-v7.3');

%% plot firing rate change
% allFields = fieldnames(spikeData.SpikeRate);
% pattern = '10Hz$';
% selFields=allFields(~cellfun('isempty', regexp(allFields, pattern, 'once')));
% noSelStim = length(selFields);
% 
% selStim = [1:3:12];
% tmpAvg = spikeData.SpikeRate.Avg(selStim,:);
% tmpSte = spikeData.SpikeRate.Ste(selStim,:);
% fig = figure;
% for stimIdx = 1:noSelStim
%     errorbar([1:3],tmpAvg(stimIdx,:),tmpSte(stimIdx,:),'o-','LineWidth', 2,'MarkerSize',5);
%     hold on;
%     tmp = split(selFields{stimIdx},'_');
%     stimName{stimIdx} = [tmp{2},' ',tmp{3}];
% end
% hold off;
% ax = gca;
% ax.XLim = [0.7, 3.3];
% ax.YLim(1) = 0;
% ax.XTick = [1,2,3];
% ax.XTickLabel = {'Pre', 'Stim', 'Post'};
% ax.YLabel.String = 'Spike Rate (Hz)';
% legend(stimName);
% title(dataName);
% outFileName = [dataFileName,'_SpikeRate.fig'];
% saveas(fig,fullfile(outpath,outFileName));

end
%% "FOR TEST" plot filtered data, peaks, and light stimumation.

% stimIdx - 1: 20ms 10Hz, 9/12: 2ms 10Hz
% figure;
% subplot(2,1,1);
% ch=goodChList(1);
% stimIdx = 9;
% dataLength = length(tmpSpkBin);
% plot(spikeData.SpikeSeg{ch}(stimIdx, :),'-');
% hold on;
% plot(spikeData.SpkTimeAmp{ch}{stimIdx}(:,1),spikeData.SpkTimeAmp{ch}{stimIdx}(:,2),'*');
% noSpk = length(spikeData.SpkTimeAmp{ch}{stimIdx}(:,1));
% plot(spikeData.SpkTimeAmp{ch}{stimIdx}(:,1),zeros(noSpk,1)-25,'|');
% 
% % plot stim
% x1 = double([spikeData.stimON{stimIdx}]-[spikeData.stimON{stimIdx}(1)])+30000;
% x2 = double([spikeData.stimOFF{stimIdx}]-[spikeData.stimON{stimIdx}(1)])+30000;
% 
% for pulse_num = 1:length(x1)
%     x = [x1(pulse_num) x1(pulse_num) x2(pulse_num) x2(pulse_num)];
%     y = [-100 100 100 -100];
%     h = fill(x,y,'r','EdgeColor','none'); % for antibiased
%     h(1).FaceColor = [0,0.4,0.9]; 
%     h(1).FaceAlpha = 0.1;
%     hold on
% end
% plot([1:240000],tmpSpkBin(:,3)*-30,'.'); %%%%%%%%%%
% ax = gca;
% xlim([1,240000]);
% ylim([-100,100]);
% ax.XTick = [0:8]*30000;
% ax.XTickLabel = {'0','1','2','3','4','5','6','7','8'};
% xlabel('Time (s)');
% ylabel('voltage (a.u.)');
% title('ChR2 - 20211110(Rec1) 2ms-10Hz');