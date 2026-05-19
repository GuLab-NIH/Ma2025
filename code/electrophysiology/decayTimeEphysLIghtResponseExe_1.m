% Time constant decay time of light-evoked ePhys resopnses by repetitive light pulses
%
% by June Hoan Kim, December 30, 2025

close all; clear;

fpath='/Users/kimj65/Documents/Data/MEC_optogenetics/AAV8-hsyn_ChR2/20211110/5MaxSNR_Ch';
dataList = dir([fpath,'/Record1*processed.mat']);

noData = length(dataList);

for dataIdx = 1%:noData
    close all;

    fname = dataList(dataIdx).name;
    load(fullfile(fpath,fname));
    
    tmp = split(fname,'.');
    prefix = tmp{1};
    tmp = split(prefix,'_');
    dataName = [tmp{1},'-',tmp{2}];
    
    tmp=split(fname,'.');
    dataFileName=tmp{1}
    
    sampleRate = 30000;
    % stimBlockPeriodSec = 1; % time period of 1 stimulation sequence (sec)
    % stimBlockPeriod = stimBlockPeriodSec * sampleRate;
    
    % allFields = fieldnames(spikeData.SpikeRate);
    % pattern = '^ch';
    % tmpChList=allFields(~cellfun('isempty', regexp(allFields, pattern, 'once')));
    % totalNoCh = 32;
    chList = 24;%1:totalNoCh;%str2num(tmpChList{1}(3:end))%
    totalNoCh = length(chList);

    stimDur = [2, 5, 10, 20]/1000;
    stimFreq = [10, 20, 40];

    noStim = length(spikeData.stimuli_info);
    spikeData.LFPfit.timeConstant = struct();
    spikeData.LFPfit.adjRsq = struct();

    for chIdx = 1:totalNoCh
        chName = ['Ch',num2str(chList(chIdx))];

        close all;

        for i=1:3
            fig(i) = figure;
            t(i) = tiledlayout(fig(i), 3,4, "TileSpacing", "tight", "Padding", "compact"); % Create a 3x5 layout with tight spacing and compact padding
            title(t(i),[dataName,' ',chName,'-',num2str(i)],'fontSize', 20);
        end
        
        spikeData.LFPfit.timeConstant.(chName) = zeros(noStim,1);
        spikeData.LFPfit.adjRsq.(chName) = zeros(noStim,1);
        tmpTimeConstant = nan(noStim,1);
        tmpAdjRsq = nan(noStim,1);

        for freqNum = 1:length(stimFreq)
            for durNum = 1:length(stimDur)
                stimIdx = find(spikeData.stimuli_info(:,1) == stimDur(durNum) & spikeData.stimuli_info(:,2) == stimFreq(freqNum));
                stimInfo = [num2str(stimDur(durNum)*1000),'ms ' , num2str(stimFreq(freqNum)),'Hz'];
                
                for stimRound = 1%:length(stimIdx)
                    
                    % stimName = ['stim',num2str(stimDur(durNum)*1000),'ms_' , num2str(stimFreq(freqNum)),'Hz_',num2str(stimRound)];
                    figure(fig(stimRound));
                    nexttile(t(stimRound),(freqNum-1)*4+durNum);
                    ax = gca;
                    tmpRange = 1:90000;
                    tmpX = (tmpRange-1)/30000;
                    tmpSeg = spikeData.LFPSeg{chList(chIdx)}(stimIdx(stimRound),tmpRange);
                    plot(tmpX,tmpSeg);
                    hold on;

                    tmpData = spikeData.LFPTimeAmp{chList(chIdx)}{stimIdx(stimRound)};
                    

                    if(length(tmpData) < 5)
                        continue;
                    end
            
                    tmpData(:,1) = tmpData(:,1) / sampleRate;
                    fitResult = expGrowthFit(tmpData);
                    if(~isnan(fitResult.Params.k))
                        
                        h = plot(ax,fitResult.Params, tmpData(:,1), tmpData(:,2));  % let cfit/plot handle data & fit together
                        set(h(1), 'Marker', 'o', 'MarkerFaceColor', 'k'); % data
                        set(h(2), 'Color', 'r');%, 'LineWidth', 2);          % fitted curve
                    
                        xlabel('time (s)');
                        ylabel('uV');
                        % title('Exponential Growth Fit:  y = A + B e^{-kx}');
                        grid on;
                        
                        ax = gca;
                        ax.YLim = [-1000, 200];
                        
                        ax.XTickLabel='';
                        ax.YTickLabel='';
                        ax.YLabel.String='';
                        ax.XLabel.String='';
                        title(stimInfo);
            
                        tmpTimeConstant(stimIdx(stimRound)) = 1/fitResult.Params.k; % sec
                        tmpAdjRsq(stimIdx(stimRound)) = fitResult.gof.adjrsquare;
                        % legend([info, ' ({\tau}=', num2str(timeConstant), ' s)'], ['fit (adjR^2 = ',num2str(fitResult.gof.adjrsquare),')'], 'Location','southeast','FontSize',15);
                        lgd = legend(ax,[ax.Children(2),ax.Children(1)],num2str(tmpTimeConstant(stimIdx(stimRound))),num2str(fitResult.gof.adjrsquare), 'Location','south');
                        lgd.Box='off';
                        lgd.IconColumnWidth = 10;
                    end
                end
            end
        end
        % tmpIdx = find (tmpAdjRsq < 0.6);
        % tmpTimeConstant(tmpIdx) = nan;
    
        spikeData.LFPfit.timeConstant.(chName) = tmpTimeConstant;
        spikeData.LFPfit.adjRsq.(chName) = tmpAdjRsq;
        for stimRound = 1:length(stimIdx)
            % saveas(fig(stimRound),fullfile(fpath,[dataName,chName,'-',num2str(stimRound),'_LFP_peaksFit.fig']));
        end
    end    
    % save(fullfile(dataList(dataIdx).folder,[dataName,'_processed_1.mat']),'spikeData', '-v7.3');
end