%% reorganize data by stim conditions

function dataStruct = reorganDataStimCondition(dataStruct, stimInfo, chList)
    noCh = numel(chList);

    PulseDur=unique(stimInfo(:,1));
    PulseFreq=unique(stimInfo(:,2));
    
    noDur = length(PulseDur);
    noFreq = length(PulseFreq);
    noStim = noDur * noFreq;

    dataStruct.Avg = zeros(noStim,3);
    dataStruct.Ste = zeros(noStim,3);
    stimIdx = 0;
    for durIdx = 1:noDur
        for freqIdx = 1:noFreq
            stimIdx = stimIdx + 1;
            stimCondition=['stim_',num2str(PulseDur(durIdx)*1000),'ms_',num2str(PulseFreq(freqIdx)),'Hz'];
            % tmpStimIdx = find(stimInfo(:,1) == PulseDur(durIdx) & stimInfo(:,2) == PulseFreq(freqIdx));
            tmpStimIdx = stimInfo(:,1) == PulseDur(durIdx) & stimInfo(:,2) == PulseFreq(freqIdx);
            dataStruct.(stimCondition) = [];
            for chIdx = 1:noCh
                chName = ['ch',num2str(chList(chIdx))];
                
                tmpData = dataStruct.(chName)(tmpStimIdx,:);
                tmpChArr = zeros(length(tmpData),1) + chList(chIdx);
                dataStruct.(stimCondition) = [dataStruct.(stimCondition);tmpChArr, tmpData];
            end
            
            dataStruct.Avg(stimIdx,:) = mean(dataStruct.(stimCondition)(:,2:4));
            dataStruct.Ste(stimIdx,:) = std(dataStruct.(stimCondition)(:,2:4)) / sqrt(length(dataStruct.(stimCondition)));
        end
    end
end