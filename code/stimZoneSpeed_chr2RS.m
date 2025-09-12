

%% chr2 mice only
bin=1; %cm
beforeAfterDis=15;%15cm before and after speed
allSpeedStim={};
allSpeedBefore={};
allSpeedAfter={};

% for mouse=1:length(stimZoneLengthAllData);
for mouse=1:15;
    disp(mouse)
    allSpeedStim{mouse}={};
    allSpeedBefore{mouse}={};
    allSpeedAfter{mouse}={};
    for run=1:length(stimZoneLengthAllData{mouse});
        disp(run)
          allSpeedStim{mouse}{run}={};
    allSpeedBefore{mouse}{run}={};
    allSpeedAfter{mouse}{run}={};
        if ~isempty(stimZoneLengthAllData{mouse}{run});%stimulation should be at least 14cm
            stimIdx=stimZoneIdxAllDataFix{mouse}{run};
            position=positionSegAllData{mouse}{run};
            time=timeSegAllData{mouse}{run};
            for nStim=1:size(stimIdx,1)
                allSpeedStim{mouse}{run}{nStim}=[];
    allSpeedBefore{mouse}{run}{nStim}=[];
    allSpeedAfter{mouse}{run}{nStim}=[];
    thisStimIdx=[stimIdx(nStim,1):1:stimIdx(nStim,2)];
    thisPosition=position(thisStimIdx)-position(thisStimIdx(1));


    thisTime=time(thisStimIdx)-time(thisStimIdx(1));

    %during stim speed
    positionBins=[0:bin:ceil(thisPosition(end))];

    % positionBins=[0:bin:15];
    [~,nbin]=histc(thisPosition,positionBins);


    for binNumber=1:max(nbin);
        positionThisBin=thisPosition(nbin==binNumber);
        timeThisBin=thisTime(nbin==binNumber);
        if max(timeThisBin)>0
        allSpeedStim{mouse}{run}{nStim}(binNumber)=(positionThisBin(end)-positionThisBin(1))/((timeThisBin(end)-timeThisBin(1))*60); %the original time already in min 
        end
    end

    %calculating before sim speed
    stimStartPosition=position(thisStimIdx(1));%no subtraction by first one
    beforePosition=stimStartPosition-beforeAfterDis;
    iInclude=intersect(find(position>=beforePosition),find(position<stimStartPosition));

    thisPositionBefore=position(iInclude);
    thisPositionBefore=thisPositionBefore-thisPositionBefore(1);
    thisTimeBefore=time(iInclude);
    thisTimeBefore=thisTimeBefore-thisTimeBefore(1);

    % positionBinsBefore=[0:bin:ceil(thisPositionBefore(end))];
    positionBinsBefore=[0:bin:beforeAfterDis];
    [~,nbinBefore]=histc(thisPositionBefore,positionBinsBefore);

     for binNumber=1:max(nbinBefore);
        positionThisBinBefore=thisPositionBefore(nbinBefore==binNumber);
        timeThisBinBefore=thisTimeBefore(nbinBefore==binNumber);
        if max(timeThisBinBefore)>0
        allSpeedBefore{mouse}{run}{nStim}(binNumber)=(positionThisBinBefore(end)-positionThisBinBefore(1))/((timeThisBinBefore(end)-timeThisBinBefore(1))*60); %the original time already in min 
        end
     end

     %calculating after stim speed

 stimEndPosition=position(thisStimIdx(end));%no subtraction by first one
    afterPosition=stimEndPosition+beforeAfterDis;
    iInclude=intersect(find(position<=afterPosition),find(position>stimEndPosition));

    thisPositionAfter=position(iInclude);
    thisPositionAfter=thisPositionAfter-thisPositionAfter(1);
    thisTimeAfter=time(iInclude);
    thisTimeAfter=thisTimeAfter-thisTimeAfter(1);
    % 
    % positionBinsAfter=[0:bin:ceil(thisPositionAfter(end))];
    positionBinsAfter=[0:bin:beforeAfterDis];
    [~,nbinAfter]=histc(thisPositionAfter,positionBinsAfter);

     for binNumber=1:max(nbinAfter);
        positionThisBinAfter=thisPositionAfter(nbinAfter==binNumber);
        timeThisBinAfter=thisTimeAfter(nbinAfter==binNumber);
        if max(timeThisBinAfter)>0
        allSpeedAfter{mouse}{run}{nStim}(binNumber)=(positionThisBinAfter(end)-positionThisBinAfter(1))/((timeThisBinAfter(end)-timeThisBinAfter(1))*60); %the original time already in min 
        end
     end

            end
        else
          allSpeedAfter{mouse}{run}=[];
          allSpeedBefore{mouse}{run}=[];
          allSpeedStim{mouse}{run}=[];
        end
    end
end


% save('allSpeedAfter.mat','allSpeedAfter')
% save('allSpeedBefore.mat','allSpeedBefore')
% save('allSpeedStim.mat','allSpeedStim')

%% merge all mice runs together

allSpeedStimAll=[];
allSpeedAfterAll=[];
allSpeedBeforeAll=[];


for mouse=1:length(allSpeedStim)
    for run=1:length(allSpeedStim{mouse});
         if ~isempty(allSpeedStim{mouse}{run});
        for stim=1:length(allSpeedStim{mouse}{run});
            if length(allSpeedStim{mouse}{run}{stim})>=15/bin;
            allSpeedStimAll=[allSpeedStimAll;allSpeedStim{mouse}{run}{stim}(1:15/bin)];
            
            if length(allSpeedBefore{mouse}{run}{stim})>=15/bin;
            allSpeedBeforeAll=[allSpeedBeforeAll;allSpeedBefore{mouse}{run}{stim}(end-15/bin+1:end)];
            end
            if length(allSpeedAfter{mouse}{run}{stim})>=15/bin;
            allSpeedAfterAll=[allSpeedAfterAll;allSpeedAfter{mouse}{run}{stim}(1:15/bin)];
            end
            end
        end
        end
     end
end

%% 

figure
x=[0:bin:15]+bin/2;
x=x(x<=15);

errorbar(x,nanmean(allSpeedBeforeAll,1),nansem(allSpeedBeforeAll,1),'k');
hold on
errorbar(15+x,nanmean(allSpeedStimAll,1),nansem(allSpeedStimAll,1),'b');
hold on
errorbar(30+x,nanmean(allSpeedAfterAll,1),nansem(allSpeedAfterAll,1),'k');
title(['beforeStimAfter',num2str(bin),'cm']);

saveas(gcf,'chr2RS_beforeStimAfter1cm.fig')

%  %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %%
% bin=2.5; %cm
% beforeAfterDis=15;%15cm before and after speed
% allSpeedStim={};
% allSpeedBefore={};
% allSpeedAfter={};
% for mouse=1:length(stimZoneLengthAllData);
%     disp(mouse)
%     allSpeedStim{mouse}={};
%     allSpeedBefore{mouse}={};
%     allSpeedAfter{mouse}={};
%     for run=1:length(stimZoneLengthAllData{mouse});
%         disp(run)
%           allSpeedStim{mouse}{run}={};
%     allSpeedBefore{mouse}{run}={};
%     allSpeedAfter{mouse}{run}={};
%         if ~isempty(stimZoneLengthAllData{mouse}{run});%stimulation should be at least 14cm
%             stimIdx=stimZoneIdxAllDataFix{mouse}{run};
%             position=positionSegAllData{mouse}{run};
%             time=timeSegAllData{mouse}{run};
%             for nStim=1:size(stimIdx,1)
%                 allSpeedStim{mouse}{run}{nStim}=[];
%     allSpeedBefore{mouse}{run}{nStim}=[];
%     allSpeedAfter{mouse}{run}{nStim}=[];
%     thisStimIdx=[stimIdx(nStim,1):1:stimIdx(nStim,2)];
%     thisPosition=position(thisStimIdx)-position(thisStimIdx(1));
%     thisTime=time(thisStimIdx)-time(thisStimIdx(1));
% 
%     %during stim speed
%     positionBins=[0:bin:ceil(thisPosition(end))];
%     % positionBins=[0:bin:15];
%     [~,nbin]=histc(thisPosition,positionBins);
% 
% 
%     for binNumber=1:max(nbin);
%         positionThisBin=thisPosition(nbin==binNumber);
%         timeThisBin=thisTime(nbin==binNumber);
%         if max(timeThisBin)>0
%         allSpeedStim{mouse}{run}{nStim}(binNumber)=(positionThisBin(end)-positionThisBin(1))/((timeThisBin(end)-timeThisBin(1))*60); %the original time already in min 
%         end
%     end
% 
%     %calculating before sim speed
%     stimStartPosition=position(thisStimIdx(1));%no subtraction by first one
%     beforePosition=stimStartPosition-beforeAfterDis;
%     iInclude=intersect(find(position>=beforePosition),find(position<stimStartPosition));
% 
%     thisPositionBefore=position(iInclude);
%     thisPositionBefore=thisPositionBefore-thisPositionBefore(1);
%     thisTimeBefore=time(iInclude);
%     thisTimeBefore=thisTimeBefore-thisTimeBefore(1);
% 
%     % positionBinsBefore=[0:bin:ceil(thisPositionBefore(end))];
%     positionBinsBefore=[0:bin:beforeAfterDis];
%     [~,nbinBefore]=histc(thisPositionBefore,positionBinsBefore);
% 
%      for binNumber=1:max(nbinBefore);
%         positionThisBinBefore=thisPositionBefore(nbinBefore==binNumber);
%         timeThisBinBefore=thisTimeBefore(nbinBefore==binNumber);
%         if max(timeThisBinBefore)>0
%         allSpeedBefore{mouse}{run}{nStim}(binNumber)=(positionThisBinBefore(end)-positionThisBinBefore(1))/((timeThisBinBefore(end)-timeThisBinBefore(1))*60); %the original time already in min 
%         end
%      end
% 
%      %calculating after stim speed
% 
%  stimEndPosition=position(thisStimIdx(end));%no subtraction by first one
%     afterPosition=stimEndPosition+beforeAfterDis;
%     iInclude=intersect(find(position<=afterPosition),find(position>stimEndPosition));
% 
%     thisPositionAfter=position(iInclude);
%     thisPositionAfter=thisPositionAfter-thisPositionAfter(1);
%     thisTimeAfter=time(iInclude);
%     thisTimeAfter=thisTimeAfter-thisTimeAfter(1);
%     % 
%     % positionBinsAfter=[0:bin:ceil(thisPositionAfter(end))];
%     positionBinsAfter=[0:bin:beforeAfterDis];
%     [~,nbinAfter]=histc(thisPositionAfter,positionBinsAfter);
% 
%      for binNumber=1:max(nbinAfter);
%         positionThisBinAfter=thisPositionAfter(nbinAfter==binNumber);
%         timeThisBinAfter=thisTimeAfter(nbinAfter==binNumber);
%         if max(timeThisBinAfter)>0
%         allSpeedAfter{mouse}{run}{nStim}(binNumber)=(positionThisBinAfter(end)-positionThisBinAfter(1))/((timeThisBinAfter(end)-timeThisBinAfter(1))*60); %the original time already in min 
%         end
%      end
%             end
%         else
%           allSpeedAfter{mouse}{run}=[];
%           allSpeedBefore{mouse}{run}=[];
%           allSpeedStim{mouse}{run}=[];
%         end
%     end
% end
% 
% % 
% % save('allSpeedAfter.mat','allSpeedAfter')
% % save('allSpeedBefore.mat','allSpeedBefore')
% % save('allSpeedStim.mat','allSpeedStim')
% 
% %% merge all mice runs together
% 
% allSpeedStimAll=[];
% allSpeedAfterAll=[];
% allSpeedBeforeAll=[];
% 
% 
% for mouse=1:length(allSpeedStim)
%     for run=1:length(allSpeedStim{mouse});
%          if ~isempty(allSpeedStim{mouse}{run});
%         for stim=1:length(allSpeedStim{mouse}{run});
%             if length(allSpeedStim{mouse}{run}{stim})>=15/bin;
%             allSpeedStimAll=[allSpeedStimAll;allSpeedStim{mouse}{run}{stim}(1:15/bin)];
%             end
%             if length(allSpeedBefore{mouse}{run}{stim})>=15/bin;
%             allSpeedBeforeAll=[allSpeedBeforeAll;allSpeedBefore{mouse}{run}{stim}(end-15/bin+1:end)];
%             end
%             if length(allSpeedAfter{mouse}{run}{stim})>=15/bin;
%             allSpeedAfterAll=[allSpeedAfterAll;allSpeedAfter{mouse}{run}{stim}(1:15/bin)];
%             end
%         end
%         end
%      end
% end
% 
% %% 
% 
% figure
% x=[0:bin:15]+bin/2;
% x=x(x<=15);
% 
% errorbar(x,nanmean(allSpeedBeforeAll,1),nansem(allSpeedBeforeAll,1),'k');
% hold on
% errorbar(15+x,nanmean(allSpeedStimAll,1),nansem(allSpeedStimAll,1),'b');
% hold on
% errorbar(30+x,nanmean(allSpeedAfterAll,1),nansem(allSpeedAfterAll,1),'k');
% title(['beforeStimAfter',num2str(bin),'cm']);
% 
% saveas(gcf,'chr2RS_beforeStimAfter2_5cm.fig')
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
