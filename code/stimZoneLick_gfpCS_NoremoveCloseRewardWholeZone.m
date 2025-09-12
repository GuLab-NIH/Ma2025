
%only use the random stim that are started from 50 and ended by 275. The
%actual random stim can go from ~5cm to 366cm
%just plot the "allPositions" below to see
% figure
% for n=1:size(allPositions,1);
%     hold on
%     line([allPositions(n,1) allPositions(n,2)],[n n]);
% end
% hold on
% line([366 366],[0 300])

%% gfp mice only
startDis=0;
endDis=400;

% gfp mice only
load('threshold.mat')
bin=1; %cm
beforeAfterDis=15;%15cm before and after speed
allLickStim={};
allLickBefore={};
allLickAfter={};
allPositions=[];
% for mouse=1:length(stimZoneLengthAllData);
mouseUse=[45:56];
for mouse=1:length(mouseUse);
    disp(mouse)
    allLickStim{mouse}={};
    allLickBefore{mouse}={};
    allLickAfter{mouse}={};
    for run=1:length(stimZoneLengthAllData{mouseUse(mouse)});
        disp(run)
          allLickStim{mouse}{run}={};
    allLickBefore{mouse}{run}={};
    allLickAfter{mouse}{run}={};
        if ~isempty(stimZoneLengthAllData{mouseUse(mouse)}{run});%stimulation should be at least 14cm
            stimIdx=stimZoneIdxAllDataFix{mouseUse(mouse)}{run};
            position=positionSegAllData{mouseUse(mouse)}{run};
            lick=lickSegAllData{mouseUse(mouse)}{run};
            lick(lick>=threshold(mouseUse(mouse)))=1;
            lick(lick<threshold(mouseUse(mouse)))=0;
            time=timeSegAllData{mouseUse(mouse)}{run};
            for nStim=1:size(stimIdx,1)
                allLickStim{mouse}{run}{nStim}=[];
    allLickBefore{mouse}{run}{nStim}=[];
    allLickAfter{mouse}{run}{nStim}=[];
    thisStimIdx=[stimIdx(nStim,1):1:stimIdx(nStim,2)];
    thisPosition=position(thisStimIdx)-position(thisStimIdx(1));

allPositions=[allPositions;[position(thisStimIdx(1)) position(thisStimIdx(end))]];
if position(thisStimIdx(1))>=startDis & position(thisStimIdx(end))<=endDis;

    thisTime=time(thisStimIdx)-time(thisStimIdx(1));
    thisLick=lick(thisStimIdx);

        allLickStim{mouse}{run}{nStim}=sum(thisLick)/(position(thisStimIdx(end))-position(thisStimIdx(1)));%normalize by distance


    % %during stim speed
    % positionBins=[0:bin:ceil(thisPosition(end))];
    % 
    % % positionBins=[0:bin:15];
    % [~,nbin]=histc(thisPosition,positionBins);
    % 
    % 
    % for binNumber=1:max(nbin);
    %     positionThisBin=thisPosition(nbin==binNumber);
    %     timeThisBin=thisTime(nbin==binNumber);
    %     lickThisBin=thisLick(nbin==binNumber);
    %     if max(timeThisBin)>0
    %     % allLickStim{mouse}{run}{nStim}(binNumber)=sum(lickThisBin)/((timeThisBin(end)-timeThisBin(1))*60); %the original time already in min 
    %             allLickStim{mouse}{run}{nStim}(binNumber)=sum(lickThisBin); %the original time already in min 
    %     end
    % end

    %calculating before sim speed
    stimStartPosition=position(thisStimIdx(1));%no subtraction by first one
    beforePosition=stimStartPosition-beforeAfterDis;
    iInclude=intersect(find(position>=beforePosition),find(position<stimStartPosition));

    thisPositionBefore=position(iInclude);
    thisPositionBefore=thisPositionBefore-thisPositionBefore(1);
    thisTimeBefore=time(iInclude);
    thisTimeBefore=thisTimeBefore-thisTimeBefore(1);
    thisLickBefore=lick(iInclude);

        allLickBefore{mouse}{run}{nStim}=sum(thisLickBefore)/(thisPositionBefore(end)-thisPositionBefore(1));


    % % positionBinsBefore=[0:bin:ceil(thisPositionBefore(end))];
    % positionBinsBefore=[0:bin:beforeAfterDis];
    % [~,nbinBefore]=histc(thisPositionBefore,positionBinsBefore);
    % 
    %  for binNumber=1:max(nbinBefore);
    %     positionThisBinBefore=thisPositionBefore(nbinBefore==binNumber);
    %     timeThisBinBefore=thisTimeBefore(nbinBefore==binNumber);
    %     lickThisBinBefore=thisLickBefore(nbinBefore==binNumber);
    %     if max(timeThisBinBefore)>0
    %     % allLickBefore{mouse}{run}{nStim}(binNumber)=sum(lickThisBinBefore)/((timeThisBinBefore(end)-timeThisBinBefore(1))*60); %the original time already in min 
    %     allLickBefore{mouse}{run}{nStim}(binNumber)=sum(lickThisBinBefore);
    %     end
    %  end

     %calculating after stim speed

 stimEndPosition=position(thisStimIdx(end));%no subtraction by first one
    afterPosition=stimEndPosition+beforeAfterDis;
    iInclude=intersect(find(position<=afterPosition),find(position>stimEndPosition));

    thisPositionAfter=position(iInclude);
    thisPositionAfter=thisPositionAfter-thisPositionAfter(1);
    thisTimeAfter=time(iInclude);
    thisTimeAfter=thisTimeAfter-thisTimeAfter(1);
    thisLickAfter=lick(iInclude);

            allLickAfter{mouse}{run}{nStim}=sum(thisLickAfter)/(thisPositionAfter(end)-thisPositionAfter(1));

    % 
    % % positionBinsAfter=[0:bin:ceil(thisPositionAfter(end))];
    % positionBinsAfter=[0:bin:beforeAfterDis];
    % [~,nbinAfter]=histc(thisPositionAfter,positionBinsAfter);
    % 
    %  for binNumber=1:max(nbinAfter);
    %     positionThisBinAfter=thisPositionAfter(nbinAfter==binNumber);
    %     timeThisBinAfter=thisTimeAfter(nbinAfter==binNumber);
    %     lickThisBinAfter=thisLickAfter(nbinAfter==binNumber);
    %     if max(timeThisBinAfter)>0
    %     % allLickAfter{mouse}{run}{nStim}(binNumber)=sum(lickThisBinAfter)/((timeThisBinAfter(end)-timeThisBinAfter(1))*60); %the original time already in min         
    %     allLickAfter{mouse}{run}{nStim}(binNumber)=sum(lickThisBinAfter);
    %     end
    %  end
     else
       allLickAfter{mouse}{run}{nStim}=[];
          allLickBefore{mouse}{run}{nStim}=[];
          allLickStim{mouse}{run}{nStim}=[];
end
            end
        else
          allLickAfter{mouse}{run}=[];
          allLickBefore{mouse}{run}=[];
          allLickStim{mouse}{run}=[];
        end
    end
end


% save('allSpeedAfter.mat','allSpeedAfter')
% save('allSpeedBefore.mat','allSpeedBefore')
% save('allSpeedStim.mat','allSpeedStim')

%% merge all mice runs together

allLickStimAll=[];
allLickAfterAll=[];
allLickBeforeAll=[];

allLickStimAllPerSession=[];
allLickAfterAllPerSession=[];
allLickBeforeAllPerSession=[];

yesOrNo=[];
load('stimZoneLengthAllData.mat')

stimLength=stimZoneLengthAllData(mouseUse);

for mouse=1:length(allLickStim)

stimSession=[];
afterSession=[];
beforeSession=[];

    for run=1:length(allLickStim{mouse});
         if ~isempty(allLickStim{mouse}{run});
        for stim=1:length(allLickStim{mouse}{run});
            if round(stimLength{mouse}{run}(stim))>10; %15 cm stimulations
                if ~isempty(allLickStim{mouse}{run}{stim})
            allLickStimAll=[allLickStimAll;allLickStim{mouse}{run}{stim}];
            yesOrNo(end+1,1)=1;
            stimSession=[stimSession;allLickStim{mouse}{run}{stim}];
            if ~isempty(allLickBefore{mouse}{run}{stim})
            allLickBeforeAll=[allLickBeforeAll;allLickBefore{mouse}{run}{stim}];
             yesOrNo(end,2)=1;
             beforeSession=[beforeSession;allLickBefore{mouse}{run}{stim}];
            else
                allLickBeforeAll=[allLickBeforeAll;nan];
                beforeSession=[beforeSession;nan];
                yesOrNo(end,2)=0;
            end
            if ~isempty(allLickAfter{mouse}{run}{stim})
            allLickAfterAll=[allLickAfterAll;allLickAfter{mouse}{run}{stim}];
            yesOrNo(end,3)=1;
            afterSession=[afterSession;allLickAfter{mouse}{run}{stim}];
            else
                allLickAfterAll=[allLickAfterAll;nan];
                afterSession=[afterSession;nan];
                yesOrNo(end,3)=0;
            end
                end
            end
        end
        end
    end
    allLickStimAllPerSession(mouse,1)=nanmean(stimSession);
allLickAfterAllPerSession(mouse,1)=nanmean(afterSession);
allLickBeforeAllPerSession(mouse,1)=nanmean(beforeSession);
end
allLickStimAll(isinf(allLickStimAll))=nan;
allLickAfterAll(isinf(allLickAfterAll))=nan;
allLickBeforeAll(isinf(allLickBeforeAll))=nan;

%find the stimulations with all components
useIdx=find(sum(yesOrNo,2)==3); 

allLickAfterAllUse=allLickAfterAll(useIdx,:);
allLickStimAllUse=allLickStimAll(useIdx,:);
allLickBeforeAllUse=allLickBeforeAll(useIdx,:);
save('gfpCSLickUseWholeZoneNoRemoveReward.mat','allLickAfterAllUse','allLickBeforeAllUse','allLickStimAllUse','allLickAfterAllPerSession','allLickBeforeAllPerSession','allLickStimAllPerSession')


%% 
% 
% figure
% x=[0:bin:15]+bin/2;
% x=x(x<=15);
% 
% errorbar(x,nanmean(allLickBeforeAll,1),nansem(allLickBeforeAll,1),'k');
% hold on
% errorbar(15+x,nanmean(allLickStimAll,1),nansem(allLickStimAll,1),'b');
% hold on
% errorbar(30+x,nanmean(allLickAfterAll,1),nansem(allLickAfterAll,1),'k');
% title(['beforeStimAfter',num2str(bin),'cm']);
% 
% saveas(gcf,'gfpRSLick_beforeStimAfter1cmRemoveReward.fig')

%  %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %%
% bin=2.5; %cm
% beforeAfterDis=15;%15cm before and after speed
% allSpeedStim={};
% allSpeedBefore={};
% allSpeedAfter={};
% mouseUse=[28:48];
% for mouse=1:length(mouseUse);
%     disp(mouse)
%     allSpeedStim{mouse}={};
%     allSpeedBefore{mouse}={};
%     allSpeedAfter{mouse}={};
%     for run=1:length(stimZoneLengthAllData{mouse});
%         disp(run)
%           allSpeedStim{mouse}{run}={};
%     allSpeedBefore{mouse}{run}={};
%     allSpeedAfter{mouse}{run}={};
%         if ~isempty(stimZoneLengthAllData{mouseUse(mouse)}{run});%stimulation should be at least 14cm
%             stimIdx=stimZoneIdxAllDataFix{mouseUse(mouse)}{run};
%             position=positionSegAllData{mouseUse(mouse)}{run};
%             time=timeSegAllData{mouseUse(mouse)}{run};
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
% saveas(gcf,'gfpCS_beforeStimAfter2_5cm.fig')
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
