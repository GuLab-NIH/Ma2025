function [predLickRBR,predLickAllRuns] = predLickingPercentileRollingWindow_RBR(allLick,allReward,allPosition,thresholdLick,beforeDistance,afterDistance)

%different from the previous ones: use before distance as window, calculate
%number of licks per window, then calculate the percentile among the
%remaining windows

%RBR predictive licking: percentage of licking among all licks

%allLick: each cell is one session, within the cell, there are n runs (10
%here): those are lick signals
%allReward: same to allLick and contains reward signals
%allPosition: same to allLick and contains position signal
%thresholdLick: threshold for detecting real licks
%beforeDistace: distance before reward to determine predictive lick.
%Normally 20cm.
%afterDistance: distance after reward to remove reward consumption lick.
%Normally 30cm.

predLickRBR=[];%each column is one session, each row is one run
predLickAllRuns=[];%each column is one session
for session=1:length(allLick);
    predLickIdxAll=[];
    otherLickAllPositions=[];
    for lap=1:length(allLick{session});
lick=allLick{session}{lap};
reward=allReward{session}{lap};
position=allPosition{session}{lap};
thresh=thresholdLick(session);

lickIdx=find(lick>=thresh);
if ~isempty(lickIdx)
rewardIdx=find(reward>0.5);
lickPosition=position(lickIdx);
rewardPosition=position(rewardIdx(1));

lickIdxBefore=lickIdx(lickIdx<rewardIdx);
lickIdxAfter=lickIdx(lickPosition>(afterDistance+rewardPosition));

lickIdxUse=[lickIdxBefore;lickIdxAfter];

if ~isempty(lickIdxUse)
lickIdxUsePosition=position(lickIdxUse);

% lickIdxUse=lickIdxUse(lickIdxUsePosition>30);
% lickIdxUsePosition=position(lickIdxUse);

predLickIdx=intersect(find(lickIdxUsePosition>=(rewardPosition-beforeDistance)),find(lickIdxUsePosition<rewardPosition));
otherLickIdx=setdiff([1:1:length(lickIdxUsePosition)],predLickIdx);
otherLickPosition=lickIdxUsePosition(otherLickIdx);

% trackZones1=[mod((rewardPosition-beforeDistance),beforeDistance):beforeDistance:(rewardPosition-beforeDistance)];%before reward
trackZones1=[];
for i=1:floor(rewardPosition-beforeDistance)-beforeDistance+1;
    trackZones1=[trackZones1;[i-1 i-1+beforeDistance]];
end

if round(max(position))-(rewardPosition+beforeDistance+afterDistance)>=0;
    trackZones2=[];%after reward zone
    for i=1:ceil(max(position))-(rewardPosition+beforeDistance+afterDistance)+1;
        trackZones2=[trackZones2;[rewardPosition+afterDistance+i-1 rewardPosition+afterDistance+i-1+beforeDistance]];
    end
else
    trackZones2=[];
end
trackZones=[trackZones1;trackZones2];

nData=[];
for i=1:size(trackZones,1);
    nData(i)=length(find(otherLickPosition>trackZones(i,1)&otherLickPosition<=trackZones(i,2)));
end


predLickRBR(session,lap)=length(find(nData<length(predLickIdx)))/(length(nData));%percentile among all
% predLickRBR(session,lap)=length(predLickIdx);

otherLickAllPositions=[otherLickAllPositions;otherLickPosition];
predLickIdxAll=[predLickIdxAll;predLickIdx];
else
    predLickRBR(session,lap)=nan;
end

else
    predLickRBR(session,lap)=nan;
end

    end

    if ~isempty([otherLickAllPositions;predLickIdxAll])

        nDataAllLap=[];
for i=1:size(trackZones,1);
    nDataAllLap(i)=length(find(otherLickAllPositions>trackZones(i,1)&otherLickAllPositions<=trackZones(i,2)));
end


        predLickAllRuns(session,1)=length(find(nDataAllLap<length(predLickIdxAll)))/(length(nDataAllLap));;
        % predLickAllRuns(session,1)=length(predLickIdxAll);
    else
        predLickAllRuns(session,1)=nan;
    end

 


end




end