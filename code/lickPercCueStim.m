function [lickPercentileBefore,lickPercentileStim,lickPercentileAfter] = lickPercCueStim(lickUsePositionBeforeRBR,lickUsePositionStimRBR,lickUsePositionAfterRBR,stimLoc,stimExp,samplingSpacing)
%calculate the percentage of licks near stim locations
%lickUsePositionBeforeRBR: every cell is a session, within: each run. all
%licks after reward 366 was removed
%stimLoc: three rows, each row is the start and end of the stimulation
%stimExp: expansion of stimulation zone: for now use 10 (to both sides)
%samplingSpacing: how dense to sample the other zones: 0.5cm, for example

%output: lickPercentile: per run

%first get other zones:
%get other zone by doing rolling average
stimLocExp=[];
stimLocExp(:,1)=stimLoc(:,1)-stimExp;
stimLocExp(:,2)=stimLoc(:,2)+stimExp;

oneStimZoneLength=nanmean(stimLoc(:,2)-stimLoc(:,1))+stimExp*2;%15cm after exppansion

otherZones=[];
otherZones(1,:)=[0 stimLocExp(1,1)];
otherZones(2,:)=[stimLocExp(1,2) stimLocExp(2,1)];
otherZones(3,:)=[stimLocExp(2,2) stimLocExp(3,1)];
otherZones(4,:)=[stimLocExp(3,2) 366];

allOtherZones=[];
for n=1:size(otherZones,1);
    NSample=(otherZones(n,2)-oneStimZoneLength)/samplingSpacing;
    for i=1:NSample+1
        allOtherZones=[allOtherZones;[otherZones(n,1)+(i-1)*samplingSpacing otherZones(n,1)+(i-1)*samplingSpacing+oneStimZoneLength]];
    end
end

lickPercentileBefore=[];%each column is one session, each row is one run
lickPercentileStim=[];%each column is one session, each row is one run
lickPercentileAfter=[];%each column is one session, each row is one run

for session=1:length(lickUsePositionBeforeRBR)
    for run=1:length(lickUsePositionBeforeRBR{session});
        beforeLick=lickUsePositionBeforeRBR{session}{run};
        stimLick=lickUsePositionStimRBR{session}{run};
        afterLick=lickUsePositionAfterRBR{session}{run};

        inZoneBefore=[];%within stim zone
        inZoneStim=[];
        inZoneAfter=[];
        for n=1:size(stimLocExp,1)
            licks=beforeLick(beforeLick>=stimLocExp(n,1)&beforeLick<=stimLocExp(n,2));
            inZoneBefore=[inZoneBefore;licks];
            
            licks=stimLick(stimLick>=stimLocExp(n,1)&stimLick<=stimLocExp(n,2));
            inZoneStim=[inZoneStim;licks];

            licks=afterLick(afterLick>=stimLocExp(n,1)&afterLick<=stimLocExp(n,2));
            inZoneAfter=[inZoneAfter;licks];
        end

        inZoneBefore=length(inZoneBefore);
        inZoneStim=length(inZoneStim);
        inZoneAfter=length(inZoneAfter);

        inZoneBeforeAllOtherZones=[];%in other zones
        inZoneStimAllOtherZones=[];
        inZoneAfterAllOtherZones=[];

        for n=1:size(allOtherZones,1)
            licks=beforeLick(beforeLick>allOtherZones(n,1)&beforeLick<allOtherZones(n,2));%did not use "=" so it is not overlapping with stimzone
            inZoneBeforeAllOtherZones(n)=length(licks);

            licks=stimLick(stimLick>allOtherZones(n,1)&stimLick<allOtherZones(n,2));%did not use "=" so it is not overlapping with stimzone
            inZoneStimAllOtherZones(n)=length(licks);  

            licks=afterLick(afterLick>allOtherZones(n,1)&afterLick<allOtherZones(n,2));%did not use "=" so it is not overlapping with stimzone
            inZoneAfterAllOtherZones(n)=length(licks);  
        end

        inZoneBeforeShuffle=[];%in other zones
        inZoneStimShuffle=[];
        inZoneAfterShuffle=[];

        for n=1:1000;%shuffle
            i=randperm(size(allOtherZones,1));
            i=i(1:3);
            inZoneBeforeShuffle(n)=sum(inZoneBeforeAllOtherZones(i));   
            inZoneStimShuffle(n)=sum(inZoneStimAllOtherZones(i)); 
            inZoneAfterShuffle(n)=sum(inZoneAfterAllOtherZones(i)); 
        end

        if isempty(find([inZoneBefore inZoneBeforeShuffle]>0))
        lickPercentileBefore(run,session)=nan;
        else
             lickPercentileBefore(run,session)=length(find(inZoneBeforeShuffle<inZoneBefore))/1000;
        end

        if isempty(find([inZoneStim inZoneStimShuffle]>0))
        lickPercentileStim(run,session)=nan;
        else
             lickPercentileStim(run,session)=length(find(inZoneStimShuffle<inZoneStim))/1000;
        end

        if isempty(find([inZoneAfter inZoneAfterShuffle]>0))
        lickPercentileAfter(run,session)=nan;
        else
             lickPercentileAfter(run,session)=length(find(inZoneAfterShuffle<inZoneStim))/1000;
        end
    end
end
        











end

