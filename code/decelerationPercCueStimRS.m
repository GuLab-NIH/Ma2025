function [decePercentileBefore,decePercentileStim,decePercentileAfter] = decelerationPercCueStimRS(speedBefore,speedStim,speedAfter,yBefore,yStim,yAfter,stimLocAll,stimExp1,stimExp2,samplingSpacing)

%this code calculate RS based on RBR changed stim zones. before and after
%stim runs just used the corresponding stim loc in the stim run (e.g.,all run 1
%in the three period used the same stim loc)

%calculate the percentage of licks near stim locations
%speedBefore: every cell is a session, within: each run. all
%licks after reward 366 was removed
%stimLocAll: ROW: session, column: run. within: three rows, each row is the start and end of the stimulation
%stimExp: expansion of stimulation zone: for now use 10 (to both sides)
%samplingSpacing: how dense to sample the other zones: 0.5cm, for example

%output: lickPercentile: per run

%first get other zones:
%get other zone by doing rolling average
% stimLocExp=[];
% stimLocExp(:,1)=stimLoc(:,1)+stimExp1; %assymetrical expansion
% stimLocExp(:,2)=stimLoc(:,2)+stimExp2;
% 
% oneStimZoneLength=nanmean(stimLoc(:,2)-stimLoc(:,1))-stimExp1+stimExp2;%15cm after exppansion
% 
% otherZones=[];
% otherZones(1,:)=[0 stimLocExp(1,1)];
% otherZones(2,:)=[stimLocExp(1,2) stimLocExp(2,1)];
% otherZones(3,:)=[stimLocExp(2,2) stimLocExp(3,1)];
% otherZones(4,:)=[stimLocExp(3,2) 366];
% 
% allOtherZones=[];
% for n=1:size(otherZones,1);
%     NSample=floor((otherZones(n,2)-oneStimZoneLength)/samplingSpacing);
%     for i=1:NSample+1
%         allOtherZones=[allOtherZones;[otherZones(n,1)+(i-1)*samplingSpacing otherZones(n,1)+(i-1)*samplingSpacing+oneStimZoneLength]];
%     end
% end

decePercentileBefore=[];%each column is one session, each row is one run
decePercentileStim=[];%each column is one session, each row is one run
decePercentileAfter=[];%each column is one session, each row is one run
window=80;
for session=1:length(speedBefore)
    for run=1:length(speedBefore{session});


stimLoc=stimLocAll{session,run};

stimLocExp=[];
stimLocExp(:,1)=stimLoc(:,1)-stimExp1; %assymetrical expansion
stimLocExp(:,2)=stimLoc(:,2)+stimExp2;

oneStimZoneLength=nanmean(stimLoc(:,2)-stimLoc(:,1))+stimExp1+stimExp2;%15cm after exppansion

otherZones=[];
otherZones(1,:)=[0 stimLocExp(1,1)];
otherZones(2,:)=[stimLocExp(1,2) stimLocExp(2,1)];
otherZones(3,:)=[stimLocExp(2,2) stimLocExp(3,1)];
otherZones(4,:)=[stimLocExp(3,2) 366];

allOtherZones=[];
for n=1:size(otherZones,1);
    NSample=floor((otherZones(n,2)-oneStimZoneLength)/samplingSpacing);
    for i=1:NSample+1
        allOtherZones=[allOtherZones;[otherZones(n,1)+(i-1)*samplingSpacing otherZones(n,1)+(i-1)*samplingSpacing+oneStimZoneLength]];
    end
end


        beforeSpeed=speedBefore{session}{run};
        beforeSpeed=smoothdata(beforeSpeed,'movmean',window);

        stimSpeed=speedStim{session}{run};
        stimSpeed=smoothdata(stimSpeed,'movmean',window);

        afterSpeed=speedAfter{session}{run};
        afterSpeed=smoothdata(afterSpeed,'movmean',window);        

        beforeY=yBefore{session}{run};
        beforeY=smoothdata(beforeY,'movmean',window);

        stimY=yStim{session}{run};
        stimY=smoothdata(stimY,'movmean',window);

        afterY=yAfter{session}{run};
        afterY=smoothdata(afterY,'movmean',window);

        inZoneBefore=[];%within stim zone
        inZoneStim=[];
        inZoneAfter=[];
        for n=1:size(stimLocExp,1)
            i=find(beforeY>=stimLocExp(n,1)&beforeY<=stimLocExp(n,2));
            s=beforeSpeed(i);
            dece=nanmean(diff(s));
            inZoneBefore=[inZoneBefore;dece];
            
            i=find(stimY>=stimLocExp(n,1)&stimY<=stimLocExp(n,2));
            s=stimSpeed(i);
            dece=nanmean(diff(s));
            inZoneStim=[inZoneStim;dece];

            i=find(afterY>=stimLocExp(n,1)&afterY<=stimLocExp(n,2));
            s=afterSpeed(i);
            dece=nanmean(diff(s));
            inZoneAfter=[inZoneAfter;dece];
        end

        inZoneBefore=nanmean(inZoneBefore);
        inZoneStim=nanmean(inZoneStim);
        inZoneAfter=nanmean(inZoneAfter);

        inZoneBeforeAllOtherZones=[];%in other zones
        inZoneStimAllOtherZones=[];
        inZoneAfterAllOtherZones=[];

        for n=1:size(allOtherZones,1)
            i=find(beforeY>allOtherZones(n,1)&beforeY<allOtherZones(n,2));
            s=beforeSpeed(i);
            dece=nanmean(diff(s));
            inZoneBeforeAllOtherZones(n)=dece;

            i=find(stimY>allOtherZones(n,1)&stimY<allOtherZones(n,2));
            s=stimSpeed(i);
            dece=nanmean(diff(s));
            inZoneStimAllOtherZones(n)=dece;

            i=find(afterY>allOtherZones(n,1)&afterY<allOtherZones(n,2));
            s=afterSpeed(i);
            dece=nanmean(diff(s));
            inZoneAfterAllOtherZones(n)=dece;
        end

        inZoneBeforeShuffle=[];%in other zones
        inZoneStimShuffle=[];
        inZoneAfterShuffle=[];

        for n=1:1000;%shuffle
            i=randperm(size(allOtherZones,1));
            i=i(1:3);
            inZoneBeforeShuffle(n)=nanmean(inZoneBeforeAllOtherZones(i));   
            inZoneStimShuffle(n)=nanmean(inZoneStimAllOtherZones(i)); 
            inZoneAfterShuffle(n)=nanmean(inZoneAfterAllOtherZones(i)); 
        end

        if isempty(find([inZoneBefore inZoneBeforeShuffle]>0))
        decePercentileBefore(run,session)=nan;
        else
             decePercentileBefore(run,session)=length(find(inZoneBeforeShuffle>inZoneBefore))/1000;
        end

        if isempty(find([inZoneStim inZoneStimShuffle]>0))
        decePercentileStim(run,session)=nan;
        else
             decePercentileStim(run,session)=length(find(inZoneStimShuffle>inZoneStim))/1000;
        end

        if isempty(find([inZoneAfter inZoneAfterShuffle]>0))
        decePercentileAfter(run,session)=nan;
        else
             decePercentileAfter(run,session)=length(find(inZoneAfterShuffle>inZoneStim))/1000;
        end
    end
end
        


end

