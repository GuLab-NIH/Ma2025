function [rewardSegAfter,lickSegAfter,positionSegAfter,timeSegAfter,rewardSegBefore,lickSegBefore,positionSegBefore,timeSegBefore,rewardSeg,lickSeg,positionSeg,timeSeg,stimZoneIdxAll,stimZoneLengthAll,idxStim] = within_ouside_calculation_all_5cm_yg3StimStartEndLickNoRemove(logPath)
%within_ouside_calculation_all_5cm_yg3StimStartEndLick included lick.
%This one did not remove the before reward stimulation

%% This script calculates the speed and licking within and ouside the stimulation
% %% (1) determine the lick threshold
% prompt = {'licking threshold'};
% dlgtitle = 'Input';
% dims = [1 35];
% definput = {'0.62'};
% answer = inputdlg(prompt,dlgtitle,dims,definput);
% lick_threshold = str2num(answer{1});
% save('lick_threshold.mat','lick_threshold');

%% (2) load data
logData = readLog(logPath);

%% (3) retrieve the stimulation trial numbers and stimulation state
state_all = logData(:,10);
TIME = (logData(:,1)-logData(1,1))*24*60;
% licks = logData(:,9);
% licks(licks<lick_threshold) = nan; % YG: changed threshold according to the rig lick signals 
% licks(licks>=lick_threshold) = 1;
% find the start point of each run
start = find(logData(:,3)==0);
segment = diff(start);
realstart = find(segment>300);
StartIdx = start(realstart);



for i = 1:length(StartIdx)-1 % didn't count the last run
    Position_Segment{i} = logData(StartIdx(i):StartIdx(i+1)-1,3);
    Reward_segment{i} = logData(StartIdx(i):StartIdx(i+1)-1,8);
    Lick_segment{i} = logData(StartIdx(i):StartIdx(i+1)-1,9);
    state_segment{i} = state_all(StartIdx(i):StartIdx(i+1)-1);
    time_segment{i} = TIME(StartIdx(i):StartIdx(i+1)-1);
    
    [M,I] = max([Position_Segment{i}]);
    Position_Segment{i} = [Position_Segment{i}(1:I)];
    Reward_segment{i} = [Reward_segment{i}(1:I)];
    Lick_segment{i} = [Lick_segment{i}(1:I)];
    state_segment{i} = [state_segment{i}(1:I)];
    time_segment{i} = [time_segment{i}(1:I)];
    
    % Calculate velocity correctly
    displacement = diff([Position_Segment{i}]);
    time_interval = diff([time_segment{i}]);
    Velocity_Segment{i} = (displacement./time_interval)/60;
    % now align at the reward location
    if find([Reward_segment{i}(1:I)]>0) 
        try
            RewardIdx(i) = find([Reward_segment{i}(1:I)]>0);
        catch
            RewardIdx(i) = find([Reward_segment{i}(1:I)]>0,1,'first');
        end
    else
        continue
    end
end

for i = 1:length(StartIdx)-1
    stimulation(i) = nanmean([state_segment{i}]);
end
stimulation(stimulation>0) =1;

%% BY YG
positionSeg=Position_Segment(find(stimulation==1));
timeSeg=time_segment(find(stimulation==1));
lickSeg=Lick_segment(find(stimulation==1));
rewardSeg=Reward_segment(find(stimulation==1));

idxStim=find(stimulation==1);

positionSegBefore=Position_Segment(idxStim(1)-10:idxStim(1)-1);
timeSegBefore=time_segment(idxStim(1)-10:idxStim(1)-1);
lickSegBefore=Lick_segment(idxStim(1)-10:idxStim(1)-1);
rewardSegBefore=Reward_segment(idxStim(1)-10:idxStim(1)-1);

positionSegAfter=Position_Segment(idxStim(end)+1:idxStim(end)+10);
timeSegAfter=time_segment(idxStim(end)+1:idxStim(end)+10);
lickSegAfter=Lick_segment(idxStim(end)+1:idxStim(end)+10);
rewardSegAfter=Reward_segment(idxStim(end)+1:idxStim(end)+10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate how many sections of no-stim and stim runs
no_stim_trial = find(stimulation<1);
stim_trial = find(stimulation>0);
last_run = 0;
if last_run>0
    no_stim_trial(no_stim_trial>last_run) = [];
    stim_trial(stim_trial>last_run) = []; % stim_trials are the trials used for analysis
end


%% (4) Calculate the speed within and outside the stim zone (Reward zone excluded)
vel_within_stim=[];
vel_outside_stim=[];

stimZoneIdxAll={};
stimZoneLengthAll={};
for i = 1:length(stim_trial)
    disp(i)
%     trial_num = stim_trial(i);
%     place = [Position_Segment{trial_num}];
%     stim_state = [state_segment{trial_num}];
%     stim_state(place>=366-60)=-1;%rempving the area too close to the reward
%     stim_state(1) = [];
%     
%     vel = [Velocity_Segment{trial_num}]';
%     if length(stim_state) - length(vel) ~=0
%         difference = length(stim_state) - length(vel);
%         stim_state(1:difference) =[];
%     end
%     vel_within_stim{i} = vel(stim_state>0);
%     vel_outside_stim{i} = vel(stim_state==0);    


 trial_num = stim_trial(i);
    place = [Position_Segment{trial_num}];
    stim_state = [state_segment{trial_num}];
    % stim_state(place>=366-60)=-1;
    stim_state(place<=4)=-1;%random never occured to 0-4cm
    stim_state(1) = [];
    
    %lick = [Lick_segment{trial_num}];
    position = [Position_Segment{trial_num}];
%     lick(isnan(lick))=0;
%     stim_location = round(position(stim_state>0));
    stim_location=position(stim_state>0);%real positino without rounding

    stim_locationRealIdx=find(stim_state>0);%yg: stim locations

    %% added by YG
    stim_LocationRealIdxZone=contiguous([1;diff(stim_locationRealIdx)],1);
     stim_LocationRealIdxZone= stim_LocationRealIdxZone{1,2};
     allIdx=[1:1:length(stim_state)];
     stim_LocationRealIdxZone2=stim_locationRealIdx(stim_LocationRealIdxZone);    
     stimZoneIdxAll{i}=stim_LocationRealIdxZone2;
%above: added by YG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    start_stim_diff = diff(stim_location);
    stim_start_idx = find(start_stim_diff>20)+1;%this 20 cm distance is ok for 5cm random because they were set to be 45 cm apart already
    stim_start_idx = [1; stim_start_idx];

    stim_start = stim_location(stim_start_idx);%this is location
    % stim_start(stim_start>=366-60) = [];
    stim_end=[];%this is location
    
    for ii=1:length(stim_start);
        if ii~=length(stim_start);
        stim_end(ii)=stim_location(stim_start_idx(ii+1)-1);
        stim_end_idx(ii)=stim_start_idx(ii+1)-1;
        else
            % stim_end(ii)=min([max(stim_location) 366-60]);
             stim_end(ii)=max(stim_location);
        end
    end
    stimZoneLengthAll{i}=stim_end'-stim_start;%by YG

    stimLength=sum(stim_end-stim_start');
    stimTime=1/60*length(find(stim_state>0));%in second, because 60hz per samble
    vel_within_stim(i)=stimLength/stimTime;

    % noStimLength=366-60-4-stimLength;
     noStimLength=366-4-stimLength;
    noStimTime=1/60*length(find(stim_state==0));
     vel_outside_stim(i)=noStimLength/noStimTime;


end

mean_vel_within = nanmean(vel_within_stim);
mean_vel_outside = nanmean(vel_outside_stim);
% disp('outside and within speed')
% [mean_vel.outside_stim mean_vel.within_stim]
% save('within_outside_speed.mat','mean_vel')


% %% (5) calculate licking within and outside the stimulation zone (reward zone excluded)
% for i = 1:length(stim_trial)
%     trial_num = stim_trial(i);
%     place = [Position_Segment{trial_num}];
%     stim_state = [state_segment{trial_num}];
%     stim_state(place>=366-60)=-1;
%         stim_state(place<=4)=-1;%random never occured to 0-4cm
%     stim_state(1) = [];
% 
%     lick = [Lick_segment{trial_num}];
%     position = [Position_Segment{trial_num}];
%     lick(isnan(lick))=0;
% %     stim_location = round(position(stim_state>0));
%     stim_location=position(stim_state>0);%real positino without rounding
%     start_stim_diff = diff(stim_location);
%     stim_start_idx = find(start_stim_diff>20)+1;%this 20 cm distance is ok for 5cm random because they were set to be 45 cm apart already
%     stim_start_idx = [1; stim_start_idx];
%     stim_start = stim_location(stim_start_idx);
%     stim_start(stim_start>=366-60) = [];
%     stim_end=[];
% 
%     for ii=1:length(stim_start);
%         if ii~=length(stim_start);
%         stim_end(ii)=stim_location(stim_start_idx(ii+1)-1);
%         else
%             stim_end(ii)=min([max(stim_location) 366-60]);
%         end
%     end
% 
% 
%     if ~isempty(find(lick>0))
%         n=contiguous(lick,[1]); %find contiguous 1
%         nn=n{1,2};
%         lickIdx=nn(:,1);%%%%%%%%%%%%%%%%%%%licking indices
%         yLick=position(lickIdx);%y position of licking
%         yLick(yLick>=366-60) = [];
%          yLick(yLick<=4) = [];
%         Lick_within_stim(i) = 0;
%         for j = 1:length(stim_start)%10
%             for k = 1:length(yLick)
%     %             if yLick(k) >vr.RanPos{stim_trial(i)}(j) && yLick(k) < vr.RanPos{stim_trial(i)}(j)+10
% %                 if yLick(k) >stim_start(j) && yLick(k) < stim_start(j)+5 %should be changed to 5 for three cue random, shoudl be 15 for BR random
%                 if yLick(k) >stim_start(j) && yLick(k) < stim_end(j) %should be changed to 5 for three cue random, shoudl be 15 for BR random
% 
% Lick_within_stim(i) = Lick_within_stim(i) +1;
%                 end
%             end
%         end
%         Lick_outside_stim(i) = length(yLick) - Lick_within_stim(i);
% 
% 
%         num_stim = find(diff(stim_location)>2);
%         num_stim1 = length(num_stim)+1;% should be length(num_stim)+1
%         totalDistansceStim=sum(stim_end'-stim_start);
%         meanLick_within_stim(i) = Lick_within_stim(i)/totalDistansceStim;%10 should be changed to 5 for the random of 5cm at cues, should be changed to 15 for the 15cm random
%         meanLick_outside_stim(i) = Lick_outside_stim(i)/(366-60-4-totalDistansceStim);
%     else
%         meanLick_within_stim(i) = 0;
%         meanLick_outside_stim(i) = 0;
%     end
% end
% mean_lick_within = nanmean(meanLick_within_stim);
% mean_lick_outside = nanmean(meanLick_outside_stim);
% 
% % disp('outside and within licks')
% % [mean_lick.outside_stim mean_lick.within_stim ]
% % save('within_outside_licks.mat','mean_lick')

end
