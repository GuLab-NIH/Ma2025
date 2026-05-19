%% Open ephys data analysis
% First run "load_open_ephys_binary.m" to extract data for each recording,
% then modify and run "Extract_and_save_20211021.m" to organize the data format 
% then load the data and stimulus timestmaps

% modified from OpenEphys_analysis1.m by June Kim Nov. 14, 2025

function segData=segEphysData(stim, voltage, time, segTimeRange, chList, sampleRate, adjShift)
    segTimeRange = segTimeRange * sampleRate;

    stimBlockPeriod = 1; % time period of 1 stimulation sequence (sec)
    stimBlockPeriod = stimBlockPeriod * sampleRate;

    totalNoCh = length(chList);

    %% Segment the pulse trains
    PulseTrain_interval = [0 diff(stim)'];
    PulseSeg = [1 find(PulseTrain_interval>5*sampleRate)]; % the sampling rate is 30k/sec, 
    % so what it does is to find the index when the train interval is larger than 5 sec, we set the train interval to be 10 sec long
    NoPulseSeg=length(PulseSeg);
    % Get the pulse duration and frequency
    number = 1;
    first_stimON = PulseSeg; % the 1st LED ON of each pulse train
    stimuli_info = zeros(NoPulseSeg,2); % the pulse duration and freque
    
    for stim_num = 1:NoPulseSeg % should be 36 stimuli
        if stim_num < NoPulseSeg
            PulseTrain_time = [PulseSeg(stim_num):PulseSeg(stim_num+1)-1];% find the timestamps within one pulse train set
            Pulse_interval = diff(stim(PulseTrain_time));
            PulseDur = round(mean(Pulse_interval(1:2:end)))/sampleRate; % the pulse duration in sec
            PulseFreq = length(PulseTrain_time)/2; % the pulse frequency in Hz
            
            stimON{stim_num} = stim(PulseTrain_time(1:2:end)); % the timestamps for the stim ON within one pulse train
            stimOFF{stim_num} = stim(PulseTrain_time(2:2:end)); % the timestamps for the stim OFF within one pulse train
            stimuli_info(number,:) = [PulseDur PulseFreq]; 
            number = number+1;
        else
            PulseTrain_time = [PulseSeg(stim_num):length(PulseTrain_interval)];% find the timestamps within the last pulse train set
            Pulse_interval = diff(stim(PulseTrain_time));
            PulseDur = round(mean(Pulse_interval(1:2:end)))/sampleRate; % the pulse duration in sec
            PulseFreq = length(PulseTrain_time)/2; % the pulse frequency in Hz
            
            stimON{stim_num} = stim(PulseTrain_time(1:2:end))'; % the timestamps for the stim ON within one pulse train
            stimOFF{stim_num} = stim(PulseTrain_time(2:2:end))'; % the timestamps for the stim OFF within one pulse train
            stimuli_info(number,:) = [PulseDur PulseFreq]; 
        end
           
    end
    
    segData.stimON=stimON;
    segData.stimOFF=stimOFF;
    segData.stimuli_info=stimuli_info;

    %% Segment the filtered signals from given channels based on stimuli

    segData.SpikeSeg = cell(totalNoCh,1);
    segData.LFPSeg = cell(totalNoCh,1);
    % segData.rawSeg = cell(totalNoCh,1);
    for chIdx = 1:totalNoCh
    
        % Filter the signals based on frequency
        raw = voltage(chList(chIdx),:);
        % rawFiltered = notchFilter60Hz(raw, sampleRate);

        % Spike_filtered = bandpass(voltage(chList(chIdx),:),[250 8000],sampleRate); % 250-8000 Hz for spikes
        % Define the bandpass filter first
        bpFilt = designfilt('bandpassiir', ...
            'FilterOrder', 4, ...                % 4th order IIR (fast and stable)
            'HalfPowerFrequency1', 250, ...      % lower cut
            'HalfPowerFrequency2', 8000, ...     % upper cut
            'SampleRate', sampleRate);
        
        % Apply zero-phase filtering (no phase shift)
        Spike_filtered = filtfilt(bpFilt, raw);
        % Spike_filtered = detrend(raw,'constant');

        % For LFP signals
        Hd = LFP_bandpass_filter1; % call the filter function
        LFP_filtered = filtfilt(Hd, raw); % 0.7-170 Hz for LFP
        
        % Segment the recording from 0.5 sec before the stim ON, 5 sec after the
        % stim OFF
        start = zeros(1,NoPulseSeg); % the timestamps 0.5 sec before each pulse train
        stop = zeros(1,NoPulseSeg); % the timestamps 5 sec after each pulse train
        
        for stim_num = 1:NoPulseSeg % should be 36 stimuli
            start(stim_num) = stim(first_stimON(stim_num)) - segTimeRange(1) - time(1) + adjShift; % sampleRate per sec for sampling rate
            stop(stim_num) = start(stim_num) + segTimeRange(1) + stimBlockPeriod + segTimeRange(2) - 1 + adjShift;% the start time plus 1 sec of stimulus and 5 sec after stimulus off
            
            % segment the filtered data
            Spike_filtered_seg(stim_num,:) = Spike_filtered(start(stim_num):stop(stim_num))*0.195; %*0.195 convert to uV
            LFP_filtered_seg(stim_num,:) = LFP_filtered(start(stim_num):stop(stim_num))*0.195 ; %* convert to uV   
            % raw_seg(stim_num,:) = rawFIltered(start(stim_num):stop(stim_num))*0.195; %*0.195 convert to uV
        end
        segData.SpikeSeg{chList(chIdx)} = Spike_filtered_seg;
        segData.LFPSeg{chList(chIdx)} = LFP_filtered_seg;
        % segData.rawSeg{chList(chIdx)} = raw_seg;
    end
end