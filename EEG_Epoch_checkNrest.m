function [Epochs_check, Epochs_rest, EEG_events] = EEG_Epoch_checkNrest(EEG_input, videoTrigger, checkName, restName)
events = EEG_input.event;

%% 
%tag events for epochs, tag events every 2s after video trigger
%0-20s: Rest
%20-40s: checkerboard

%find the latencies of the video starts
Latencies_triggers = [];
triggerIdx = 1;
for eventIdx = 1:size(events,2)
   event = events(eventIdx).type;
   if eventIdx == 101
       thing = 1;
   end
   if strcmp(events(eventIdx).type, videoTrigger)
       Latencies_triggers(triggerIdx) = events(eventIdx).latency;
       triggerIdx = triggerIdx + 1;
   end
end
%find events for the first 20s (from 0-18s) after video trigger every 2 seconds for rest epochs
restDur = 18 * EEG_input.srate;
epochDur = 2 * EEG_input.srate;

EEG_events = EEG_input;
for i = 1:size(Latencies_triggers, 2)
    EEG_events = EEG_Epoch_Periodic(EEG_events, ...
        Latencies_triggers(i), ...
        Latencies_triggers(i) + restDur, ...
        epochDur, restName);
end

%find events for the second 20s (20s-38s) after video trigger every 2 seconds for checkerboard epochs
checkOffset_start = restDur + (2 * EEG_input.srate);
checkOffset_end = checkOffset_start + (18 * EEG_input.srate);
epochDur = 2 * EEG_input.srate;
for i = 1:size(Latencies_triggers, 2)
    EEG_events = EEG_Epoch_Periodic(EEG_events, ...
        Latencies_triggers(i) + checkOffset_start, ...
        Latencies_triggers(i) + checkOffset_end, ...
        epochDur, checkName);
end

%% Epoch rest events from 0-2s before and after each rest event

Epochs_rest = pop_epoch( EEG_events, {restName}, [0  2]);
Epochs_check = pop_epoch( EEG_events, {checkName}, [0  2]);



    

