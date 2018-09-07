inputDir = '~/Documents/NKI/EEG_Data/Processed_Data/';


%% get epochs for inside scanner
inputName = 'EEG_fMRI_20180830_01_Checkerboard_Flash_Inside_04_bcg_bp.set';
fullfile_input = [inputDir, inputName];

EEG_Processed = pop_loadset(fullfile_input);

[Epochs_Check, Epochs_Rest, EEG_events] = ...
    EEG_Epoch_checkNrest(EEG_Processed, 'S  1', 'Check', 'Rest');


numEpochs = size(Epochs_Check.data,3);

%Get channels with ssvep response for average
%             O1, O2, P7, P8, P1, P2, PO3, PO4, PO7, PO8 
SSVEP_Chan = [ 9, 10, 15, 16, 37, 38,  45,  46,  59,  60]; 

%% make figures for checkboard inside scanner
figure
subplot(2,2,1)
title('Checkerboard Inside Scanner')
xlabel('time(s)')
ylabel('uV')
axis([0, 2, -8, 8])
hold

%find average of epochs for each channel
for chanIdx = 1:size(SSVEP_Chan, 2)
    sumEpochs = Epochs_Check.data(SSVEP_Chan(chanIdx),:,1);
    for epochIdx = 2:numEpochs
        sumEpochs = sumEpochs +  ...
            Epochs_Check.data(SSVEP_Chan(chanIdx), :, epochIdx);
    end
    
    EpochAvg(chanIdx, :) = sumEpochs ./ numEpochs;
    
    %plot each Channel average
    t = 0:1/5000:2-(1/5000);
    plot(t, EpochAvg(chanIdx, :))

    
end


%find average of each channel
sumChans = EpochAvg(1, :);
for chanIdx = 2:size(EpochAvg, 1)
    sumChans = sumChans + EpochAvg(chanIdx, :);
end
subplot(2,2,2)
title('Average of O1, O2, P7, P8, P1, P2, PO3, PO4, PO7, PO8')
xlabel('time(s)')
ylabel('uV')
axis([0, 2, -8, 8])
hold
chanAvg = sumChans ./ size(EpochAvg, 1);
plot(t, chanAvg)

%plot freqency domain of average chanAvg and of multiple epochs
subplot(2,2,3)
title('Frequency')
xlabel('Frequency(Hz)')
ylabel('uV')
axis([0, 60, 0, 1])
hold

Fs = 5000;
dt = 1/Fs;
N = size(chanAvg, 2);
dF = Fs/N;
f = -Fs/2:dF:Fs/2 - dF;

for chanIdx = 1:size(EpochAvg, 1)
    chanfft = fft(EpochAvg(chanIdx,:));
    chanfft = fftshift(chanfft);
    %plot from 0 to 60 Hz, 
    %NOTE: since fft is centered, 0 Hz = N/2 samples & bins are 0.5Hz
    semilogy(f(N/2 + 1:N/2+(60*2)),abs(chanfft(N/2 + 1:N/2 + (60*2)))/N)
end

subplot(2,2,4)
title('Frequency')
xlabel('Frequency(Hz)')
ylabel('uV')
axis([0, 20, 0, 1])
hold
for chanIdx = 1:size(EpochAvg, 1)
    chanfft = fft(EpochAvg(chanIdx,:));
    chanfft = fftshift(chanfft);
    %plot from 0 to 60 Hz, 
    %NOTE: since fft is centered, 0 Hz = N/2 samples & bins are 0.5Hz
    plot(f(N/2 + 1:N/2+(60*2)),abs(chanfft(N/2 + 1:N/2 + (60*2)))/N)
end
%% Get figures for rest inside scanner

figure
subplot(2,2,1)
title('Rest Inside Scanner')
xlabel('time(s)')
ylabel('uV')
axis([0, 2, -8, 8])
hold

%find average of epochs for each channel
for chanIdx = 1:size(SSVEP_Chan, 2)
    sumEpochs = Epochs_Rest.data(SSVEP_Chan(chanIdx),:,1);
    for epochIdx = 2:numEpochs
        sumEpochs = sumEpochs +  ...
            Epochs_Rest.data(SSVEP_Chan(chanIdx), :, epochIdx);
    end
    
    EpochAvg(chanIdx, :) = sumEpochs ./ numEpochs;
    
    %plot each Channel average
    t = 0:1/5000:2-(1/5000);
    plot(t, EpochAvg(chanIdx, :))

    
end


%find average of each channel
sumChans = EpochAvg(1, :);
for chanIdx = 2:size(EpochAvg, 1)
    sumChans = sumChans + EpochAvg(chanIdx, :);
end
subplot(2,2,2)
title('Average of O1, O2, P7, P8, P1, P2, PO3, PO4, PO7, PO8')
xlabel('time(s)')
ylabel('uV')
axis([0, 2, -8, 8])
hold
chanAvg = sumChans ./ size(EpochAvg, 1);
plot(t, chanAvg)

%plot freqency domain of average chanAvg and of multiple epochs
subplot(2,2,3)
title('Frequency')
xlabel('Frequency(Hz)')
ylabel('uV')
axis([0, 60, 0, 1])
hold

Fs = 5000;
dt = 1/Fs;
N = size(chanAvg, 2);
dF = Fs/N;
f = -Fs/2:dF:Fs/2 - dF;

for chanIdx = 1:size(EpochAvg, 1)
    chanfft = fft(EpochAvg(chanIdx,:));
    chanfft = fftshift(chanfft);
    %plot from 0 to 60 Hz, 
    %NOTE: since fft is centered, 0 Hz = N/2 samples & bins are 0.5Hz
    plot(f(N/2 + 1:N/2+(60*2)),abs(chanfft(N/2 + 1:N/2 + (60*2)))/N)
end

subplot(2,2,4)
title('Frequency')
xlabel('Frequency(Hz)')
ylabel('uV')
axis([0, 20, 0, 1])
hold
for chanIdx = 1:size(EpochAvg, 1)
    chanfft = fft(EpochAvg(chanIdx,:));
    chanfft = fftshift(chanfft);
    %plot from 0 to 60 Hz, 
    %NOTE: since fft is centered, 0 Hz = N/2 samples & bins are 0.5Hz
    plot(f(N/2 + 1:N/2+(60*2)),abs(chanfft(N/2 + 1:N/2 + (60*2)))/N)
end

%% get epochs for outside scanner
inputName = 'EEG_fMRI_20180830_Checkerboard_Flash_Outside_02_bcg_bp.set';
fullfile_input = [inputDir, inputName];

EEG_Processed = pop_loadset(fullfile_input);

[Epochs_Check, Epochs_Rest, EEG_events] = ...
    EEG_Epoch_checkNrest(EEG_Processed, 'S  1', 'Check', 'Rest');


numEpochs = size(Epochs_Check.data,3);

%Get channels with ssvep response for average
%             O1, O2, P7, P8, P1, P2, PO3, PO4, PO7, PO8 
SSVEP_Chan = [ 9, 10, 15, 16, 37, 38,  45,  46,  59,  60]; 

%% make figures for checkboard inside scanner
figure
subplot(2,2,1)
title('Checkerboard Outside Scanner')
xlabel('time(s)')
ylabel('uV')
axis([0, 2, -8, 8])
hold

%find average of epochs for each channel
for chanIdx = 1:size(SSVEP_Chan, 2)
    sumEpochs = Epochs_Check.data(SSVEP_Chan(chanIdx),:,1);
    for epochIdx = 2:numEpochs
        sumEpochs = sumEpochs +  ...
            Epochs_Check.data(SSVEP_Chan(chanIdx), :, epochIdx);
    end
    
    EpochAvg(chanIdx, :) = sumEpochs ./ numEpochs;
    
    %plot each Channel average
    t = 0:1/5000:2-(1/5000);
    plot(t, EpochAvg(chanIdx, :))

    
end


%find average of each channel
sumChans = EpochAvg(1, :);
for chanIdx = 2:size(EpochAvg, 1)
    sumChans = sumChans + EpochAvg(chanIdx, :);
end
subplot(2,2,2)
title('Average of O1, O2, P7, P8, P1, P2, PO3, PO4, PO7, PO8')
xlabel('time(s)')
ylabel('uV')
axis([0, 2, -8, 8])
hold
chanAvg = sumChans ./ size(EpochAvg, 1);
plot(t, chanAvg)

%plot freqency domain of average chanAvg and of multiple epochs
subplot(2,2,3)
title('Frequency')
xlabel('Frequency(Hz)')
ylabel('uV')
axis([0, 60, 0, 1])
hold

Fs = 5000;
dt = 1/Fs;
N = size(chanAvg, 2);
dF = Fs/N;
f = -Fs/2:dF:Fs/2 - dF;

for chanIdx = 1:size(EpochAvg, 1)
    chanfft = fft(EpochAvg(chanIdx,:));
    chanfft = fftshift(chanfft);
    %plot from 0 to 60 Hz, 
    %NOTE: since fft is centered, 0 Hz = N/2 samples & bins are 0.5Hz
    plot(f(N/2 + 1:N/2+(60*2)),abs(chanfft(N/2 + 1:N/2 + (60*2)))/N)
end

subplot(2,2,4)
title('Frequency')
xlabel('Frequency(Hz)')
ylabel('uV')
axis([0, 20, 0, 1])
hold
for chanIdx = 1:size(EpochAvg, 1)
    chanfft = fft(EpochAvg(chanIdx,:));
    chanfft = fftshift(chanfft);
    %plot from 0 to 60 Hz, 
    %NOTE: since fft is centered, 0 Hz = N/2 samples & bins are 0.5Hz
    plot(f(N/2 + 1:N/2+(60*2)),abs(chanfft(N/2 + 1:N/2 + (60*2)))/N)
end
%% Get figures for rest inside scanner

figure
subplot(2,2,1)
title('Rest Outside Scanner')
xlabel('time(s)')
ylabel('uV')
axis([0, 2, -8, 8])
hold

%find average of epochs for each channel
for chanIdx = 1:size(SSVEP_Chan, 2)
    sumEpochs = Epochs_Rest.data(SSVEP_Chan(chanIdx),:,1);
    for epochIdx = 2:numEpochs
        sumEpochs = sumEpochs +  ...
            Epochs_Rest.data(SSVEP_Chan(chanIdx), :, epochIdx);
    end
    
    EpochAvg(chanIdx, :) = sumEpochs ./ numEpochs;
    
    %plot each Channel average
    t = 0:1/5000:2-(1/5000);
    plot(t, EpochAvg(chanIdx, :))

    
end


%find average of each channel
sumChans = EpochAvg(1, :);
for chanIdx = 2:size(EpochAvg, 1)
    sumChans = sumChans + EpochAvg(chanIdx, :);
end
subplot(2,2,2)
title('Average of O1, O2, P7, P8, P1, P2, PO3, PO4, PO7, PO8')
xlabel('time(s)')
ylabel('uV')
axis([0, 2, -8, 8])
chanAvg = sumChans ./ size(EpochAvg, 1);
plot(t, chanAvg)

%plot freqency domain of average chanAvg and of multiple epochs
subplot(2,2,3)
title('Frequency')
xlabel('Frequency(Hz)')
ylabel('uV')
axis([0, 60, 0, 1])
hold

Fs = 5000;
dt = 1/Fs;
N = size(chanAvg, 2);
dF = Fs/N;
f = -Fs/2:dF:Fs/2 - dF;

for chanIdx = 1:size(EpochAvg, 1)
    chanfft = fft(EpochAvg(chanIdx,:));
    chanfft = fftshift(chanfft);
    %plot from 0 to 60 Hz, 
    %NOTE: since fft is centered, 0 Hz = N/2 samples & bins are 0.5Hz
    plot(f(N/2 + 1:N/2+(60*2)),abs(chanfft(N/2 + 1:N/2 + (60*2)))/N)
end

subplot(2,2,4)
title('Frequency')
xlabel('Frequency(Hz)')
ylabel('uV')
axis([0, 20, 0, 1])
hold
for chanIdx = 1:size(EpochAvg, 1)
    chanfft = fft(EpochAvg(chanIdx,:));
    chanfft = fftshift(chanfft);
    %plot from 0 to 60 Hz, 
    %NOTE: since fft is centered, 0 Hz = N/2 samples & bins are 0.5Hz
    plot(f(N/2 + 1:N/2+(60*2)),abs(chanfft(N/2 + 1:N/2 + (60*2)))/N)
end
