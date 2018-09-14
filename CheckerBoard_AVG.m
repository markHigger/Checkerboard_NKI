
fileIdx = 1;

%server processed
inputDir = '~/Documents/NKI/EEG_Data/EEG_Raw/Server_Processed/';

inputName = 'EEG_fMRI_20180830_01_Checkerboard_Flash_Inside_bcg.set';
files(fileIdx) = {[inputDir, inputName]};
Headers(fileIdx, 1) = {'Check 01 Server'};
Headers(fileIdx, 2) = {'Rest 01 Server'};
fileIdx = fileIdx + 1;

inputName = 'EEG_fMRI_20180830_03_Checkerboard_Flash_Inside_bcg.set';
files(fileIdx) = {[inputDir, inputName]};
Headers(fileIdx, 1) = {'Check 03 Server'};
Headers(fileIdx, 2) = {'Rest 03 Server'};
fileIdx = fileIdx + 1;

inputName = 'EEG_fMRI_20180830_04_Checkerboard_Flash_Inside_bcg.set';
files(fileIdx) = {[inputDir, inputName]};
Headers(fileIdx, 1) = {'Check 04 Server'};
Headers(fileIdx, 2) = {'Rest 04 Server'};
fileIdx = fileIdx + 1;

inputName = 'EEG_fMRI_20180830_06_Checkerboard_Flash_Inside_bcg.set';
files(fileIdx) = {[inputDir, inputName]};
Headers(fileIdx, 1) = {'Check 06 Server'};
Headers(fileIdx, 2) = {'Rest 06 Server'};
fileIdx = fileIdx + 1;

%***Bad Data
%inputName = 'EEG_fMRI_20180830_Checkerboard_Flash_Outside_bcg.set';

inputName = 'EEG_fMRI_20180830_Checkerboard_Flash_Outside_02_bcg.set';
files(fileIdx) = {[inputDir, inputName]};
Headers(fileIdx, 1) = {'Check OS Server'};
Headers(fileIdx, 2) = {'Rest  OS Server'};
fileIdx = fileIdx + 1;

%Locally processed
inputDir = '~/Documents/NKI/EEG_Data/Processed_Data/';

inputName = 'EEG_fMRI_20180830_01_Checkerboard_Flash_Inside_04_bcg_bp.set';
files(fileIdx) = {[inputDir, inputName]};
Headers(fileIdx, 1) = {'Check 01 Local'};
Headers(fileIdx, 2) = {'Rest 01 Local'};
fileIdx = fileIdx + 1;

inputName ='EEG_fMRI_20180830_Checkerboard_Flash_Outside_02_bcg_bp.set';
files(fileIdx) = {[inputDir, inputName]};
Headers(fileIdx, 1) = {'Check OS Local'};
Headers(fileIdx, 2) = {'Rest OS Local'};
fileIdx = fileIdx + 1;

%% set Channels and Axis
%             P3, P4, O1, O2, P7, P8, Pz, Oz, POz, P1, P2, PO3, PO4, P5, P6, PO7, PO8 
SSVEP_Chan = [ 7,  8,  9, 10, 15, 16, 19, 20,  31, 37, 38,  45,  46, 51, 52,  59,  60];
timeAxis = [0, 2, -10, 10];
freqAxis_big = [0, 60, 0, 1.5];
freqAxis_small = [0, 20, 0, 1.5];

for idx = 1:fileIdx - 1
    %load in EEG
    EEG_Processed = pop_loadset(files{idx});

    %Bandpass bcg data
    EEG_Processed = EEG_Bandpass_Matlab(EEG_Processed, 0.5, 70, 2);

    %find checkerboard and rest epochs
    [Epochs_Check, Epochs_Rest, EEG_events] = ...
        EEG_Epoch_checkNrest(EEG_Processed, 'S  1', 'Check', 'Rest');
    
    plotChans(Epochs_Check, SSVEP_Chan, Headers{idx, 1}, timeAxis, freqAxis_big, freqAxis_small)
    
    plotChans(Epochs_Rest, SSVEP_Chan, Headers{idx, 2}, timeAxis, freqAxis_big, freqAxis_small)
end

%Get channels with ssvep response for average
 


%plotChans(Epochs, Chans, Header, timeAxis, freqAxis_big, freqAxis_small)


%% plot functions 
function plotChans(Epochs, Chans, Header, timeAxis, freqAxis_big, freqAxis_small)
    
    figure()
    numEpochs = size(Epochs.data,3);

    subplot(2,2,1)
    title(Header)
    xlabel('time(s)')
    ylabel('uV')
    axis(timeAxis)
    hold

    %find average of epochs for each channel
    for chanIdx = 1:size(Chans, 2)
        sumEpochs = Epochs.data(Chans(chanIdx),:,1);
        for epochIdx = 2:numEpochs
            sumEpochs = sumEpochs +  ...
                Epochs.data(Chans(chanIdx), :, epochIdx);
        end

        EpochAvg(chanIdx, :) = sumEpochs ./ numEpochs;

        %plot each Channel average
        t = 0:1/5000:2-(1/5000);
        plot(t, EpochAvg(chanIdx, :))


    end
    
    sumChans = EpochAvg(1, :);
    for chanIdx = 2:size(EpochAvg, 1)
        sumChans = sumChans + EpochAvg(chanIdx, :);
    end
    subplot(2,2,2)
    title('Average of Occipital & Perital')
    xlabel('time(s)')
    ylabel('uV')
    axis(timeAxis)
    hold
    chanAvg = sumChans ./ size(EpochAvg, 1);
    plot(t, chanAvg)
    
    %plot freqency domain of average chanAvg and of multiple epochs
    subplot(2,2,3)
    title('Frequency')
    xlabel('Frequency(Hz)')
    ylabel('uV')
    axis(freqAxis_big)
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
    axis(freqAxis_small)
    hold
    for chanIdx = 1:size(EpochAvg, 1)
        chanfft = fft(EpochAvg(chanIdx,:));
        chanfft = fftshift(chanfft);
        %plot from 0 to 60 Hz, 
        %NOTE: since fft is centered, 0 Hz = N/2 samples & bins are 0.5Hz
        plot(f(N/2 + 1:N/2+(60*2)),abs(chanfft(N/2 + 1:N/2 + (60*2)))/N)
    end
end
%% Legacy code

% %% make figures for checkboard inside scanner
% figure
% subplot(2,2,1)
% title(Header)
% xlabel('time(s)')
% ylabel('uV')
% axis(timeAxis)
% hold
% 
% %find average of epochs for each channel
% for chanIdx = 1:size(SSVEP_Chan, 2)
%     sumEpochs = Epochs_Check.data(SSVEP_Chan(chanIdx),:,1);
%     for epochIdx = 2:numEpochs
%         sumEpochs = sumEpochs +  ...
%             Epochs_Check.data(SSVEP_Chan(chanIdx), :, epochIdx);
%     end
%     
%     EpochAvg(chanIdx, :) = sumEpochs ./ numEpochs;
%     
%     %plot each Channel average
%     t = 0:1/5000:2-(1/5000);
%     plot(t, EpochAvg(chanIdx, :))
% 
%     
% end
% 
% 
% %find average of each channel
% sumChans = EpochAvg(1, :);
% for chanIdx = 2:size(EpochAvg, 1)
%     sumChans = sumChans + EpochAvg(chanIdx, :);
% end
% subplot(2,2,2)
% title('Average of Occipital & Perital')
% xlabel('time(s)')
% ylabel('uV')
% axis(timeAxis)
% hold
% chanAvg = sumChans ./ size(EpochAvg, 1);
% plot(t, chanAvg)
% 
% %plot freqency domain of average chanAvg and of multiple epochs
% subplot(2,2,3)
% title('Frequency')
% xlabel('Frequency(Hz)')
% ylabel('uV')
% axis(freqAxis_big)
% hold
% 
% Fs = 5000;
% dt = 1/Fs;
% N = size(chanAvg, 2);
% dF = Fs/N;
% f = -Fs/2:dF:Fs/2 - dF;
% 
% for chanIdx = 1:size(EpochAvg, 1)
%     chanfft = fft(EpochAvg(chanIdx,:));
%     chanfft = fftshift(chanfft);
%     %plot from 0 to 60 Hz, 
%     %NOTE: since fft is centered, 0 Hz = N/2 samples & bins are 0.5Hz
%     plot(f(N/2 + 1:N/2+(60*2)),abs(chanfft(N/2 + 1:N/2 + (60*2)))/N)
% end
% 
% subplot(2,2,4)
% title('Frequency')
% xlabel('Frequency(Hz)')
% ylabel('uV')
% axis(freqAxis_small)
% hold
% for chanIdx = 1:size(EpochAvg, 1)
%     chanfft = fft(EpochAvg(chanIdx,:));
%     chanfft = fftshift(chanfft);
%     %plot from 0 to 60 Hz, 
%     %NOTE: since fft is centered, 0 Hz = N/2 samples & bins are 0.5Hz
%     plot(f(N/2 + 1:N/2+(60*2)),abs(chanfft(N/2 + 1:N/2 + (60*2)))/N)
% end
%% Get figures for rest inside scanner

% figure
% subplot(2,2,1)
% title('Rest Inside Scanner')
% xlabel('time(s)')
% ylabel('uV')
% axis(timeAxis)
% hold
% 
% %find average of epochs for each channel
% for chanIdx = 1:size(SSVEP_Chan, 2)
%     sumEpochs = Epochs_Rest.data(SSVEP_Chan(chanIdx),:,1);
%     for epochIdx = 2:numEpochs
%         sumEpochs = sumEpochs +  ...
%             Epochs_Rest.data(SSVEP_Chan(chanIdx), :, epochIdx);
%     end
%     
%     EpochAvg(chanIdx, :) = sumEpochs ./ numEpochs;
%     
%     %plot each Channel average
%     t = 0:1/5000:2-(1/5000);
%     plot(t, EpochAvg(chanIdx, :))
% 
%     
% end
% 
% 
% %find average of each channel
% sumChans = EpochAvg(1, :);
% for chanIdx = 2:size(EpochAvg, 1)
%     sumChans = sumChans + EpochAvg(chanIdx, :);
% end
% subplot(2,2,2)
% title('Average of O1, O2, P7, P8, P1, P2, PO3, PO4, PO7, PO8')
% xlabel('time(s)')
% ylabel('uV')
% axis(timeAxis)
% hold
% chanAvg = sumChans ./ size(EpochAvg, 1);
% plot(t, chanAvg)
% 
% %plot freqency domain of average chanAvg and of multiple epochs
% subplot(2,2,3)
% title('Frequency')
% xlabel('Frequency(Hz)')
% ylabel('uV')
% axis(freqAxis_big)
% hold
% 
% Fs = 5000;
% dt = 1/Fs;
% N = size(chanAvg, 2);
% dF = Fs/N;
% f = -Fs/2:dF:Fs/2 - dF;
% 
% for chanIdx = 1:size(EpochAvg, 1)
%     chanfft = fft(EpochAvg(chanIdx,:));
%     chanfft = fftshift(chanfft);
%     %plot from 0 to 60 Hz, 
%     %NOTE: since fft is centered, 0 Hz = N/2 samples & bins are 0.5Hz
%     plot(f(N/2 + 1:N/2+(60*2)),abs(chanfft(N/2 + 1:N/2 + (60*2)))/N)
% end
% 
% subplot(2,2,4)
% title('Frequency')
% xlabel('Frequency(Hz)')
% ylabel('uV')
% axis(freqAxis_small)
% hold
% for chanIdx = 1:size(EpochAvg, 1)
%     chanfft = fft(EpochAvg(chanIdx,:));
%     chanfft = fftshift(chanfft);
%     %plot from 0 to 60 Hz, 
%     %NOTE: since fft is centered, 0 Hz = N/2 samples & bins are 0.5Hz
%     plot(f(N/2 + 1:N/2+(60*2)),abs(chanfft(N/2 + 1:N/2 + (60*2)))/N)
% end




% %% get epochs for outside scanner
% inputName = 'EEG_fMRI_20180830_Checkerboard_Flash_Outside_02_bcg_bp.set';
% fullfile_input = [inputDir, inputName];
% 
% EEG_Processed = pop_loadset(fullfile_input);
% 
% [Epochs_Check, Epochs_Rest, EEG_events] = ...
%     EEG_Epoch_checkNrest(EEG_Processed, 'S  1', 'Check', 'Rest');
% 
% 
% numEpochs = size(Epochs_Check.data,3);
% 
% %Get channels with ssvep response for average
% %             O1, O2, P7, P8, P1, P2, PO3, PO4, PO7, PO8 
% SSVEP_Chan = [ 9, 10, 15, 16, 37, 38,  45,  46,  59,  60]; 
% 
% %% make figures for checkboard outside scanner
% figure
% subplot(2,2,1)
% title('Checkerboard Outside Scanner')
% xlabel('time(s)')
% ylabel('uV')
% axis([0, 2, -8, 8])
% hold
% 
% %find average of epochs for each channel
% for chanIdx = 1:size(SSVEP_Chan, 2)
%     sumEpochs = Epochs_Check.data(SSVEP_Chan(chanIdx),:,1);
%     for epochIdx = 2:numEpochs
%         sumEpochs = sumEpochs +  ...
%             Epochs_Check.data(SSVEP_Chan(chanIdx), :, epochIdx);
%     end
%     
%     EpochAvg(chanIdx, :) = sumEpochs ./ numEpochs;
%     
%     %plot each Channel average
%     t = 0:1/5000:2-(1/5000);
%     plot(t, EpochAvg(chanIdx, :))
% 
%     
% end
% 
% 
% %find average of each channel
% sumChans = EpochAvg(1, :);
% for chanIdx = 2:size(EpochAvg, 1)
%     sumChans = sumChans + EpochAvg(chanIdx, :);
% end
% subplot(2,2,2)
% title('Average of O1, O2, P7, P8, P1, P2, PO3, PO4, PO7, PO8')
% xlabel('time(s)')
% ylabel('uV')
% axis([0, 2, -8, 8])
% hold
% chanAvg = sumChans ./ size(EpochAvg, 1);
% plot(t, chanAvg)
% 
% %plot freqency domain of average chanAvg and of multiple epochs
% subplot(2,2,3)
% title('Frequency')
% xlabel('Frequency(Hz)')
% ylabel('uV')
% axis([0, 60, 0, 1])
% hold
% 
% Fs = 5000;
% dt = 1/Fs;
% N = size(chanAvg, 2);
% dF = Fs/N;
% f = -Fs/2:dF:Fs/2 - dF;
% 
% for chanIdx = 1:size(EpochAvg, 1)
%     chanfft = fft(EpochAvg(chanIdx,:));
%     chanfft = fftshift(chanfft);
%     %plot from 0 to 60 Hz, 
%     %NOTE: since fft is centered, 0 Hz = N/2 samples & bins are 0.5Hz
%     plot(f(N/2 + 1:N/2+(60*2)),abs(chanfft(N/2 + 1:N/2 + (60*2)))/N)
% end
% 
% subplot(2,2,4)
% title('Frequency')
% xlabel('Frequency(Hz)')
% ylabel('uV')
% axis([0, 20, 0, 1])
% hold
% for chanIdx = 1:size(EpochAvg, 1)
%     chanfft = fft(EpochAvg(chanIdx,:));
%     chanfft = fftshift(chanfft);
%     %plot from 0 to 60 Hz, 
%     %NOTE: since fft is centered, 0 Hz = N/2 samples & bins are 0.5Hz
%     plot(f(N/2 + 1:N/2+(60*2)),abs(chanfft(N/2 + 1:N/2 + (60*2)))/N)
% end
% %% Get figures for rest outside scanner
% 
% figure
% subplot(2,2,1)
% title('Rest Outside Scanner')
% xlabel('time(s)')
% ylabel('uV')
% axis([0, 2, -8, 8])
% hold
% 
% %find average of epochs for each channel
% for chanIdx = 1:size(SSVEP_Chan, 2)
%     sumEpochs = Epochs_Rest.data(SSVEP_Chan(chanIdx),:,1);
%     for epochIdx = 2:numEpochs
%         sumEpochs = sumEpochs +  ...
%             Epochs_Rest.data(SSVEP_Chan(chanIdx), :, epochIdx);
%     end
%     
%     EpochAvg(chanIdx, :) = sumEpochs ./ numEpochs;
%     
%     %plot each Channel average
%     t = 0:1/5000:2-(1/5000);
%     plot(t, EpochAvg(chanIdx, :))
% 
%     
% end
% 
% 
% %find average of each channel
% sumChans = EpochAvg(1, :);
% for chanIdx = 2:size(EpochAvg, 1)
%     sumChans = sumChans + EpochAvg(chanIdx, :);
% end
% subplot(2,2,2)
% title('Average of O1, O2, P7, P8, P1, P2, PO3, PO4, PO7, PO8')
% xlabel('time(s)')
% ylabel('uV')
% axis([0, 2, -8, 8])
% chanAvg = sumChans ./ size(EpochAvg, 1);
% plot(t, chanAvg)
% 
% %plot freqency domain of average chanAvg and of multiple epochs
% subplot(2,2,3)
% title('Frequency')
% xlabel('Frequency(Hz)')
% ylabel('uV')
% axis([0, 60, 0, 1])
% hold
% 
% Fs = 5000;
% dt = 1/Fs;
% N = size(chanAvg, 2);
% dF = Fs/N;
% f = -Fs/2:dF:Fs/2 - dF;
% 
% for chanIdx = 1:size(EpochAvg, 1)
%     chanfft = fft(EpochAvg(chanIdx,:));
%     chanfft = fftshift(chanfft);
%     %plot from 0 to 60 Hz, 
%     %NOTE: since fft is centered, 0 Hz = N/2 samples & bins are 0.5Hz
%     plot(f(N/2 + 1:N/2+(60*2)),abs(chanfft(N/2 + 1:N/2 + (60*2)))/N)
% end
% 
% subplot(2,2,4)
% title('Frequency')
% xlabel('Frequency(Hz)')
% ylabel('uV')
% axis([0, 20, 0, 1])
% hold
% for chanIdx = 1:size(EpochAvg, 1)
%     chanfft = fft(EpochAvg(chanIdx,:));
%     chanfft = fftshift(chanfft);
%     %plot from 0 to 60 Hz, 
%     %NOTE: since fft is centered, 0 Hz = N/2 samples & bins are 0.5Hz
%     plot(f(N/2 + 1:N/2+(60*2)),abs(chanfft(N/2 + 1:N/2 + (60*2)))/N)
% end
