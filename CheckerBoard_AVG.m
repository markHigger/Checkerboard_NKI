
fileIdx = 1;

%             P3, P4, O1, O2, P7, P8, Pz, Oz, POz, P1, P2, PO3, PO4, P5, P6, PO7, PO8 
SSVEP_Chan = [ 7,  8,  9, 10, 15, 16, 19, 20,  31, 37, 38,  45,  46, 51, 52,  59,  60];

%server processed
inputDir = '~/Documents/NKI/EEG_Data/EEG_Raw/Server_Processed/';

inputName = 'EEG_fMRI_20180830_01_Checkerboard_Flash_Inside_bcg.set';
files(fileIdx) = {[inputDir, inputName]};
Headers(fileIdx, 1, :) = {'Checkerboard | Run 01 ', 'Epoch: 2s | Analysis: Server'};
Headers(fileIdx, 2, :) = {'Rest | Run 01', 'Epoch: 2s | Analysis: Server'};
fileDesc{fileIdx} = 'Run_01';
fileIdx = fileIdx + 1;

 inputName = 'EEG_fMRI_20180830_03_Checkerboard_Flash_Inside_bcg.set';
 files(fileIdx) = {[inputDir, inputName]};
Headers(fileIdx, 1, :) = {'Checkerboard | Run 03 ', 'Epoch: 2s | Analysis: Server'};
Headers(fileIdx, 2, :) = {'Rest | Run 03', 'Epoch: 2s | Analysis: Server'};
fileDesc{fileIdx} = 'Run_03';
 fileIdx = fileIdx + 1;

inputName = 'EEG_fMRI_20180830_04_Checkerboard_Flash_Inside_bcg.set';
files(fileIdx) = {[inputDir, inputName]};
Headers(fileIdx, 1, :) = {'Checkerboard | Run 04 ', 'Epoch: 2s | Analysis: Server'};
Headers(fileIdx, 2, :) = {'Rest | Run 04', 'Epoch: 2s | Analysis: Server'};
fileDesc{fileIdx} = 'Run_04';
fileIdx = fileIdx + 1;


inputName = 'EEG_fMRI_20180830_06_Checkerboard_Flash_Inside_bcg.set';
files(fileIdx) = {[inputDir, inputName]};
Headers(fileIdx, 1, :) = {'Checkerboard | Run 06 ', 'Epoch: 2s | Analysis: Server'};
Headers(fileIdx, 2, :) = {'Rest | Run 06', 'Epoch: 2s | Analysis: Server'};
fileDesc{fileIdx} = 'Run_06';
fileIdx = fileIdx + 1;

% % %***Bad Data
% % %inputName = 'EEG_fMRI_20180830_Checkerboard_Flash_Outside_bcg.set';
% 
inputName = 'EEG_fMRI_20180830_Checkerboard_Flash_Outside_02_bcg.set';
files(fileIdx) = {[inputDir, inputName]};
Headers(fileIdx, 1, :) = {'Checkerboard | Outside Scanner ', 'Epoch: 2s | Analysis: Server'};
Headers(fileIdx, 2, :) = {'Rest | Outside Scanner', 'Epoch: 2s | Analysis: Server'};
fileDesc{fileIdx} = 'Run_Outside';
fileIdx = fileIdx + 1;

%Locally processed
inputDir = '~/Documents/NKI/EEG_Data/Processed_Data/';

inputName = 'EEG_fMRI_20180830_01_Checkerboard_Flash_Inside_04_bcg_bp.set';
files(fileIdx) = {[inputDir, inputName]};
Headers(fileIdx, 1, :) = {'Checkerboard | Run 1 ', 'Epoch: 2s | Analysis: Local'};
Headers(fileIdx, 2, :) = {'Rest | Run 01', 'Epoch: 2s | Analysis: Local'};
fileDesc{fileIdx} = 'Run_01_local';
fileIdx = fileIdx + 1;

inputName ='EEG_fMRI_20180830_Checkerboard_Flash_Outside_02_bcg_bp.set';
files(fileIdx) = {[inputDir, inputName]};
Headers(fileIdx, 1, :) = {'Checkerboard | Out of Scanner', 'Epoch: 2s | Analysis: Local'};
Headers(fileIdx, 2, :) = {'Rest | Outside Scanner', 'Epoch: 2s | Analysis: Local'};
fileDesc{fileIdx} = 'Run_Outside_Local';
fileIdx = fileIdx + 1;

%% set Channels and Axis

timeAxis = [0, 2, -10, 10];
timeAxis_STD = [0, 2, -20, 20];
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
    
    %plotTimeFreq(Epochs_Check, Epochs_Rest, SSVEP_Chan, Headers(idx, :, :), timeAxis, freqAxis_big, freqAxis_small, 1);
    %plotFreqs(Epochs_Check, Epochs_Rest, SSVEP_Chan, Headers(idx, :, :), freqAxis_big, freqAxis_small, 1);
    plotLocs(Epochs_Check, Epochs_Rest, Headers(idx, :, :), 1, fileDesc{idx})
    %plotAllMean(Epochs_Check, Epochs_Rest, Headers(idx, :, :), timeAxis, timeAxis_STD, 1, fileDesc{idx})
    %close all
    %plotLocFreqs(Epochs_Check, Epochs_Rest, Headers(idx, :, :), freqAxis_big, 1, fileDesc{idx})
end

%% plot functions 
function plotTimeFreq(checkEpochs, restEpochs, Chans, Headers, timeAxis, freqAxis_big, freqAxis_small, saveFigs)
    
%% Plot Time Axis for Checkerboard
    fig = figure();
    numEpochs = size(restEpochs.data,3);
    %figure()
    subplot(2,2,1)
    title(Headers(:,1,:))
    xlabel('time(s)')
    ylabel('uV')
    axis(timeAxis)
    hold

    %find Time axis for Check
    for chanIdx = 1:size(Chans, 2)
        EpochDataT(:,:) = checkEpochs.data(Chans(chanIdx),:,:);
        EpochData = EpochDataT';
        EpochAvg_Check(chanIdx, :) = mean(EpochData);
        %Mean stadard average = sig/sqrt(N)
        EpochMSA_Check(chanIdx, :) = std(EpochData) / sqrt(length(EpochData));
        %plot each Channel average
        t = 0:1/5000:2-(1/5000);
        plot(t, EpochAvg_Check(chanIdx, :))
        %errorbar(t, EpochAvg_Check(chanIdx, :), EpochMSA_Check(chanIdx,:))
        
        
        
    end
    
%% Find time axis for rest
    %find time axis of Rest
    subplot(2,2,2)
    title(Headers(:,2,:))
    xlabel('time(s)')
    ylabel('uV')
    axis(timeAxis)
    hold

    %find average of epochs for each channel
    for chanIdx = 1:size(Chans, 2)
        EpochDataT(:,:) = restEpochs.data(Chans(chanIdx),:,:);
        EpochData = EpochDataT';
        EpochAvg_Rest(chanIdx, :) = mean(EpochData);
        %Mean stadard average = sig/sqrt(N)
        EpochMSA_Rest(chanIdx, :) = std(EpochData) / sqrt(length(EpochData));
        %plot each Channel average
        t = 0:1/5000:2-(1/5000);
        plot(t, EpochAvg_Rest(chanIdx, :));
        %errorbar(t, EpochAvg_Rest(chanIdx, :), EpochMSA_Rest(chanIdx,:))
        
        
    end
  %     fprintf('%s: Max MSA = %d \n', Header, max(max(EpochMSA)));  
    
    Fs = 5000;
    N = size(EpochAvg_Check, 2);
    dF = Fs/N;
    f = -Fs/2:dF:Fs/2 - dF;


%     %plot freqency domain of average chanAvg and of multiple epochs
%     subplot(2,2,3)
%     title('Frequency')
%     xlabel('Frequency(Hz)')
%     ylabel('uV')
%     axis(freqAxis_big)
%     hold
% 

%     for chanIdx = 1:size(EpochAvg, 1)
%         chanfft = fft(EpochAvg(chanIdx,:));
%         chanfft = fftshift(chanfft);
%         %plot from 0 to 60 Hz, 
%         %NOTE: since fft is centered, 0 Hz = N/2 samples & bins are 0.5Hz
%         plot(f(N/2 + 1:N/2+(60*2)),abs(chanfft(N/2 + 1:N/2 + (60*2)))/N)
%     end

%% frequency for Check
    subplot(2,2,3)
    
    title('Frequency')
    xlabel('Frequency(Hz)')
    ylabel('uV')
    axis(freqAxis_small)
    hold
    for chanIdx = 1:size(EpochAvg_Check, 1)
        chanfft = fft(EpochAvg_Check(chanIdx,:));
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
    for chanIdx = 1:size(EpochAvg_Rest, 1)
        chanfft = fft(EpochAvg_Rest(chanIdx,:));
        chanfft = fftshift(chanfft);
        %plot from 0 to 60 Hz, 
        %NOTE: since fft is centered, 0 Hz = N/2 samples & bins are 0.5Hz
        plot(f(N/2 + 1:N/2+(60*2)),abs(chanfft(N/2 + 1:N/2 + (60*2)))/N)
    end
    str = [regexprep(Headers{1,1}, ' ', '_'), '.mat'];
    save(str, 'EpochAvg_Check', 'EpochMSA_Check', 'EpochAvg_Rest', 'EpochMSA_Rest');
    str = [regexprep(Headers{1,1}, ' ', '_'), '.png'];
    if saveFigs
        saveas(fig, str)
    end
end

function plotFreqs(checkEpochs, restEpochs, Chans, Headers, freqAxis_big, freqAxis_small, saveFigs)
%% Plot spectrum of channels
    for chanIdx = 1:length(Chans)
        checkDataT(:,:) = checkEpochs.data(Chans(chanIdx),:,:);
        checkData = checkDataT';
        restDataT(:,:) = restEpochs.data(Chans(chanIdx),:,:);
        restData = restDataT';
        EpochAvg_Check(chanIdx, :) = mean(checkData);
        EpochAvg_Rest(chanIdx, :) = mean(restData);
    end
    Fs = 5000;
    N = size(EpochAvg_Check, 2);
    dF = Fs/N;
    f = -Fs/2:dF:Fs/2 - dF;
    
    fig = figure();
    %plot freqency domain of average chanAvg and of multiple epochs
     subplot(2,2,1)
     title(Headers(:,1,:))
     xlabel('Frequency(Hz)')
     ylabel('uV')
     axis(freqAxis_big)
     hold
 

     for chanIdx = 1:size(EpochAvg_Check, 1)
         chanfft = fft(EpochAvg_Check(chanIdx,:));
         chanfft = fftshift(chanfft);
         %plot from 0 to 60 Hz, 
         %NOTE: since fft is centered, 0 Hz = N/2 samples & bins are 0.5Hz
         plot(f(N/2 + 1:N/2+(60*2)),abs(chanfft(N/2 + 1:N/2 + (60*2)))/N)
     end
     
     subplot(2,2,2)
     title('Frequency')
     xlabel('Frequency(Hz)')
     ylabel('uV')
     axis(freqAxis_small)
     hold
 

     for chanIdx = 1:size(EpochAvg_Check, 1)
         chanfft = fft(EpochAvg_Check(chanIdx,:));
         chanfft = fftshift(chanfft);
         %plot from 0 to 60 Hz, 
         %NOTE: since fft is centered, 0 Hz = N/2 samples & bins are 0.5Hz
         plot(f(N/2 + 1:N/2+(60*2)),abs(chanfft(N/2 + 1:N/2 + (60*2)))/N)
     end
     
     %plot freqency domain of average chanAvg and of multiple epochs
     subplot(2,2,3)
     title(Headers(:,2,:))
     xlabel('Frequency(Hz)')
     ylabel('uV')
     axis(freqAxis_big)
     hold
 

     for chanIdx = 1:size(EpochAvg_Check, 1)
         chanfft = fft(EpochAvg_Rest(chanIdx,:));
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
 

     for chanIdx = 1:size(EpochAvg_Check, 1)
         chanfft = fft(EpochAvg_Rest(chanIdx,:));
         chanfft = fftshift(chanfft);
         %plot from 0 to 60 Hz, 
         %NOTE: since fft is centered, 0 Hz = N/2 samples & bins are 0.5Hz
         plot(f(N/2 + 1:N/2+(60*2)),abs(chanfft(N/2 + 1:N/2 + (60*2)))/N)
     end
     
    str = [regexprep(Headers{1,1}, ' ', '_'), 'Freqs.png'];
    if saveFigs
        saveas(fig, str)
    end
     
     

end

function plotLocs(checkEpochs, restEpochs, Headers, saveFigs, fileDesc)
%% find average and standard error of check
    for chanIdx = 1:64
        EpochDataT(:,:) = checkEpochs.data(chanIdx,:,:);
        EpochData = EpochDataT';
        EpochAvg_Check(chanIdx, :) = mean(EpochData);
        EpochSTD_Check(chanIdx, :) = std(EpochData);
        %Mean stadard Error = sig/sqrt(N)
        EpochMSE_Check(chanIdx, :) = std(EpochData) / sqrt(length(EpochData));  
    end

%% put Average and avg +/- Error into eeglab topograph
    topoData(:, :, 1) = EpochAvg_Check;
    topoData(:, :, 2) = EpochAvg_Check + EpochMSE_Check;
    topoData(:, :, 3) = EpochAvg_Check - EpochMSE_Check;
    topoData = topoData / checkEpochs.srate;
    fig = figure;
    %get header into format readable by eeglab current format is ..
    %   {lin1Col1|lin1Col2}{lin2Col1|lin2Col2}... Must be in form
    %   lin1Col1|lin1Col2|lin1Col3|lin1Col4
    header = Headers{:,1,1};
    header = [header, '| Standard Error'];
    for headerCellIdx = 2:size(Headers, 3)
        header = [  header, ' | ', Headers{:,1,headerCellIdx}];
    end
    limits   = [0 2000 -0.004 0.004];
    plottopo(topoData, ...
        'chanlocs', restEpochs.chanlocs, ...
        'frames', 10000, ...
        'title', header, ...
        'limits', limits, ...
        'colors', {'k', 'r', 'r'});

    if saveFigs
        dir = '~/Documents/NKI/FreqsbyChan/';
        fileName = [fileDesc];
        filepng = [dir, fileName, 'Check_Error_Locs.png'];
        saveas(fig, filepng)
    end
%% put Average and avg +/- STD into eeglab topograph
    topoData(:, :, 1) = EpochAvg_Check;
    topoData(:, :, 2) = EpochAvg_Check + EpochSTD_Check;
    topoData(:, :, 3) = EpochAvg_Check - EpochSTD_Check;
    topoData = topoData / checkEpochs.srate;
    fig = figure;
    %get header into format readable by eeglab current format is ..
    %   {lin1Col1|lin1Col2}{lin2Col1|lin2Col2}... Must be in form
    %   lin1Col1|lin1Col2|lin1Col3|lin1Col4
    header = Headers{:,1,1};
    header = [header, '| Standard Deviation'];
    for headerCellIdx = 2:size(Headers, 3)
        header = [  header, ' | ', Headers{:,1,headerCellIdx}];
    end
    
    limits   = [0 2000 -0.004 0.004];
    plottopo(topoData, ...
        'chanlocs', restEpochs.chanlocs, ...
        'frames', 10000, ...
        'title', header, ...
        'limits', limits, ...
        'colors', {'k', 'r', 'r'});

    if saveFigs
        dir = '~/Documents/NKI/FreqsbyChan/';
        fileName = [fileDesc];
        filepng = [dir, fileName, 'Check_STD_Locs.png'];
        saveas(fig, filepng)
    end
 %% find average and standard error of Rest
    for chanIdx = 1:64
        EpochDataT(:,:) = restEpochs.data(chanIdx,:,:);
        EpochData = EpochDataT';
        EpochAvg_Rest(chanIdx, :) = mean(EpochData);
        EpochSTD_Rest(chanIdx, :) = std(EpochData);
        %Mean stadard Error = sig/sqrt(N)
        EpochMSE_Rest(chanIdx, :) = std(EpochData) / sqrt(length(EpochData));
    end

%% put Average and avg +/- Error into eeglab topograph
    topoData(:, :, 1) = EpochAvg_Rest;
    topoData(:, :, 2) = EpochAvg_Rest + EpochMSE_Rest;
    topoData(:, :, 3) = EpochAvg_Rest - EpochMSE_Rest;
    topoData = topoData / restEpochs.srate;
    fig = figure;
    %get header into format readable by eeglab current format is ..
    %   {lin1Col1|lin1Col2}{lin2Col1|lin2Col2}... Must be in form
    %   lin1Col1|lin1Col2|lin1Col3|lin1Col4
    header = Headers{:,2,1};
    header = [header, '| Standard Error'];
    for headerCellIdx = 2:size(Headers, 3)
        header = [  header, ' | ', Headers{:,2,headerCellIdx}];
    end
     limits   = [0 2000 -0.004 0.004];
    plottopo(topoData, ...
        'chanlocs', restEpochs.chanlocs, ...
        'frames', 10000, ...
        'title', header, ...
        'limits', limits, ...
        'colors', {'k', 'r', 'r'});


    if saveFigs
        dir = '~/Documents/NKI/FreqsbyChan/';
        fileName = [fileDesc];
        filepng = [dir, fileName, 'Rest_Error_Locs.png'];
        saveas(fig, filepng)
    end
    
%% put Average and avg +/- STD into eeglab topograph
    topoData(:, :, 1) = EpochAvg_Rest;
    topoData(:, :, 2) = EpochAvg_Rest + EpochSTD_Rest;
    topoData(:, :, 3) = EpochAvg_Rest - EpochSTD_Rest;
    topoData = topoData / restEpochs.srate;
    fig = figure;
    %get header into format readable by eeglab current format is ..
    %   {lin1Col1|lin1Col2}{lin2Col1|lin2Col2}... Must be in form
    %   lin1Col1|lin1Col2|lin1Col3|lin1Col4
    header = Headers{:,2,1};
    header = [header, '| Standard Deviation'];
    for headerCellIdx = 2:size(Headers, 3)
        header = [  header, ' | ', Headers{:,2,headerCellIdx}];
    end
    limits   = [0 2000 -0.004 0.004];
    plottopo(topoData, ...
        'chanlocs', restEpochs.chanlocs, ...
        'frames', 10000, ...
        'title', header, ...
        'limits', limits, ...
        'colors', {'k', 'r', 'r'});

    if saveFigs
        dir = '~/Documents/NKI/FreqsbyChan/';
        fileName = [fileDesc];
        filepng = [dir, fileName, 'Rest_STD_Locs.png'];
        saveas(fig, filepng)
    end
end

function plotAllMean(checkEpochs, restEpochs, Headers, timeAxis, timeAxis_STD, saveFigs, fileDesc)
%% plot means of all chans fo both check and rest

%plot error and std for checkerboard
    for chanIdx = 1:64
        %Transpose data to more workable format
        EpochDataT(:,:) = checkEpochs.data(chanIdx,:,:);
        EpochData = EpochDataT';
        
        %find mean, standard deviation and mean standard eror
        EpochAvg_Check(chanIdx, :) = mean(EpochData);
        EpochSTD_Check(chanIdx, :) = std(EpochData);
        %Mean stadard Error = sig/sqrt(N)
        EpochMSE_Check(chanIdx, :) = std(EpochData) / sqrt(length(EpochData));  
        
        %plot figure with Mean Error
        fig = figure();
      
        %Add Line indicating channel to plot
        chanHeader = Headers;
        %add line with 'Ch{channum}:{ChanName}'
        chanInfo = sprintf('Ch%d:%s | Mean Standard Error', chanIdx, checkEpochs.chanlocs(chanIdx).labels);
        chanHeader(:, :, size(Headers, 3)+1) = {chanInfo};
        title(chanHeader(:,1,:))
        
        xlabel('time(s)')
        ylabel('uV')
        axis(timeAxis)
        hold
        t = 0:1/5000:2-(1/5000);
        data = double(EpochAvg_Check(chanIdx, :));
        err = double(EpochMSE_Check(chanIdx, :));
        boundedline(t, data, err)
%         plot(t, EpochAvg_Check(chanIdx, :))
%         plot(t, EpochAvg_Check(chanIdx, :) - EpochMSE_Check(chanIdx, :))
%         plot(t, EpochAvg_Check(chanIdx, :) + EpochMSE_Check(chanIdx, :))
        
        if saveFigs
            dir = '~/Documents/NKI/MSE_Chan/';
            fileName = [fileDesc, '_' ,checkEpochs.chanlocs(chanIdx).labels];
            filepng = [dir, fileName, '_Check_MSE.png'];
            saveas(fig, filepng)
        end
        
        %plot figure with standard deviation
        fig = figure();
      
        %Add Line indicating channel to plot
        chanHeader = Headers;
        %add line with 'Ch{channum}:{ChanName}'
        chanInfo = sprintf('Ch%d:%s | Standard Deviation', chanIdx, checkEpochs.chanlocs(chanIdx).labels);
        chanHeader(:, :, size(Headers, 3)+1) = {chanInfo};
        title(chanHeader(:,1,:))
        
        xlabel('time(s)')
        ylabel('uV')
        axis(timeAxis_STD)
        hold
        t = 0:1/5000:2-(1/5000);
        data = double(EpochAvg_Check(chanIdx, :));
        err = double(EpochSTD_Check(chanIdx, :));
        boundedline(t, data, err)
%         plot(t, EpochAvg_Check(chanIdx, :))
%         plot(t, EpochAvg_Check(chanIdx, :) - EpochSTD_Check(chanIdx, :))
%         plot(t, EpochAvg_Check(chanIdx, :) + EpochSTD_Check(chanIdx, :))
        
        if saveFigs
            dir = '~/Documents/NKI/STD_Chan/';
            fileName = [fileDesc, '_' ,checkEpochs.chanlocs(chanIdx).labels];
            filepng = [dir, fileName, '_Check_STD.png'];
            saveas(fig, filepng)
        end
%% Find Mean at rest
        EpochDataT(:,:) = restEpochs.data(chanIdx,:,:);
        EpochData = EpochDataT';
        %find mean, standard deviation and mean standard eror
        EpochAvg_Rest(chanIdx, :) = mean(EpochData);
        EpochSTD_Rest(chanIdx, :) = std(EpochData);
        %Mean stadard Error = sig/sqrt(N)
        EpochMSE_Rest(chanIdx, :) = std(EpochData) / sqrt(length(EpochData));
    
        %plot figure with Mean Error
        fig = figure();
      
        %Add Line indicating channel to plot
        chanHeader = Headers;
        %add line with 'Ch{channum}:{ChanName}'
        chanInfo = sprintf('Ch%d:%s | Mean Standard Error', chanIdx, checkEpochs.chanlocs(chanIdx).labels);
        chanHeader(:, :, size(Headers, 3)+1) = {chanInfo};
        title(chanHeader(:,2,:))
        
        xlabel('time(s)')
        ylabel('uV')
        axis(timeAxis)
        hold
        t = 0:1/5000:2-(1/5000);
        data = double(EpochAvg_Rest(chanIdx, :));
        err = double(EpochMSE_Rest(chanIdx, :));
        boundedline(t, data, err)
%         plot(t, EpochAvg_Rest(chanIdx, :))
%         plot(t, EpochAvg_Rest(chanIdx, :) - EpochMSE_Rest(chanIdx, :))
%         plot(t, EpochAvg_Rest(chanIdx, :) + EpochMSE_Rest(chanIdx, :))
        
        if saveFigs
            dir = '~/Documents/NKI/MSE_Chan/';
            fileName = [fileDesc, '_' ,checkEpochs.chanlocs(chanIdx).labels];
            filepng = [dir, fileName, '_Rest_MSE.png'];
            saveas(fig, filepng)
        end
        
        
        %plot figure with standard deviation
        fig = figure();
      
        %Add Line indicating channel to plot
        chanHeader = Headers;
        %add line with 'Ch{channum}:{ChanName}'
        chanInfo = sprintf('Ch%d:%s | Standard Deviation', chanIdx, checkEpochs.chanlocs(chanIdx).labels);
        chanHeader(:, :, size(Headers, 3)+1) = {chanInfo};
        title(chanHeader(:,2,:))
        
        xlabel('time(s)')
        ylabel('uV')
        axis(timeAxis_STD)
        hold
        t = 0:1/5000:2-(1/5000);
        data = double(EpochAvg_Rest(chanIdx, :));
        err = double(EpochSTD_Rest(chanIdx, :));
        boundedline(t, data, err)
%         plot(t, EpochAvg_Rest(chanIdx, :))
%         plot(t, EpochAvg_Rest(chanIdx, :) - EpochSTD_Rest(chanIdx, :))
%         plot(t, EpochAvg_Rest(chanIdx, :) + EpochSTD_Rest(chanIdx, :))
        
        if saveFigs
            dir = '~/Documents/NKI/STD_Chan/';
            fileName = [fileDesc, '_' ,checkEpochs.chanlocs(chanIdx).labels];
            filepng = [dir, fileName, '_Rest_STD.png'];
            saveas(fig, filepng)
        end
    end
end

function plotLocFreqs(checkEpochs, restEpochs, Headers, freqAxis, saveFigs, fileDesc)
%% find average and standard error of check
    for chanIdx = 1:64
        EpochDataT(:,:) = checkEpochs.data(chanIdx,:,:);
        EpochData = EpochDataT';
        EpochAvg_Check(chanIdx, :) = mean(EpochData);
        Checkfft(chanIdx,:) = fft(EpochAvg_Check(chanIdx,:));
        Checkfft(chanIdx,:) = abs(fftshift(Checkfft(chanIdx,:)));
        %plot from 0 to 60 Hz, 
        %NOTE: since fft is centered, 0 Hz = N/2 samples & bins are 0.5Hz
 
    end

    %% Get Freq response for each channel
    
    Fs = 5000;
    N = size(EpochAvg_Check, 2);
    dF = Fs/N;
    f = -Fs/2:dF:Fs/2 - dF;
    
    topoData = Checkfft(:, N/2 + 1:N/2 + (60*2))/N;

    
    
    fig = figure;
    
    %get header into format readable by eeglab current format is ..
    %   {lin1Col1|lin1Col2}{lin2Col1|lin2Col2}... Must be in form
    %   lin1Col1|lin1Col2|lin1Col3|lin1Col4
    header = Headers{:,1,1};
    header = [header, '| Freq'];
    for headerCellIdx = 2:size(Headers, 3)
        header = [  header, ' | ', Headers{:,1,headerCellIdx}];
    end

    
    limits   = freqAxis;
    plottopo(topoData, ...
        'chanlocs', restEpochs.chanlocs, ...
        'frames', 120, ...
        'title', header, ...
        'limits', limits);

        if saveFigs
            dir = '~/Documents/NKI/FreqsbyChan/';
            fileName = [fileDesc];
            filepng = [dir, fileName, '_Check_FreqsbyChan.png'];
            saveas(fig, filepng)
        end
        
    for chanIdx = 1:64
        EpochDataT(:,:) = restEpochs.data(chanIdx,:,:);
        EpochData = EpochDataT';
        EpochAvg_Rest(chanIdx, :) = mean(EpochData);
        Restfft(chanIdx,:) = fft(EpochAvg_Rest(chanIdx,:));
        Restfft(chanIdx,:) = abs(fftshift(Restfft(chanIdx,:)));
        %plot from 0 to 60 Hz, 
        %NOTE: since fft is centered, 0 Hz = N/2 samples & bins are 0.5Hz
 
    end

    %% Get Freq response for each channel
    
    Fs = 5000;
    N = size(EpochAvg_Rest, 2);
    dF = Fs/N;
    f = -Fs/2:dF:Fs/2 - dF;
    
    topoData = Restfft(:, N/2 + 1:N/2 + (60*2))/N;

    
    
    fig = figure;
    
    %get header into format readable by eeglab current format is ..
    %   {lin1Col1|lin1Col2}{lin2Col1|lin2Col2}... Must be in form
    %   lin1Col1|lin1Col2|lin1Col3|lin1Col4
    header = Headers{:,2,1};
    header = [header, '| Freq'];
    for headerCellIdx = 2:size(Headers, 3)
        header = [  header, ' | ', Headers{:,2,headerCellIdx}];
    end

    
    limits   = freqAxis;
    plottopo(topoData, ...
        'chanlocs', restEpochs.chanlocs, ...
        'frames', 120, ...
        'title', header, ...
        'limits', limits);

        if saveFigs
            dir = '~/Documents/NKI/FreqsbyChan/';
            fileName = [fileDesc];
            filepng = [dir, fileName, '_Rest_FreqsbyChan.png'];
            saveas(fig, filepng)
        end
end