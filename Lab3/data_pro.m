clc;
clear;
close all;
%% variable constants:
% Radar and acquisition parameters.
c = physconst('LightSpeed');
eps = 1e-5;
N_plots = 6;

%% reading files
[file,path] = uigetfile("*.bin","Select the data file","MultiSelect", "on");
%file = char(file_string)
%path = 'C:\Users\User\OneDrive\Documenten\MATLAB\lab3_files\'


for i = 1:length(file)
    %% detrending
    file(i)
    [BScans,SweepTime,CentreFreq,Bandwidth,T] = ReadSiversBIN(char(file(i)),path);
    % setting radar parameters from read file
    lambda = c / CentreFreq;
    [Nsamp, NSweeps] = size(BScans);
    Fs = Nsamp / SweepTime;
    fastTime = (0:Nsamp-1).' / Fs;
    slowTime = (0:NSweeps-1) * (T / NSweeps);

    % detrend
    BScans_detrend = detrend(BScans,'linear');

    % plot raw data and detrended
    figure('Name', 'Raw And Detrended Sweep', 'Color', 'w');
    subplot(N_plots, 1, 1);
    plot(fastTime * 1e3, BScans(:, 1), 'LineWidth', 1.0);
    grid on;
    xlabel('Fast time (ms)');
    ylabel('Amplitude');
    title('First Raw Sweep');
    
    subplot(N_plots, 1, 2);
    plot(fastTime * 1e3, BScans_detrend(:, 1), 'LineWidth', 1.0);
    grid on;
    xlabel('Fast time (ms)');
    ylabel('Amplitude');
    title('First Sweep After DC / Trend Removal');

    win = hamming(Nsamp);
    x = BScans_detrend .* win;
%% FMCW radar
    % to get the range profile use fft over fast time
    if Bandwidth>0 % check if CW or FMCW
        % range detection:
        x = x - mean(x, 2);
        Y = fft(x, [], 1);
        P2 = abs(Y/Nsamp);
        P1 = P2(1:Nsamp/2+1, :);
        P1(2:end-1, :) = 2*P1(2:end-1, :);
        
        f = Fs*(0:(Nsamp/2))/Nsamp;
        range = (c * f) / ((2 * (Bandwidth / SweepTime))+eps);
        % range bigger then 1 and smaller then 6:
        rangeWindow = find(range>1 & range<6);
        % Find max range bin per sweep (column-wise)
        [~, r_idx_range] = max(P1(rangeWindow,:), [], 1);   % 1 × NSweeps
        
        % Convert indices to actual range values
        range_max = range(rangeWindow(r_idx_range));      % 1 × NSweeps
        range_axis_max = 6; %(c * Fs) / ((2 * (Bandwidth / SweepTime))+eps);
        idx = range <= range_axis_max;
        % plotting
        subplot(N_plots, 1, 3);
        imagesc(slowTime, range(idx), 20*log10(P1(idx,:)+eps))
        grid on;
        axis xy
        xlabel('Slow Time (s)')
        ylabel('Range (m)')
        title('Range-Time Map')
        colorbar
        subplot(N_plots, 1, 4);
        plot(slowTime, range_max)
        grid on;
        axis xy
        xlabel('Slow Time (s)')
        ylabel('Range (m)')
        title('Range-Time Map')

        % CW radar breaks here, but can measure velocity over slow time
        % FMCW shows range correctly, can improve result by removing static
        % objects
        % There is also the static tests, which show distance
        
        % velocity detection:
        Y_range = fft(BScans_detrend, [], 1);
        Y_range = Y_range - mean(Y_range, 2);
        Y_range = Y_range(1:Nsamp/2+1, :);
        
        RD = fft(Y_range, [], 2);
        %RD = fftshift(RD, 2);
        PRF = NSweeps / T;   % chirp rate
        
        % --- pick strongest range bin ---
        RD_mag = abs(RD);

        % Only consider selected range window
        RD_window = RD_mag(rangeWindow, :);
        
        % Find strongest range bin inside window
        [~, r_idx_local] = max(sum(RD_window, 2));
        
        % Map back to full range index
        range_idx = find(rangeWindow);
        r_idx = range_idx(r_idx_local);
        
        % Extract slow-time signal at that range
        slow_signal = mean(Y_range(rangeWindow, :), 1);
        
        % --- STFT parameters ---
        winLen = 128;
        hop = 2;
        nFrames = floor((NSweeps - winLen)/hop) + 1;
        
        vel_time = zeros(1, nFrames);
        time_axis = zeros(1, nFrames);
        RD_time = zeros(winLen, nFrames);

        for k = 1:nFrames
            idx_win = (k-1)*hop + (1:winLen);
            segment = slow_signal(idx_win);
        
            segment = segment .* hamming(winLen).';
        
            spec = fft(segment);
            mag = abs(spec)/winLen;
        
            RD_time(:,k) = mag;   % store full spectrum
        
            f_doppler_win = (0:floor(winLen)-1) * (PRF / winLen);
            vel_axis = (f_doppler_win * lambda) / 2;
        
            [~, id] = max(mag);
            vel_time(k) = vel_axis(id);
        
            time_axis(k) = slowTime(idx_win(round(winLen/2)));
        end
        subplot(N_plots,1,5)
        imagesc(time_axis, vel_axis, 20*log10(RD_time + eps))
        axis xy
        xlabel('Time (s)')
        ylabel('Velocity (m/s)')
        title('Velocity-Time Map')
        colorbar
        subplot(N_plots,1,6)
        plot(time_axis, vel_time, 'LineWidth', 1)
        grid on
        xlabel('Time (s)')
        ylabel('Velocity (m/s)')
        title('Maximum Velocity vs Time')

    else
        % velocity estimation using fft over fast time
        cwSignal = BScans_detrend(:);
        cwSignal = cwSignal - mean(cwSignal);   % remove DC
        
        % --- STFT parameters ---
        winLen = 256;                 % window length (tune this!)
        hop = 32;                     % step size
        NfftDoppler = 512;            % FFT size
        
        nFrames = floor((length(cwSignal)-winLen)/hop) + 1;
        
        RD = zeros(NfftDoppler, nFrames);
        
        % --- Sliding FFT ---
        for k = 1:nFrames
            idx = (k-1)*hop + (1:winLen);
            segment = cwSignal(idx);
        
            % windowing
            segment = segment .* hann(winLen);
        
            % FFT
            spec = fft(segment, NfftDoppler);
        
            RD(:,k) = abs(spec);
        end
        % Doppler frequency axis
        half = 1:floor(NfftDoppler/2)+1;
        RD_oneSided = RD(half, :);
        fdAxis = (0:floor(NfftDoppler/2))' * (Fs / NfftDoppler);
        velAxis = fdAxis * lambda / 2;
        vel_max = 5;
        idx = abs(velAxis) <= vel_max;
        % Time axis (center of each window)
        timeAxis = ((0:nFrames-1)*hop + winLen/2) / Fs;
        subplot(N_plots,1,5)
        imagesc(timeAxis, velAxis(idx), 20*log10(RD(idx,:) + eps))
        axis xy
        xlabel('Time (s)')
        ylabel('Velocity (m/s)')
        title('Velocity-Time Map (CW Doppler)')
        colorbar
        subplot(N_plots,1,6)
        imagesc(timeAxis, velAxis(idx), 20*log10(RD(idx,:) + eps))
        axis xy
        xlabel('Time (s)')
        ylabel('Velocity (m/s)')
        title('Velocity-Time Map (CW Doppler)')
        colorbar
    end
end