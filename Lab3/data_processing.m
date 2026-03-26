clc;
clear;
close all;

% Load one of the Lab 3 binary captures.
[BScans, SweepTime, CentreFreq, Bandwidth, T] = ReadSiversBIN;

% Radar and acquisition parameters.
c = physconst('LightSpeed');
lambda = c / CentreFreq;
[Nsamp, NSweeps] = size(BScans);
Fs = Nsamp / SweepTime;
fastTime = (0:Nsamp-1).' / Fs;
slowTime = (0:NSweeps-1) * (T / NSweeps);

% Remove the DC component and linear trend from each sweep.
rawData = double(BScans);
procData = detrend(rawData, 'linear');

fprintf('Samples per sweep : %d\n', Nsamp);
fprintf('Number of sweeps  : %d\n', NSweeps);
fprintf('Sweep time        : %.6f s\n', SweepTime);
fprintf('Acquisition time  : %.3f s\n', T);
fprintf('Centre frequency  : %.3f GHz\n', CentreFreq / 1e9);
fprintf('Bandwidth         : %.3f MHz\n', Bandwidth / 1e6);
fprintf('Sample rate       : %.3f kHz\n', Fs / 1e3);

figure('Name', 'Raw And Detrended Sweep', 'Color', 'w');
subplot(2, 1, 1);
plot(fastTime * 1e3, rawData(:, 1), 'LineWidth', 1.0);
grid on;
xlabel('Fast time (ms)');
ylabel('Amplitude');
title('First Raw Sweep');

subplot(2, 1, 2);
plot(fastTime * 1e3, procData(:, 1), 'LineWidth', 1.0);
grid on;
xlabel('Fast time (ms)');
ylabel('Amplitude');
title('First Sweep After DC / Trend Removal');

if Bandwidth > 0
    % FMCW processing: FFT along fast time gives the range profile.
    sweepSlope = Bandwidth / SweepTime;
    NfftRange = 4 * 2^nextpow2(Nsamp);
    rangeWindow = hann(Nsamp, 'periodic');
    rangeFFT = fft(procData .* rangeWindow, NfftRange, 1);
    rangeFFT = rangeFFT(1:NfftRange/2, :);

    beatFreq = (0:NfftRange/2-1).' * (Fs / NfftRange);
    rangeAxis = c * beatFreq / (2 * sweepSlope);
    rangeMap = abs(rangeFFT);
    rangeMapdB = 20 * log10(rangeMap / max(rangeMap(:)) + eps);

    meanProfile = mean(rangeMap, 2);
    peakMask = false(size(meanProfile));
    peakMask(2:end-1) = meanProfile(2:end-1) >= meanProfile(1:end-2) & ...
        meanProfile(2:end-1) >= meanProfile(3:end);
    peakIdx = find(peakMask);
    if isempty(peakIdx)
        [~, peakIdx] = max(meanProfile);
    end
    [~, order] = sort(meanProfile(peakIdx), 'descend');
    peakIdx = peakIdx(order(1:min(2, numel(order))));
    estRanges = sort(rangeAxis(peakIdx));

    figure('Name', 'FMCW Range Processing', 'Color', 'w');
    subplot(2, 1, 1);
    plot(rangeAxis, 20 * log10(meanProfile / max(meanProfile) + eps), ...
        'LineWidth', 1.2);
    grid on;
    xlabel('Range (m)');
    ylabel('Magnitude (dB)');
    title('Average Range Profile');
    xlim([0, max(rangeAxis)]);
    hold on;
    for k = 1:numel(estRanges)
        xline(estRanges(k), '--r', sprintf('%.2f m', estRanges(k)), ...
            'LabelVerticalAlignment', 'bottom');
    end

    subplot(2, 1, 2);
    imagesc(slowTime, rangeAxis, rangeMapdB);
    axis xy;
    colormap turbo;
    colorbar;
    caxis([-45 0]);
    xlabel('Slow time (s)');
    ylabel('Range (m)');
    title('Range-Time Intensity Map');

    % Doppler from the strongest range bin.
    [~, strongestBin] = max(meanProfile);
    rangeBinSignal = squeeze(rangeFFT(strongestBin, :)).';
    NfftDoppler = 2^nextpow2(NSweeps);
    dopplerSpec = fftshift(fft(rangeBinSignal .* hann(NSweeps).', NfftDoppler));
    fdAxis = (-NfftDoppler/2:NfftDoppler/2-1) * (1 / (T / NSweeps) / NfftDoppler);
    velAxis = fdAxis * lambda / 2;

    figure('Name', 'FMCW Doppler At Strongest Range Bin', 'Color', 'w');
    plot(velAxis, 20 * log10(abs(dopplerSpec) / max(abs(dopplerSpec)) + eps), ...
        'LineWidth', 1.2);
    grid on;
    xlabel('Radial velocity (m/s)');
    ylabel('Magnitude (dB)');
    title(sprintf('Doppler Spectrum At %.2f m', rangeAxis(strongestBin)));

    fprintf('Estimated target range(s): %s m\n', num2str(estRanges.', '%.3f '));
else
    % CW processing: no range information, so analyse the Doppler content.
    cwSignal = procData(:);
    totalTime = (0:numel(cwSignal)-1).' / Fs;
    NfftDoppler = 4 * 2^nextpow2(numel(cwSignal));
    dopplerSpec = fftshift(fft(cwSignal .* hann(numel(cwSignal)), NfftDoppler));
    fdAxis = (-NfftDoppler/2:NfftDoppler/2-1).' * (Fs / NfftDoppler);
    velAxis = fdAxis * lambda / 2;
    dopplerMag = abs(dopplerSpec);
    [~, peakIdx] = max(dopplerMag);
    peakVelocity = velAxis(peakIdx);

    figure('Name', 'CW Doppler Processing', 'Color', 'w');
    subplot(2, 1, 1);
    plot(totalTime, cwSignal, 'LineWidth', 1.0);
    grid on;
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('CW Baseband Signal After DC / Trend Removal');

    subplot(2, 1, 2);
    plot(velAxis, 20 * log10(dopplerMag / max(dopplerMag) + eps), ...
        'LineWidth', 1.2);
    grid on;
    xlabel('Radial velocity (m/s)');
    ylabel('Magnitude (dB)');
    title('CW Doppler Spectrum');

    fprintf('Strongest Doppler velocity estimate: %.4f m/s\n', peakVelocity);
end
