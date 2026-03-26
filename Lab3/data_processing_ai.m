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
% For sweep-by-sweep processing, the slow-time sampling interval is the
% sweep repetition time.
slowTime = (0:NSweeps-1) * SweepTime;

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

    % Track the pendulum by looking for the strongest reflector inside a
    % physically plausible range gate around the dominant target.
    validRange = rangeAxis >= 0.3 & rangeAxis <= 5.0;
    rangeAxisTrack = rangeAxis(validRange);
    rangeMapTrackAll = rangeMap(validRange, :);
    meanProfileTrack = mean(rangeMapTrackAll, 2);
    [~, dominantTrackIdx] = max(meanProfileTrack);
    dominantRange = rangeAxisTrack(dominantTrackIdx);
    trackHalfWidth = 1.2;
    validTrackGate = rangeAxisTrack >= max(0.3, dominantRange - trackHalfWidth) & ...
        rangeAxisTrack <= min(5.0, dominantRange + trackHalfWidth);
    rangeAxisGate = rangeAxisTrack(validTrackGate);
    rangeMapGate = rangeMapTrackAll(validTrackGate, :);

    [peakStrength, peakIdxTrack] = max(rangeMapGate, [], 1);
    pendulumRangeTrack = rangeAxisGate(peakIdxTrack);
    noiseFloor = median(rangeMapGate, 1);
    globalThreshold = 0.03 * max(rangeMapGate(:));
    detectionThreshold = max(2.0 * noiseFloor, globalThreshold);
    weakDetections = peakStrength < detectionThreshold;
    pendulumRangeTrack(weakDetections) = NaN;
    pendulumRangeTrack = fillmissing(pendulumRangeTrack, 'previous');
    pendulumRangeTrack = fillmissing(pendulumRangeTrack, 'nearest');

    figure('Name', 'FMCW Range Migration', 'Color', 'w');
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
    caxis([-20 0]);
    xlabel('Slow time (s)');
    ylabel('Range (m)');
    title('Normalised Reflectivity (dB)');
    ylim([0 5]);

    figure('Name', 'Pendulum Motion From FMCW', 'Color', 'w');
    plot(slowTime, pendulumRangeTrack, 'LineWidth', 1.2);
    grid on;
    xlabel('Slow time (s)');
    ylabel('Range (m)');
    title('Pendulum Motion');

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
    % CW processing: FFT each sweep into Doppler, then track it in slow time.
    NfftDoppler = 4 * 2^nextpow2(Nsamp);
    dopplerWindow = hann(Nsamp, 'periodic');
    dopplerFFT = fft(procData .* dopplerWindow, NfftDoppler, 1);
    dopplerFFT = dopplerFFT(1:NfftDoppler/2, :);

    fdAxis = (0:NfftDoppler/2-1).' * (Fs / NfftDoppler);
    velAxis = fdAxis * lambda / 2;
    dopplerMap = abs(dopplerFFT);
    dopplerMapdB = 20 * log10(dopplerMap / max(dopplerMap(:)) + eps);

    % Ignore the DC bin and unrealistic high-velocity bins when tracking
    % the pendulum peak; otherwise isolated noise spikes dominate.
    vmaxTrack = min(7, velAxis(end));
    validBins = velAxis >= 0.15 & velAxis <= vmaxTrack;
    dopplerMapTrack = dopplerMap(validBins, :);
    velAxisTrack = velAxis(validBins);
    [peakStrength, peakIdx] = max(dopplerMapTrack, [], 1);
    peakVelocityTrack = velAxisTrack(peakIdx);

    % Keep the original per-sweep peak search, but reject detections whose
    % spectral magnitude is too weak compared with the background.
    noiseFloor = median(dopplerMapTrack, 1);
    globalThreshold = 0.03 * max(dopplerMapTrack(:));
    detectionThreshold = max(1.5 * noiseFloor, globalThreshold);
    peakVelocityTrack(peakStrength < detectionThreshold) = 0;

    % Oscillation frequency from the velocity envelope in slow time.
    peakVelocityZeroMean = peakVelocityTrack - mean(peakVelocityTrack);
    NfftOsc = 2^nextpow2(NSweeps);
    oscSpec = fft(peakVelocityZeroMean .* hann(NSweeps).', NfftOsc);
    oscFreqAxis = (0:NfftOsc/2-1) * ((1 / SweepTime) / NfftOsc);
    oscMag = abs(oscSpec(1:NfftOsc/2));

    figure('Name', 'CW Doppler Processing', 'Color', 'w');
    subplot(3, 1, 1);
    plot(fastTime * 1e3, procData(:, 1), 'LineWidth', 1.0);
    grid on;
    xlabel('Fast time (ms)');
    ylabel('Amplitude');
    title('First CW Sweep After DC / Trend Removal');

    subplot(3, 1, 2);
    imagesc(slowTime, velAxis, dopplerMapdB);
    axis xy;
    colormap turbo;
    colorbar;
    caxis([-45 0]);
    ylim([0 vmaxTrack]);
    xlabel('Slow time (s)');
    ylabel('V_r (m/s)');
    title('Pendulum Doppler Spectrum');

    subplot(3, 1, 3);
    plot(slowTime, peakVelocityTrack, 'LineWidth', 1.2);
    grid on;
    xlabel('Slow time (s)');
    ylabel('V_r (m/s)');
    title('Maximum Of Doppler Spectrum');

    figure('Name', 'Pendulum Oscillation Frequency', 'Color', 'w');
    plot(oscFreqAxis, 20 * log10(oscMag / max(oscMag) + eps), ...
        'LineWidth', 1.2);
    grid on;
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    title('FFT Of Peak Doppler Velocity Track');

    [~, peakIdxOsc] = max(oscMag(2:end));
    peakIdxOsc = peakIdxOsc + 1;
    pendulumFrequency = oscFreqAxis(peakIdxOsc);

    fprintf('Maximum observed Doppler velocity: %.4f m/s\n', ...
        max(peakVelocityTrack));
    fprintf('Estimated pendulum oscillation frequency: %.4f Hz\n', ...
        pendulumFrequency);
end
