clc;
clear;
close all;

% Minimal FMCW range migration processing for Lab 3.
% 1. Load one FMCW capture
% 2. Remove DC and linear trend from each sweep
% 3. FFT along the fast-time dimension
% 4. Keep only the left half of the spectrum because the data are real
% 5. Convert beat-frequency bins into range
% 6. Plot range versus slow time and track the maximum reflectivity

[BScans, SweepTime, ~, Bandwidth, T] = ReadSiversBIN;

if Bandwidth <= 0
    error('Please select an FMCW file with non-zero bandwidth.');
end

c = physconst('LightSpeed');
data = detrend(double(BScans), 'linear');
[Nsamp, NSweeps] = size(data);

Fs = Nsamp / SweepTime;
slope = Bandwidth / SweepTime;
slowTimeStep = T / NSweeps;
slowTime = (0:NSweeps-1) * slowTimeStep;
Nfft = Nsamp;

rangeFFT = fft(data, Nfft, 1);
rangeFFT = rangeFFT(1:Nfft/2, :); % Use only left half for real data.

beatFreqAxis = (0:Nfft/2-1).' * (Fs / Nfft);
rangeAxis = beatFreqAxis * c / (2 * slope);

rangeMap = abs(rangeFFT);
rangeMapdB = 20 * log10(rangeMap / max(rangeMap(:)) + eps);

% Track the strongest reflector inside a reasonable pendulum range interval.
validRange = rangeAxis >= 0.3 & rangeAxis <= 5.0;
rangeMapTrack = rangeMap(validRange, :);
rangeAxisTrack = rangeAxis(validRange);

[peakStrength, peakIdx] = max(rangeMapTrack, [], 1);
pendulumRangeTrack = rangeAxisTrack(peakIdx);

% Reject weak detections to reduce random jumps in the range track.
noiseFloor = median(rangeMapTrack, 1);
detectionThreshold = max(2 * noiseFloor, 0.03 * max(rangeMapTrack(:)));
pendulumRangeTrack(peakStrength < detectionThreshold) = NaN;
pendulumRangeTrack = fillmissing(pendulumRangeTrack, 'previous');
pendulumRangeTrack = fillmissing(pendulumRangeTrack, 'nearest');

figure('Name', 'FMCW Range Migration', 'Color', 'w');
subplot(2, 1, 1);
imagesc(slowTime, rangeAxis, rangeMapdB);
axis xy;
colormap turbo;
colorbar;
caxis([-20 0]);
ylim([0 5]);
xlabel('Slow time (s)');
ylabel('Range (m)');
title('Normalised Reflectivity (dB)');

subplot(2, 1, 2);
plot(slowTime, pendulumRangeTrack, 'LineWidth', 1.2);
grid on;
xlabel('Slow time (s)');
ylabel('Range (m)');
title('Pendulum Motion');
