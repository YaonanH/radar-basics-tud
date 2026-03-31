clc;
clear;
close all;

% Minimal Lab 3 script:
% 1. Load one FMCW capture
% 2. Perform range compression with different windows
% 3. Keep only the left half of the real-valued spectrum
% 4. Convert beat-frequency bins into range
% 5. Compare the range profiles in dB

[BScans, SweepTime, ~, Bandwidth, ~] = ReadSiversBIN;

if Bandwidth <= 0
    error('Please select an FMCW file with non-zero bandwidth.');
end

c = physconst('LightSpeed');
data = detrend(double(BScans), 'linear');
[Nsamp, ~] = size(data);
Fs = Nsamp / SweepTime;
slope = Bandwidth / SweepTime;

Nfft = Nsamp;

beatFreqAxis = (0:Nfft/2-1).' * (Fs / Nfft);
rangeAxis = beatFreqAxis * c / (2 * slope);

windowNames = {'Rectangular', 'Hamming', 'Blackman'};
windowList = {
    ones(Nsamp, 1), ...
    hamming(Nsamp, 'periodic'), ...
    blackman(Nsamp, 'periodic')
    };

figure('Name', 'Range Compression Comparison', 'Color', 'w');
hold on;
grid on;

for k = 1:numel(windowList)
    win = windowList{k};
    spectrum = fft(data .* win, Nfft, 1);
    spectrum = spectrum(1:Nfft/2, :); % Use only left half for real data.

    profile = mean(abs(spectrum), 2);
    profiledB = 20 * log10(profile / max(profile) + eps);

    plot(rangeAxis, profiledB, 'LineWidth', 1.2, ...
        'DisplayName', windowNames{k});
end

xlabel('Range (m)');
ylabel('Magnitude (dB)');
title('Range Profiles For Different Window Functions');
legend('Location', 'best');
xlim([0, max(rangeAxis)]);
ylim([-60, 5]);
