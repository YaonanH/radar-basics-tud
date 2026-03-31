clc;
clear;
close all;

% Minimal CW Doppler processing for Lab 3.
% 1. Load one CW capture
% 2. Remove DC and linear trend from each sweep
% 3. FFT along the first dimension
% 4. Keep only the left half of the spectrum because the data are real
% 5. Convert Doppler frequency bins into absolute radial velocity
% 6. Plot Doppler spectrum versus slow time and the maximum-velocity track

[BScans, SweepTime, CentreFreq, Bandwidth, T] = ReadSiversBIN;

if Bandwidth ~= 0
    error('Please select a CW file with zero bandwidth.');
end

c = physconst('LightSpeed');
lambda = c / CentreFreq;

data = detrend(double(BScans), 'linear');
[Nsamp, NSweeps] = size(data);

Fs = Nsamp / SweepTime;
slowTimeStep = T / NSweeps;
slowTime = (0:NSweeps-1) * slowTimeStep;
Nfft = Nsamp;

dopplerFFT = fft(data, Nfft, 1);
dopplerFFT = dopplerFFT(1:Nfft/2, :); % Use only left half for real data.

dopplerFreqAxis = (0:Nfft/2-1).' * (Fs / Nfft);
velocityAxis = dopplerFreqAxis * lambda / 2;

dopplerMap = abs(dopplerFFT);
dopplerMapdB = 20 * log10(dopplerMap / max(dopplerMap(:)) + eps);

% Track the strongest Doppler component in a realistic velocity region.
validVelocity = velocityAxis >= 0.1 & velocityAxis <= 10;
dopplerMapTrack = dopplerMap(validVelocity, :);
velocityAxisTrack = velocityAxis(validVelocity);

[peakStrength, peakIdx] = max(dopplerMapTrack, [], 1);
peakVelocityTrack = velocityAxisTrack(peakIdx);

% Reject weak detections so the line plot is less sensitive to noise.
noiseFloor = median(dopplerMapTrack, 1);
detectionThreshold = max(2 * noiseFloor, 0.03 * max(dopplerMapTrack(:)));
peakVelocityTrack(peakStrength < detectionThreshold) = 0;

% Optional: estimate pendulum oscillation frequency from the velocity track.
peakVelocityTrack = peakVelocityTrack(:);
trackZeroMean = peakVelocityTrack - mean(peakVelocityTrack);
NfftOsc = 2^nextpow2(NSweeps);
oscWindow = hann(NSweeps);
oscSignal = trackZeroMean .* oscWindow;
oscSpec = fft(oscSignal, NfftOsc, 1);
oscFreqAxis = (0:NfftOsc/2-1).' * ((1 / slowTimeStep) / NfftOsc);
oscMag = abs(oscSpec(1:NfftOsc/2));
oscMag(1) = 0;
maxOscMag = max(oscMag);
oscMagdB = 20 * log10(oscMag / maxOscMag + eps);

figure('Name', 'CW Doppler Processing', 'Color', 'w');
subplot(3, 1, 1);
plot(slowTime, peakVelocityTrack, 'LineWidth', 1.2);
grid on;
xlabel('Slow time (s)');
ylabel('V_r (m/s)');
title('Maximum Of Doppler Spectrum');

subplot(3, 1, 2);
imagesc(slowTime, velocityAxis, dopplerMapdB);
axis xy;
colormap turbo;
colorbar;
clim([-40 0]);
ylim([0 5]);
xlabel('Slow time (s)');
ylabel('V_r (m/s)');
title('Pendulum Doppler Spectrum');

subplot(3, 1, 3);
plot(oscFreqAxis, oscMagdB, 'LineWidth', 1.2);
grid on;
ylim([-40, 0]);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('FFT Of Peak Doppler Velocity Track');
