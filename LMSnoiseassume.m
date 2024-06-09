function LMSnoiseassume(audio)
% Load the audio signal
% [audio, fs] = audioread('audio.wav');

% Use a simple adaptive filter for noise estimation
M = 42; % Filter order
mu = 5e-6; % Step size

% Initialize variables
n = length(audio);
w = zeros(M, 1);
y = zeros(n, 1);
e = zeros(n, 1);

% Apply LMS algorithm
for i = M:n
    x = audio(i:-1:i-M+1);
    y(i) = w' * x;
    e(i) = audio(i) - y(i);
    w = w + mu * x * e(i);
end

% Estimate noise
noise = audio - e;

% Estimate the power of the noise
noisePower = mean(noise.^2);

% Estimate the power of the signal
signalPower = mean(audio.^2);

% Calculate SNR
snrValue = 10 * log10(signalPower / noisePower);
disp(['Estimated SNR: ', num2str(snrValue), ' dB']);

% % Calculate PSNR
% peakSignal = max(abs(audio));
% psnrValue = 10 * log10(peakSignal^2 / noisePower);
% disp(['Estimated PSNR: ', num2str(psnrValue), ' dB']);
end
