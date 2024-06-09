function BINnoiseassume(audio)

% 计算音频信号的直方图
num_bins = 110; % 直方图的bin数
[counts, edges] = histcounts(audio, num_bins);

% 估计噪声水平
% 我们假设噪声水平在直方图中最频繁出现的值附近
[~, max_count_idx] = max(counts);
noise_level = (edges(max_count_idx) + edges(max_count_idx + 1)) / 2;

% 计算信噪比
signal_power = mean(audio.^2);
noise_power = noise_level^2;
snr = 10 * log10(signal_power / noise_power);

% 输出信噪比
fprintf('Estimated SNR: %.2f dB\n', snr);
end