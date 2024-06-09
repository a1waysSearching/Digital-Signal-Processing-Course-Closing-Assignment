%% 评价相似性
function evaluation(a1, a2, fs)
s_in = a1;
filtered_s = a2;

% 时域欧氏距离
% Ensure both signals have the same length
minLength = min(length(s_in), length(filtered_s));
audio1 = s_in(1:minLength);
audio2 = filtered_s(1:minLength);
% Calculate Euclidean distance
euclideanDistance = sqrt(sum((audio1 - audio2).^2));
disp(['Euclidean Distance: ', num2str(euclideanDistance)]); 

% 梅尔频率倒谱系数的欧氏距离
% Compute MFCC
mfcc1 = mfcc(s_in, fs);
mfcc2 = mfcc(filtered_s, fs);
% Calculate Euclidean distance between MFCCs
mfccDistance = sqrt(sum((mfcc1(:) - mfcc2(:)).^2));
disp(['MFCC Distance: ', num2str(mfccDistance)]);

% Mel频谱对比
% Compute Mel spectrograms
melSpectrogram1 = melSpectrogram(s_in, fs);
melSpectrogram2 = melSpectrogram(filtered_s, fs);
% Calculate Euclidean distance between Mel spectrograms
melDistance = sqrt(sum((melSpectrogram1(:) - melSpectrogram2(:)).^2));
disp(['Mel Spectrogram Distance: ', num2str(melDistance)]);

% % 动态时间规整
% % Dynamic Time Warping
% [dist, ~, ~] = dtw(mfcc1', mfcc2');
% disp(['DTW Distance: ', num2str(dist)]);

% % 频谱相关性
% % Compute spectral features (e.g., spectral centroid, spectral flux)
% spectralCentroid1 = spectralCentroid(s_in, fs);
% spectralCentroid2 = spectralCentroid(filtered_s, fs);
% % Compare spectral features
% spectralCentroidDistance = sqrt(sum((spectralCentroid1 - spectralCentroid2).^2));
% disp(['Spectral Centroid Distance: ', num2str(spectralCentroidDistance)]);

% % 互相关
% % Cross-correlation
% [correlation, lags] = xcorr(s_in, filtered_s);
% % Find the maximum correlation value
% [maxCorr, idx] = max(correlation);
% lagAtMaxCorr = lags(idx);
% disp(['Maximum Cross-Correlation: ', num2str(maxCorr), ' at lag ', num2str(lagAtMaxCorr)]);


end
