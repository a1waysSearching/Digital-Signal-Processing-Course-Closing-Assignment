% 清除工作区变量
clear
% 清除命令行窗口
clc

% 各种音效的参数设置
% darth vader: 0.6, simple average
% chipmunk: 2.5, single
% robot: 1.6, simple average
% changing: 1.6, bufgain + modulate
% chorus: 1.2, +s_in 

% 环形缓冲区大小
B = 1024;

% 测试音高转换，读取音频文件
[s_in, fs] = audioread('去噪测试音频样本1.m4a');
s_in = s_in +1e-20; % 防止计算出现NAN
% 采样频率
% fs=44100;
% x = linspace(0,1,fs);
% s_in = sin(2*pi*440*x);

% 定义环形缓冲区
ringbuffer1 = zeros(B, 1);
ringbuffer2 = zeros(B, 1);
ringbuffer3 = zeros(B, 1);
ringbuffer4 = zeros(B, 1);

% 初始化输出信号
s_avg = zeros(length(s_in), 1);
s_no_cross = zeros(length(s_in), 1);
s_chord1 = zeros(length(s_in), 1);
s_chord2 = zeros(length(s_in), 1);

% 定义音高转换的增量
delta1 = 0.25^(-3/12); % 改变音高，负幂为降低音高，正幂为升高音高
delta2 = 0.25^(-6/12);

% RB的索引地址
idx1 = 0;
idx2 = 0;
idx3 = 0;
idx4 = 0;

for tt = 1:length(s_in)
    % 写入 环形缓冲区1
    w_addr1 = mod(tt, B) + 1;
    ringbuffer1(w_addr1) = s_in(tt);
    
    % 写入 环形缓冲区2
    w_addr2 = mod(tt + round(B / 2), B) + 1;
    ringbuffer2(w_addr2) = s_in(tt);   
    
    % 从环形缓冲区1 读取的地址
    idx1 = idx1 + delta1;
    r_addr1 = mod(round(idx1) - 1, B) + 1;
   
    % 从环形缓冲区2 读取的地址
    idx2 = idx2 + delta1;
    r_addr2 = mod(round(idx2) - 1, B) + 1;  
    
    % 写入 环形缓冲区3
    w_addr3 = mod(tt, B) + 1;
    ringbuffer3(w_addr3) = s_in(tt);
    
    % 写入 环形缓冲区4
    w_addr4 = mod(tt + round(B / 2), B) + 1;
    ringbuffer4(w_addr4) = s_in(tt);   
    
    % 从环形缓冲区3 读取的地址
    idx3 = idx3 + delta2;
    r_addr3 = mod(round(idx3) - 1, B) + 1;
   
    % 从环形缓冲区4 读取的地址
    idx4 = idx4 + delta2;
    r_addr4 = mod(round(idx4) - 1, B) + 1;  
    
    % 调制部分（未使用）
    % delta = 1 + 0.5 * sin(1 * pi * tt / fs); 
    
    % 平滑
    % 平均交叉淡化
    s_avg(tt, 1) = 0.5 * (ringbuffer1(r_addr1, 1) + ringbuffer2(r_addr2, 1));

    s_chord1(tt, 1) = 0.25 * (ringbuffer1(r_addr1, 1) + ringbuffer2(r_addr2, 1));
    s_chord2(tt, 1) = 0.25 * (ringbuffer3(r_addr3, 1) + ringbuffer4(r_addr4, 1));
    
    % 无交叉淡化
    s_no_cross(tt, 1) = ringbuffer1(r_addr1, 1);
end
out = 0.5 * s_in(:, 1) + s_chord1 + s_chord2;
% 全权重叠加
% soundsc(0.5 * s_in, fs)
% pause(0.005)
% soundsc(s_chord1, fs)
% pause(0.005)
% soundsc(s_chord2, fs)
% 
% pause 
%% 巴特沃斯低通滤波器去除3500Hz以上高频噪声
Hd = bbutter; 
filtered_s = filter(Hd, out);

%% 播放最终的音频
% soundsc(0.5 * s_in(:, 1) + s_chord1 + s_chord2, fs)
% soundsc(filtered_s, fs)

%% 评价
fprintf('----相似性计算----\n')
evaluation(s_in(:, 1), filtered_s, fs);

fprintf('----LMS计算信噪比----\n')
fprintf('单个缓冲区结果的信噪比 ')
LMSnoiseassume(s_chord1);
fprintf('双缓冲区结果的信噪比 ')
LMSnoiseassume(out);
fprintf('双缓冲区加低通滤波器结果的信噪比 ')
LMSnoiseassume(filtered_s);

fprintf('----直方图计算信噪比----\n')
fprintf('单个缓冲区结果的信噪比 ')
BINnoiseassume(s_chord1);
fprintf('双缓冲区结果的信噪比 ')
BINnoiseassume(out);
fprintf('双缓冲区加低通滤波器结果的信噪比 ')
BINnoiseassume(filtered_s);

%% 绘图
b_sound_time = s_in; % 滤波前时域表示
b_sound_fft = fft(b_sound_time); % 傅里叶变换
b_f = (0:length(b_sound_fft)-1)*fs/length(b_sound_fft); % 频率轴
signal = ifft(b_sound_fft); % IFFT操作 

a_sound_time = filtered_s; % 滤波后时域表示
a_sound_fft = fft(a_sound_time); % 傅里叶变换
a_f = (0:length(a_sound_fft)-1)*fs/length(a_sound_fft); % 频率轴
reconstructed_signal = ifft(a_sound_fft); % IFFT操作 

subplot(231);plot(s_in);xlabel('时间');ylabel('幅度');title('输入男声波形');
subplot(234);plot(filtered_s);xlabel('时间');ylabel('幅度');title('输出女声波形');
subplot(232);plot(b_f, abs(b_sound_fft));xlabel('频率 (Hz)');ylabel('频谱');title('输入男声信号的频谱');
subplot(235);plot(a_f, abs(a_sound_fft));xlabel('频率 (Hz)');ylabel('频谱');title('输出女声的频谱');
subplot(233);myspectrogram(s_in(:, 1),fs);colormap(jet);time=(0:length(s_in(:, 1))-1)/fs;axis([0 max(time*1000) 0 8000]);title('输入男声的时谱图');xlabel('时间');ylabel('频率');
subplot(236);myspectrogram(filtered_s(:, 1),fs);colormap(jet);time=(0:length(filtered_s(:, 1))-1)/fs;axis([0 max(time*1000) 0 8000]);title('输出女声的时谱图');xlabel('时间');ylabel('频率');


% figure
% plot(s_in, 'r')
% hold on
% plot(s_avg, 'b')
% set(gca, 'FontSize', 16);
% title('avg cross fading')
% legend('input sine wave', 'shifted', 'Location', 'Best')
% 
% s_med = medfilt1(s_avg, 5);
% 
% figure
% plot(s_avg, 'r')
% hold on
% plot(s_med, 'b')
% set(gca, 'FontSize', 16);
% title('Adding a median filter')
% legend('unfiltered', 'median filtered', 'Location', 'Best')
% 
% fc = 3500;
% [b,a] = butter(6, fc/(fs/2));
% s_med_butter = filter(b, a, s_med);
% s_butter = filter(b, a, s_avg);
% 
% figure
% plot(s_med, 'r')
% hold on
% plot(s_med_butter, 'b')
% hold on
% plot(s_butter, 'g')
% set(gca, 'FontSize', 16);
% title('Testing different filters')
% legend('median', 'median and butter', 'butter', 'Location', 'Best')

