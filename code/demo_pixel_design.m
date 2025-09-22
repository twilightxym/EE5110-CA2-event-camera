% DEMO: Pixel Design using X4K1000FPS clip
clear; clc;


%% ---- config ----
seq_dir   = '../test_video/Type1/TEST01_003_f0433_2k.mp4';   % set to your frames folder
fps       = 1000;               % X4K1000FPS frame rate

resize_to = [540 1024];          % downscale from 2K(2048*1080) to speed up (optional)
C         = 0.005;               % contrast threshold (log-intensity)
eps_      = 1e-3;               % avoid log(0)


% ---- load frames ----
[Y, t_us] = load_x4k_frames(seq_dir, fps, resize_to);
[H,W,N]   = size(Y);
fprintf('Loaded %dx%dx%d frames @ %d fps\n', H,W,N,fps);

% ---- choose a pixel (e.g., center or a high-gradient point) ----
y = round(H/4);
x = round(W/4);
% You can also pick a pixel with stronger motion/edges:
% G = imgradient(Y(:,:,1)); [~,lin] = max(G(:)); [y,x] = ind2sub(size(G),lin);

% ---- build log-intensity sequence for this pixel ----
L_log = log( squeeze(Y(y,x,:)) + eps_ );

%% advanced simulator
if 1
opts = struct;

% 基线
opts.refractory_us = 0;     % 先关不应期
opts.Lref0         = [];    % 默认用首帧

% 阈值噪声
opts.threshold_sigma_px = 0.02;   % 像素级一次抽样（关掉则=0）
opts.threshold_jitter   = 0.00;   % 逐事件抖动（关掉则=0）

% 光感受器带宽（可选）
opts.pr.enable        = true;     % 打开
opts.pr.tau_dark_us   = 5000;     % 暗处带宽小（时间常数大）
opts.pr.tau_bright_us = 500;      % 亮处带宽大（时间常数小）

% Leak events（可选）
opts.leak.enable          = true;
opts.leak.rate_dark_hz    = 1;    % 暗处每秒 ~1 个
opts.leak.rate_bright_hz  = 10;   % 亮处每秒 ~10 个
opts.leak.on_prob_base    = 0.5;  % ON/OFF 50/50
opts.leak.on_prob_slope   = 0.0;  % 不随亮度变
opts.leak.update_Lref     = false;% 建议 false（不扰动参考面）

% 时间抖动（可选）
opts.timing_jitter.enable   = true;
opts.timing_jitter.dark_us  = 50; % 暗处时间噪声更大
opts.timing_jitter.bright_us= 5;  % 亮处更小

% 直接调用：
[t_e, p_e] = pixel_events_pro(L_log, t_us, C, opts);
end


%% basic simulator
if 0
opts.Lref0            = [];     % default: first sample
opts.refractory_us    = 0;      % try 2~5 for deglitching if needed
opts.threshold_jitter = 0;      % try 0.01 for realism

% ---- generate events for this pixel ----
[t_e, p_e] = pixel_events(L_log, t_us, C, opts);
fprintf('Pixel (%d,%d): %d events generated.\n', y,x, numel(t_e));
end

%% ---- quick visualization: intensity vs events ----
figure; 
subplot(2,1,1);
plot(double(t_us)*1e-6, L_log, '-o'); grid on;
xlabel('time (s)'); ylabel('log-intensity');
title(sprintf('Pixel (%d,%d) log-intensity', y,x));

subplot(2,1,2);
stem(double(t_e)*1e-6, p_e, 'filled');
xlabel('time (s)'); ylabel('polarity');
title('Events (+1/-1)'); grid on;

%% Save pixel data for full frame processing
pixel_data.Y = Y;
pixel_data.t_us = t_us;
pixel_data.H = H;
pixel_data.W = W;
pixel_data.N = N;
pixel_data.C = C;
pixel_data.opts = opts;
pixel_data.selected_pixel = [x, y];
pixel_data.selected_events = [t_e, p_e];

save('pixel_analysis_data.mat', 'pixel_data');
fprintf('Pixel data saved for full frame processing.\n');