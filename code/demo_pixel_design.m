% DEMO: Pixel Design using X4K1000FPS clip
clear; clc;

%% ---- load video ----
seq_dir   = '../test_video/Type1/TEST01_003_f0433_2k.mp4';   % path to your test video
fps       = 1000;               % true frame rate of X4K1000FPS
resize_to = [540 1024];         % optional downscale from 2K (2048x1080) for faster processing
id = load_x4k_frames(seq_dir, fps, resize_to);

Y = id.Y;
t_us = id.t_us;
H = id.H;
W = id.W;
N = id.N;

% ---- choose a pixel (e.g., center or a high-gradient point) ----
y = round(H/2);
x = round(W/2);
% Alternative: pick a pixel with stronger motion/edges
% G = imgradient(Y(:,:,1)); [~,lin] = max(G(:)); [y,x] = ind2sub(size(G),lin);

% ---- build log-intensity sequence for this pixel ----
eps_ = 1e-3;               % to avoid log(0)
L = squeeze(Y(y,x,:));
L_log = log( L + eps_ );


%% 1. basic simulator
if 0  %(deactivated)
C = 0.05;              % contrast threshold (log-intensity)
opts = struct;

opts.Lref0            = [];     % default: first sample
opts.time_resolution_us = 1;    % timestamp resolution
opts.refractory_us    = 0;      % try 2–5 for deglitching if needed
opts.threshold_jitter = 0;      % try 0.01 for realism

% Generate events for this pixel
[t_e, p_e] = pixel_events(L_log, t_us, C, opts);
fprintf('Pixel (%d,%d): %d events generated.\n', y,x, numel(t_e));
end


%% advanced simulator
if 1 %(activated)
C = 0.05;              % contrast threshold (log-intensity)
opts = struct;

% Baseline
opts.time_resolution_us = 1; % timestamp resolution
opts.refractory_us = 0;      % disable refractory period
opts.Lref0         = [];     % default = first frame

% Threshold noise
opts.threshold_sigma_px = 0.02;   % per-pixel Gaussian sample (0 = disable)
opts.threshold_jitter   = 0.00;   % per-event jitter (0 = disable)

% Photoreceptor bandwidth (optional)
opts.pr.enable        = true;
opts.pr.tau_dark_us   = 500;     % larger time constant in dark regions
opts.pr.tau_bright_us = 50;      % smaller time constant in bright regions

% Leak events (optional)
opts.leak.enable          = true;
opts.leak.rate_dark_hz    = 1;    % ~1 leak event/sec in dark
opts.leak.rate_bright_hz  = 10;   % ~10 leak events/sec in bright
opts.leak.on_prob_base    = 0.5;  % ON/OFF probability balance
opts.leak.on_prob_slope   = 0.0;  % not dependent on brightness
opts.leak.update_Lref     = false;% recommended: false (do not disturb reference)

% Timing jitter (optional)
opts.timing_jitter.enable    = true;
opts.timing_jitter.dark_us   = 50; % larger noise in dark regions
opts.timing_jitter.bright_us = 5;  % smaller noise in bright regions


% Run advanced simulator
start_time = tic;

% [t_e, p_e] = pixel_events_pro(L_log, t_us, C, opts);         % first version
[t_e, p_e] = pixel_events_pro_fast(L_log, t_us, C, opts, L); % second version

processing_time = toc(start_time);
fprintf('Processing completed in %.5f seconds\n', processing_time);
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