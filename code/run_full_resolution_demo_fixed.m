function [video_data, sorted_events, stats] = run_full_resolution_demo_fixed(seq_dir, fps, resize_to, C, opts)
% DEMO: Pixel Array Construction
% [video_data, sorted_events, stats] = run_full_resolution_demo_fixed(seq_dir, fps, resize_to, C, opts)

    arguments
        seq_dir   (1,:) char
        fps       (1,1) double {mustBePositive}
        resize_to (1,2) double {mustBePositive} = [540 1024]
        C         (1,1) double {mustBePositive} = 0.5
        opts struct = struct()
    end

    fprintf('=== Full Resolution Pixel Array Construction ===\n');

    %% ---- Load Video ----
    id = load_x4k_frames(seq_dir, fps, resize_to);

    Y = id.Y;
    t_us = id.t_us;
    H = id.H;
    W = id.W;
    N = id.N;

    %% ---- Configuration ----
    if ~isfield(opts, 'time_resolution_us'); opts.time_resolution_us = 10; end
    if ~isfield(opts, 'Lref0'); opts.Lref0 = []; end
    if ~isfield(opts, 'refractory_us'); opts.refractory_us = 1; end
    if ~isfield(opts, 'threshold_jitter'); opts.threshold_jitter = 0; end
    if ~isfield(opts, 'threshold_sigma_px'); opts.threshold_sigma_px = 0; end

    if ~isfield(opts, 'pr'); opts.pr = struct(); end
    if ~isfield(opts.pr, 'enable'); opts.pr.enable = true; end
    if ~isfield(opts.pr, 'tau_dark_us'); opts.pr.tau_dark_us = 5000; end
    if ~isfield(opts.pr, 'tau_bright_us'); opts.pr.tau_bright_us = 500; end

    if ~isfield(opts, 'leak'); opts.leak = struct(); end
    if ~isfield(opts.leak, 'enable'); opts.leak.enable = false; end

    if ~isfield(opts, 'timing_jitter'); opts.timing_jitter = struct(); end
    if ~isfield(opts.timing_jitter, 'enable'); opts.timing_jitter.enable = false; end

    fprintf('Configuration: C=%.2f, refractory=%dμs\n', C, opts.refractory_us);

    %% ---- Parallel Processing of Entire Image (Fixed Version) ----
    fprintf('\n=== Starting Parallel Processing ===\n');

    event_cells = cell(H, 1);
    start_time = tic;

    parfor y = 1:H
        row_events = struct('t', [], 'x', [], 'y', [], 'p', []);
        for x = 1:W
            L = squeeze(Y(y, x, :));
            L_log = log(L + 1e-3);
            % [t_events, p_events] = pixel_events_pro(L_log, t_us, C, opts);
            [t_events, p_events] = pixel_events_pro_fast(L_log, t_us, C, opts, L);

            if ~isempty(t_events)
                row_events.t = [row_events.t; double(t_events)];
                row_events.x = [row_events.x; repmat(x, length(t_events), 1)];
                row_events.y = [row_events.y; repmat(y, length(t_events), 1)];
                row_events.p = [row_events.p; p_events];
            end
        end
        event_cells{y} = row_events;
    end

    %% ---- Merge Events from All Rows ----
    fprintf('Merging events from all rows...\n');
    all_events = struct('t', [], 'x', [], 'y', [], 'p', []);
    for y = 1:H
        if ~isempty(event_cells{y}.t)
            all_events.t = [all_events.t; event_cells{y}.t];
            all_events.x = [all_events.x; event_cells{y}.x];
            all_events.y = [all_events.y; event_cells{y}.y];
            all_events.p = [all_events.p; event_cells{y}.p];
        end
    end

    processing_time = toc(start_time);
    fprintf('Processing completed in %.1f seconds\n', processing_time);

    %% ---- Results / Stats (no plotting here) ----
    total_events = length(all_events.t);
    if total_events > 0
        [sorted_t, sort_idx] = sort(all_events.t);
        sorted_events = struct(...
            't', sorted_t, ...
            'x', all_events.x(sort_idx), ...
            'y', all_events.y(sort_idx), ...
            'p', all_events.p(sort_idx) ...
        );
        event_rate = total_events / (double(t_us(end)) * 1e-6);
        event_intervals = diff(sorted_events.t);
        polarity_counts = [sum(sorted_events.p == 1), sum(sorted_events.p == -1)];
    else
        sorted_events = struct('t',[],'x',[],'y',[],'p',[]);
        event_rate = 0;
        event_intervals = [];
        polarity_counts = [0,0];
    end

    stats = struct();
    stats.total_events     = total_events;
    stats.avg_per_pixel    = total_events / (W*H);
    stats.event_rate_hz    = event_rate;
    stats.processing_time  = processing_time;
    stats.H                = H;
    stats.W                = W;
    stats.N                = N;
    stats.polarity_counts  = polarity_counts;
    stats.mean_interval_us = iff(~isempty(event_intervals), mean(event_intervals), NaN);
    stats.min_interval_us  = iff(~isempty(event_intervals), min(event_intervals), NaN);
    stats.max_interval_us  = iff(~isempty(event_intervals), max(event_intervals), NaN);
    stats.t_us             = t_us;
    stats.C                = C;
    stats.opts             = opts;

    video_data = id;

    fprintf('\n=== SUMMARY REPORT ===\n');
    fprintf('Image Resolution: %d x %d pixels\n', W, H);
    fprintf('Total Events Generated: %d\n', total_events);
    fprintf('Average Events per Pixel: %.3f\n', stats.avg_per_pixel);
    fprintf('Overall Event Rate: %.0f events/second\n', event_rate);
    fprintf('Processing Time: %.1f seconds\n', processing_time);
    fprintf('Configuration: C=%.2f, refractory=%dμs\n', C, opts.refractory_us);
    if total_events > 0
        fprintf('ON: %d, OFF: %d\n', polarity_counts(1), polarity_counts(2));
    end
    fprintf('=== END ===\n');
end

function out = iff(cond, a, b)
    if cond, out = a; else, out = b; end
end


% % ===== Example (optional, for standalone quick test) =====
% if ~isdeployed && ~ismcc && ~isfolder('skip_this_block')
%     % quick demo (remove or guard as needed)
%     seq_dir   = '../test_video/Type1/TEST01_003_f0433_2k.mp4';
%     fps       = 1000;
%     resize_to = [540 1024];
%     C         = 0.5;
%     opts      = struct();
%     [vid, evs, st] = run_full_resolution_demo(seq_dir, fps, resize_to, C, opts);
%     % Now you could call:
%     visualize_event_overlay_formalapp(ax, vid, evs, seq_dir, 0.01);
% end