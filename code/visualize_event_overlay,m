function visualize_event_overlay(video_data, sorted_events, video_path, accum_time)
% VISUALIZE_EVENT_OVERLAY 
% Overlay ON/OFF events on original color video frames, with timing aligned
% to downsampled video from load_x4k_frames.
%
% Inputs:
%   video_data    - struct from load_x4k_frames (fields: t_us, H, W, N)
%   sorted_events - struct with fields t, x, y, p
%   video_path    - path to original color video
%   accum_time    - accumulation window in seconds (e.g., 0.01 = 10 ms)

    % --- Unpack video data ---
    t_us_video = video_data.t_us; % timestamps in µs
    H          = video_data.H;
    W          = video_data.W;
    N          = video_data.N;

    % --- Event timestamps (in µs) ---
    event_time_us = double(sorted_events.t);
    accum_time_us = accum_time * 1e6; % convert seconds → µs

    % --- Open original color video ---
    v = VideoReader(video_path);

    % --- Figure ---
    figure('Position', [100,100,800,600]);

    % --- Loop over frames ---
    for k = 1:N
        if ~hasFrame(v), break; end
        frame_color = readFrame(v);

        % Resize to match video_data if needed
        if size(frame_color,1) ~= H || size(frame_color,2) ~= W
            frame_color = imresize(frame_color, [H W]);
        end

        % Time window for events
        t_start_us = t_us_video(k);
        t_end_us   = min(t_start_us + accum_time_us, t_us_video(end));

        % Select events in this window
        idx = (event_time_us >= t_start_us) & (event_time_us < t_end_us);
        evs.x = double(sorted_events.x(idx));
        evs.y = double(sorted_events.y(idx));
        evs.p = sorted_events.p(idx);

        % --- Show video frame ---
        imshow(frame_color); hold on;

        % --- Overlay ON/OFF events ---
        onIdx  = evs.p == 1;
        offIdx = evs.p == -1;
        plot(evs.x(onIdx),  evs.y(onIdx),  'r.', 'MarkerSize', 5); % ON events
        plot(evs.x(offIdx), evs.y(offIdx), 'b.', 'MarkerSize', 5); % OFF events

        title(sprintf('Frame %d | Events %.2f–%.2f ms', ...
              k, double(t_start_us)/1000, double(t_end_us)/1000));

        hold off;
        drawnow;
        pause(0.01);
    end
end
