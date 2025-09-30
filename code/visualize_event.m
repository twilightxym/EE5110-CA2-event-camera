function visualize_event(video_data, sorted_events, accum_time, video_path, show_edges)
% VISUALIZE_EVENT
% Left: original color video frame (from video_path)
% Right: event frame synced to video_data timestamps
% Optional third panel: event edges
%
% Inputs:
%   video_data    - struct from load_x4k_frames (grayscale frames, timestamps)
%   sorted_events - struct with fields t, x, y, p
%   accum_time    - accumulation window in seconds (e.g., 0.01 = 10 ms)
%   video_path    - path to original color video
%   show_edges    - (optional) true/false, if true shows 3 panels
%
% Example:
%   visualize_event(id, sorted_events, 0.01, 'myvideo.mp4');
%   visualize_event(id, sorted_events, 0.01, 'myvideo.mp4', true);

    if nargin < 5
        show_edges = false; % default = only 2 panels
    end

    % --- Unpack video data ---
    t_us_video = video_data.t_us; % timestamps in microseconds
    H          = video_data.H;
    W          = video_data.W;
    N          = video_data.N;

    % --- Event timestamps (keep in microseconds) ---
    event_time_us = double(sorted_events.t);

    % --- Accumulation window in microseconds ---
    accum_time_us = accum_time * 1e6;

    % --- Open original color video ---
    v = VideoReader(video_path);

    % --- Figure setup ---
    if show_edges
        fig = figure('Position', [100,100,1600,600]); % 3 panels
    else
        fig = figure('Position', [100,100,1200,600]); % 2 panels
    end

    % --- Loop over frames ---
    for k = 1:N
        if ~hasFrame(v), break; end
        frame_color = readFrame(v); % original RGB frame

        % Time window for events (clamped to video end)
        t_start_us = t_us_video(k);
        t_end_us   = min(t_start_us + accum_time_us, t_us_video(end));

        % Select events within this window
        idx = (event_time_us >= t_start_us) & (event_time_us < t_end_us);
        events_x = double(sorted_events.x(idx));
        events_y = double(sorted_events.y(idx));
        events_p = sorted_events.p(idx);

        % --- Build Event Frame ---
        event_frame = zeros(H, W);
        for i = 1:length(events_x)
            % Clamp coordinates inside frame
            x = min(max(events_x(i)+1, 1), W);
            y = min(max(events_y(i)+1, 1), H);
            if events_p(i) == 1
                event_frame(y, x) = 1;   % ON event
            else
                event_frame(y, x) = -1;  % OFF event
            end
        end

        % --- Visualization ---
        subplot(1, show_edges+2, 1);
        imshow(frame_color); % always show original color frame
        title(sprintf('Video Frame %d (t = %.2f ms)', ...
              k, double(t_start_us)/1000));

        subplot(1, show_edges+2, 2);
        imagesc(event_frame, [-1 1]); % fix scale to -1..1
        colormap(gca, [0 0 1; 1 1 1; 1 0 0]); % Blue=OFF, White=none, Red=ON
        axis image; axis off;
        title(sprintf('Event Frame (%.2f–%.2f ms)', ...
              double(t_start_us)/1000, double(t_end_us)/1000));

        if show_edges
            edge_mask = edge(event_frame ~= 0, 'canny');
            subplot(1,3,3);
            imshow(edge_mask);
            title('Event Edges (Canny)');
        end

        drawnow;
        pause(0.01);
    end
end
