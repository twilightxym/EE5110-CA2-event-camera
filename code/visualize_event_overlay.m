function visualize_event_overlay(video_data, sorted_events, video_path, accum_time)
% VISUALIZE_EVENT_OVERLAY 
% Overlay ON/OFF events on original color video frames, with timing aligned
% to downsampled video from load_x2k_frames.
% Adaptive accumulation logic: 
% - accumtime < framedt: Use [tk, tk+accumtime) 
% - accumtime >= framedt: Use [tk, t_{k+1}] to avoid future events ahead of time
%
% Inputs:
%   video_data    - struct from load_x2k_frames (fields: t_us, H, W, N)
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
    ex = double(sorted_events.x(:));
    ey = double(sorted_events.y(:));
    ep = sorted_events.p(:);
    M  = numel(event_time_us);

    accum_time_us = accum_time * 1e6; % convert seconds → µs
    
    % --- Per-frame dt ---
    if N >= 2
        frame_dt_us_vec = double(diff(t_us_video));
    else
        frame_dt_us_vec = [];
    end

    % --- Open original color video ---
    v = VideoReader(video_path);

    % --- Figure ---
    fig = figure('Position', [100,100,900,650]);
    ax  = axes('Parent',fig);
    imH = []; hOn = []; hOff = [];

    i0 = 1; i1 = 1;
    for k = 1:N
        if ~hasFrame(v), break; end
        frame_color = readFrame(v);
        if size(frame_color,1) ~= H || size(frame_color,2) ~= W
            frame_color = imresize(frame_color,[H W]);
        end

        t_start_frame = double(t_us_video(k));

        if k < N
            t_next_frame  = double(t_us_video(k+1));
            frame_dt_us   = frame_dt_us_vec(k);
            step_us = accum_time_us;

            % accum_time_us >= framedt, then [t_k, t_{k+1})
            if step_us >= frame_dt_us
                num_sub = 1;
                sub_edges = [t_start_frame, t_next_frame];
            else
                num_sub = ceil(frame_dt_us / step_us);
                sub_edges = zeros(num_sub+1,1);
                for s = 0:num_sub
                    sub_edges(s+1) = min(t_start_frame + s*step_us, t_next_frame);
                end
            end
        else
            t_next_frame = min(t_start_frame + accum_time_us, double(t_us_video(end)));
            num_sub = 1;
            sub_edges = [t_start_frame, t_next_frame];
        end

        if isempty(imH)
            imH = imshow(frame_color,'Parent',ax); hold(ax,'on');
            hOn  = plot(ax, NaN, NaN, 'r.', 'MarkerSize', 1);
            hOff = plot(ax, NaN, NaN, 'b.', 'MarkerSize', 1);
            axis(ax,'image');
        else
            set(imH, 'CData', frame_color);
        end

        for s = 1:num_sub
            t0 = sub_edges(s);
            t1 = sub_edges(s+1);

            while i0 <= M && event_time_us(i0) < t0
                i0 = i0 + 1;
            end

            if i1 < i0, i1 = i0; end

            while i1 <= M && event_time_us(i1) < t1
                i1 = i1 + 1;
            end

            if i0 < i1
                idx = i0:(i1-1);
                xw = ex(idx); yw = ey(idx); pw = ep(idx);
            else
                xw = []; yw = []; pw = [];
            end

            set(hOn,  'XData', xw(pw== 1), 'YData', yw(pw== 1));
            set(hOff, 'XData', xw(pw==-1), 'YData', yw(pw==-1));

            set(ax, 'Title', text(0.5, 1.02, ...
                sprintf('Frame %d  |  Sub %d/%d  |  %.3f–%.3f ms  |  events=%d', ...
                        k, s, num_sub, t0/1000, t1/1000, numel(xw)), ...
                'Units','normalized','HorizontalAlignment','center'));

            drawnow;
            % pause(1);
        end

        hold(ax,'off'); hold(ax,'on'); 
    end
end
