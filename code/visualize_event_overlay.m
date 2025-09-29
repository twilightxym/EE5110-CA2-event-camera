function visualize_event_overlay(videoFile, sorted_events, fps, accum_time)
% VISUALIZE_EVENT_OVERLAY Overlay ON/OFF events on original video frames.
%
% Inputs:
%   videoFile      - path to the original video (string)
%   sorted_events  - struct with fields t, x, y, p (from demo_array_design.m)
%   fps            - video frame rate
%   accum_time     - accumulation window in seconds (e.g., 0.05 = 50 ms)

    v = VideoReader(videoFile);
    frameCount = 0;

    figure('Position', [100,100,800,600]);

    while hasFrame(v)
        frame = readFrame(v);
        t_frame = frameCount / fps;
        t_start = t_frame;
        t_end   = t_frame + accum_time;

        % Select events that fall into this window
        idx = (double(sorted_events.t)*1e-6 >= t_start) & ...
              (double(sorted_events.t)*1e-6 < t_end);

        evs.x = sorted_events.x(idx);
        evs.y = sorted_events.y(idx);
        evs.p = sorted_events.p(idx);

        % Show frame
        imshow(frame); hold on;

        % Plot ON and OFF events
        onIdx  = evs.p == 1;
        offIdx = evs.p == -1;
        plot(evs.x(onIdx), evs.y(onIdx), 'r.', 'MarkerSize', 5);
        plot(evs.x(offIdx), evs.y(offIdx), 'b.', 'MarkerSize', 5);

        title(sprintf('Frame %d | Events in %.0f ms', frameCount, accum_time*1000));
        drawnow;

        hold off;
        frameCount = frameCount + 1;
    end
end
visualize_event_overlay('../Type2/TEST10_172_f1905_2k.mp4', ...
                        sorted_events, fps, 0.05);
