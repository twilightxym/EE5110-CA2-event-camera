function ax = visualize_event_overlay_formalapp(ax, video_data, sorted_events, video_path, accum_time, opts)
% VISUALIZE_EVENT_OVERLAY
% Overlay ON/OFF events on original color video frames, drawing into a target axes.
%
% USAGE:
%   % In App (mlapp), draw into app.UIAxes:
%   ax = visualize_event_overlay(app.UIAxes, video_data, sorted_events, seq_dir, 0.01);
%
%   % In script (standalone), let it create its own figure:
%   ax = visualize_event_overlay([], video_data, sorted_events, seq_dir, 0.01);
%
% INPUTS:
%   ax            - target axes to draw on; if empty, a new figure+axes will be created.
%   video_data    - struct from load_x4k_frames (fields: t_us, H, W, N)
%   sorted_events - struct with fields t, x, y, p (already time-sorted)
%   video_path    - path to original color video
%   accum_time    - accumulation window in seconds (e.g., 0.01 = 10 ms)
%   opts          - (optional) struct for runtime options:
%                   .onMarker       (default '.')
%                   .offMarker      (default '.')
%                   .onSize         (default 5)
%                   .offSize        (default 5)
%                   .onColor        (default [1 0 0])   % red
%                   .offColor       (default [0 0 1])   % blue
%                   .pause_s        (default 0)         % useful only in standalone
%                   .showTitle      (default true)
%                   .clearEachFrame (default true)
%
% OUTPUT:
%   ax            - the axes actually used for drawing (same handle as input ax, or newly created)
%
% NOTE:
%   - This function does not call 'figure' if a valid axes is provided.
%   - Designed to be App-friendly (no blocking pause unless opts.pause_s>0).

arguments
    ax {mustBeA(ax,{'matlab.graphics.axis.Axes','matlab.ui.control.UIAxes','matlab.graphics.axis.GeographicAxes'})}
    video_data (1,1) struct
    sorted_events (1,1) struct
    video_path (1,:) char
    accum_time (1,1) double {mustBePositive}
    opts.onMarker (1,1) string = "."
    opts.offMarker (1,1) string = "."
    opts.onSize (1,1) double {mustBePositive} = 5
    opts.offSize (1,1) double {mustBePositive} = 5
    opts.onColor (1,3) double {mustBeGreaterThanOrEqual(opts.onColor,0),mustBeLessThanOrEqual(opts.onColor,1)} = [1 0 0]
    opts.offColor (1,3) double {mustBeGreaterThanOrEqual(opts.offColor,0),mustBeLessThanOrEqual(opts.offColor,1)} = [0 0 1]
    opts.pause_s (1,1) double {mustBeGreaterThanOrEqual(opts.pause_s,0)} = 0.1
    opts.showTitle (1,1) logical = true
    opts.clearEachFrame (1,1) logical = true
end

% --- Unpack video data ---
t_us_video = video_data.t_us; % timestamps in µs
H          = video_data.H;
W          = video_data.W;
N          = video_data.N;

% --- Event timestamps (in µs) ---
event_time_us = double(sorted_events.t);
accum_time_us = accum_time * 1e6; % convert seconds → µs

% --- Prepare target axes ---
if isempty(ax) || ~ishandle(ax)
    % standalone fallback: create a new figure and axes
    fh = figure('Name','Event Overlay','Position',[100,100,900,600]);
    ax = axes('Parent', fh);
end

% ensure axis settings
axis(ax, 'ij');           % image coordinate: (1,1) at top-left
xlim(ax, [1 W]); ylim(ax, [1 H]);
set(ax, 'XTick',[], 'YTick',[]);
hold(ax, 'on');

% --- Open original color video ---
v = VideoReader(video_path);

% Pre-create base image and two scatter graphics for speed
baseImg = imagesc(ax, zeros(H,W,3,'uint8'));
hOn  = scatter(ax, NaN, NaN, opts.onSize, opts.onColor, opts.onMarker, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'flat');
hOff = scatter(ax, NaN, NaN, opts.offSize, opts.offColor, opts.offMarker, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'flat');

% % --- Loop over frames ---
% for k = 1:N
%     if ~hasFrame(v), break; end
%     frame_color = readFrame(v);
%
%     % Resize to match video_data if needed
%     if size(frame_color,1) ~= H || size(frame_color,2) ~= W
%         frame_color = imresize(frame_color, [H W]);
%     end
%
%     if opts.clearEachFrame
%         cla(ax);
%         % redraw base image and scatter holders
%         image(ax, frame_color); axis(ax,'image','ij'); hold(ax,'on');
%         hOn  = scatter(ax, NaN, NaN, opts.onSize,  opts.onColor,  opts.onMarker,  'filled');
%         hOff = scatter(ax, NaN, NaN, opts.offSize, opts.offColor, opts.offMarker, 'filled');
%     else
%         % update existing image CData
%         imHandle = findobj(ax, 'Type','image');
%         if isempty(imHandle)
%             image(ax, frame_color); axis(ax,'image','ij'); hold(ax,'on');
%         else
%             set(imHandle(1), 'CData', frame_color);
%         end
%     end
%
%     % Time window for events
%     t_start_us = t_us_video(k);
%     t_end_us   = min(t_start_us + accum_time_us, t_us_video(end));
%
%     % Select events in this window
%     idx = (event_time_us >= t_start_us) & (event_time_us < t_end_us);
%     if any(idx)
%         evx = double(sorted_events.x(idx));
%         evy = double(sorted_events.y(idx));
%         evp = sorted_events.p(idx);
%
%         onIdx  = evp ==  1;
%         offIdx = evp == -1;
%
%         if any(onIdx)
%             set(hOn, 'XData', evx(onIdx), 'YData', evy(onIdx));
%         else
%             set(hOn, 'XData', NaN, 'YData', NaN);
%         end
%
%         if any(offIdx)
%             set(hOff, 'XData', evx(offIdx), 'YData', evy(offIdx));
%         else
%             set(hOff, 'XData', NaN, 'YData', NaN);
%         end
%     else
%         % no events in window
%         set(hOn,  'XData', NaN, 'YData', NaN);
%         set(hOff, 'XData', NaN, 'YData', NaN);
%     end
%
%     if opts.showTitle
%         title(ax, sprintf('Frame %d | Events %.2f–%.2f ms', ...
%             k, double(t_start_us)/1000, double(t_end_us)/1000));
%     end
%
%     drawnow limitrate;    % App-friendly刷新
%     if opts.pause_s > 0   % Standalone展示时可用
%         pause(opts.pause_s);
%     end
% end

% 先清空一次，并把坐标系固定好
cla(ax);
set(ax,'YDir','reverse'); axis(ax,'image'); hold(ax,'on');
xlim(ax,[1 W]); ylim(ax,[1 H]);
set(ax,'XTick',[],'YTick',[]);

% 底图句柄
imH = image(ax, zeros(H,W,3,'uint8'));  % 先放个黑底，后面更新 CData

% 用 line 作为“点层”，更稳（比 scatter 更不挑剔）
hOn  = line(ax, NaN, NaN, 'LineStyle','none', 'Marker','.', 'MarkerSize', 1, 'Color', [1 0 0]);
hOff = line(ax, NaN, NaN, 'LineStyle','none', 'Marker','.', 'MarkerSize', 1, 'Color', [0 0 1]);
uistack(hOn,'top'); uistack(hOff,'top');  % 确保在最上层

% --- 循环 ---
for k = 1:N
    if ~hasFrame(v), break; end
    frame_color = readFrame(v);

    % Resize to match
    if size(frame_color,1) ~= H || size(frame_color,2) ~= W
        frame_color = imresize(frame_color, [H W]);
    end

    % 更新底图（不要 cla，这样不会把点层删掉）
    set(imH, 'CData', frame_color);

    % 时间窗
    t_start_us = t_us_video(k);
    t_end_us   = min(t_start_us + accum_time_us, t_us_video(end));

    % 选事件
    idx = (event_time_us >= t_start_us) & (event_time_us < t_end_us);

    if any(idx)
        evx = double(sorted_events.x(idx));
        evy = double(sorted_events.y(idx));
        evp = sorted_events.p(idx);

        % 越界过滤（防止 0 基坐标等）
        valid = (evx >= 1 & evx <= W & evy >= 1 & evy <= H);
        evx = evx(valid); evy = evy(valid); evp = evp(valid);

        % 极性容错：>0 当 ON，<=0 当 OFF
        onIdx  = (evp > 0);
        offIdx = ~onIdx;

        if any(onIdx)
            set(hOn,  'XData', evx(onIdx),  'YData', evy(onIdx));
        else
            set(hOn,  'XData', NaN, 'YData', NaN);
        end

        if any(offIdx)
            set(hOff, 'XData', evx(offIdx), 'YData', evy(offIdx));
        else
            set(hOff, 'XData', NaN, 'YData', NaN);
        end

        % 确保点层始终在最上
        uistack(hOn,'top'); uistack(hOff,'top');
    else
        set(hOn,'XData',NaN,'YData',NaN);
        set(hOff,'XData',NaN,'YData',NaN);
    end

    if opts.showTitle
        title(ax, sprintf('Frame %d | Events %.2f–%.2f ms', ...
            k, double(t_start_us)/1000, double(t_end_us)/1000));
    end

    % 调试：每 20 帧输出一次事件数（可删）
    % if mod(k,20)==0
    %     fprintf('[overlay] k=%d, events in window=%d\n', k, sum(idx));
    % end

    drawnow limitrate;
    if opts.pause_s>0, pause(opts.pause_s); end
end

hold(ax,'off');
end