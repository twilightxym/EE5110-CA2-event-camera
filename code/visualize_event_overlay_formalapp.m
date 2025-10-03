function ax = visualize_event_overlay_formalapp(ax, video_data, sorted_events, video_path, accum_time, opts)
% VISUALIZE_EVENT_OVERLAY_FORMALAPP
% Overlay ON/OFF events on original color video frames into a target axes.
% Adaptive accumulation window with sub-window refresh if needed.
%
% USAGE (in mlapp):
%   ax = visualize_event_overlay_formalapp(app.UIAxes, video_data, sorted_events, video_path, 1e-5);
%
% USAGE (standalone script):
%   ax = visualize_event_overlay_formalapp([], video_data, sorted_events, video_path, 1e-5);
%
% INPUTS:
%   ax            - target axes; if empty, create a new figure
%   video_data    - struct (fields: t_us, H, W, N)
%   sorted_events - struct with fields t, x, y, p (time-sorted)
%   video_path    - path to original color video
%   accum_time    - accumulation window in seconds
%   opts          - struct with optional fields:
%       .onMarker, .offMarker  (default '.')
%       .onSize, .offSize      (default 5)
%       .onColor, .offColor    (default [1 0 0], [0 0 1])
%       .pause_s               (default 0, only used standalone)
%       .showTitle             (default true)
%
% OUTPUT:
%   ax            - handle of axes used for drawing
%
% RULE:
%   - If accum_time < frame_dt: multiple sub-windows per frame, background fixed
%   - If accum_time >= frame_dt: [t_k, t_{k+1}), no future events shown

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
    opts.onColor (1,3) double = [1 0 0]
    opts.offColor (1,3) double = [0 0 1]
    opts.pause_s (1,1) double {mustBeNonnegative} = 0.1 % visualization pause time
    opts.showTitle (1,1) logical = true
end

% --- unpack video data ---
t_us_video = double(video_data.t_us(:));
H = video_data.H; W = video_data.W; N = video_data.N;

% --- events ---
event_time_us = double(sorted_events.t(:));
ex = double(sorted_events.x(:));
ey = double(sorted_events.y(:));
ep = sorted_events.p(:);
M = numel(event_time_us);

accum_time_us = accum_time*1e6;
if N>1, frame_dt_us_vec = diff(t_us_video); else, frame_dt_us_vec = []; end

% --- prepare axes ---
if isempty(ax) || ~ishandle(ax)
    fh = figure('Name','Event Overlay','Position',[100,100,900,600]);
    ax = axes('Parent',fh);
end
cla(ax);
axis(ax,'ij'); xlim(ax,[1 W]); ylim(ax,[1 H]);
set(ax,'XTick',[],'YTick',[]);

% --- video reader ---
v = VideoReader(video_path);

% --- graphics handles (reuse) ---
imH  = image(ax, zeros(H,W,3,'uint8'));
hOn  = line(ax, NaN, NaN, 'LineStyle','none', 'Marker',opts.onMarker, ...
            'MarkerSize',opts.onSize, 'Color', opts.onColor);
hOff = line(ax, NaN, NaN, 'LineStyle','none', 'Marker',opts.offMarker, ...
            'MarkerSize',opts.offSize, 'Color', opts.offColor);
uistack(hOn,'top'); uistack(hOff,'top');

% --- double pointers for events ---
i0=1; i1=1;

for k = 1:N
    if ~hasFrame(v), break; end
    frame_color = readFrame(v);
    if size(frame_color,1)~=H || size(frame_color,2)~=W
        frame_color = imresize(frame_color,[H W]);
    end
    set(imH,'CData',frame_color);

    t_start = t_us_video(k);
    if k<N
        frame_dt = frame_dt_us_vec(k);
        t_next   = t_us_video(k+1);

        if accum_time_us < frame_dt
            % subdivide current frame into multiple sub-windows
            num_sub = ceil(frame_dt/accum_time_us);
            edges = t_start + (0:num_sub)*accum_time_us;
            edges(end) = t_next; % align exactly
        else
            num_sub = 1;
            edges = [t_start t_next];
        end
    else
        % last frame
        num_sub = 1;
        edges = [t_start min(t_start+accum_time_us, t_us_video(end))];
    end

    % ---- sub-window loop ----
    for s=1:num_sub
        t0 = edges(s); t1 = edges(s+1);

        % advance pointers
        while i0<=M && event_time_us(i0)<t0, i0=i0+1; end
        if i1<i0, i1=i0; end
        while i1<=M && event_time_us(i1)<t1, i1=i1+1; end

        if i0<i1
            idx = i0:(i1-1);
            xw=ex(idx); yw=ey(idx); pw=ep(idx);
        else
            xw=[]; yw=[]; pw=[];
        end

        set(hOn,'XData',xw(pw==1),'YData',yw(pw==1));
        set(hOff,'XData',xw(pw==-1),'YData',yw(pw==-1));

        if opts.showTitle
            title(ax,sprintf('Frame %d | Sub %d/%d | %.3f–%.3f ms | events=%d',...
                k,s,num_sub,t0/1000,t1/1000,numel(xw)));
        end

        drawnow limitrate;
        if opts.pause_s>0, pause(opts.pause_s); end
    end
end
end

% function ax = visualize_event_overlay_formalapp(ax, video_data, sorted_events, video_path, accum_time, opts)
% % VISUALIZE_EVENT_OVERLAY
% % Overlay ON/OFF events on original color video frames, drawing into a target axes.
% %
% % USAGE:
% %   % In App (mlapp), draw into app.UIAxes:
% %   ax = visualize_event_overlay(app.UIAxes, video_data, sorted_events, seq_dir, 0.01);
% %
% %   % In script (standalone), let it create its own figure:
% %   ax = visualize_event_overlay([], video_data, sorted_events, seq_dir, 0.01);
% %
% % INPUTS:
% %   ax            - target axes to draw on; if empty, a new figure+axes will be created.
% %   video_data    - struct from load_x4k_frames (fields: t_us, H, W, N)
% %   sorted_events - struct with fields t, x, y, p (already time-sorted)
% %   video_path    - path to original color video
% %   accum_time    - accumulation window in seconds (e.g., 0.01 = 10 ms)
% %   opts          - (optional) struct for runtime options:
% %                   .onMarker       (default '.')
% %                   .offMarker      (default '.')
% %                   .onSize         (default 5)
% %                   .offSize        (default 5)
% %                   .onColor        (default [1 0 0])   % red
% %                   .offColor       (default [0 0 1])   % blue
% %                   .pause_s        (default 0)         % useful only in standalone
% %                   .showTitle      (default true)
% %                   .clearEachFrame (default true)
% %
% % OUTPUT:
% %   ax            - the axes actually used for drawing (same handle as input ax, or newly created)
% %
% % NOTE:
% %   - This function does not call 'figure' if a valid axes is provided.
% %   - Designed to be App-friendly (no blocking pause unless opts.pause_s>0).
% 
% arguments
%     ax {mustBeA(ax,{'matlab.graphics.axis.Axes','matlab.ui.control.UIAxes','matlab.graphics.axis.GeographicAxes'})}
%     video_data (1,1) struct
%     sorted_events (1,1) struct
%     video_path (1,:) char
%     accum_time (1,1) double {mustBePositive}
%     opts.onMarker (1,1) string = "."
%     opts.offMarker (1,1) string = "."
%     opts.onSize (1,1) double {mustBePositive} = 5
%     opts.offSize (1,1) double {mustBePositive} = 5
%     opts.onColor (1,3) double {mustBeGreaterThanOrEqual(opts.onColor,0),mustBeLessThanOrEqual(opts.onColor,1)} = [1 0 0]
%     opts.offColor (1,3) double {mustBeGreaterThanOrEqual(opts.offColor,0),mustBeLessThanOrEqual(opts.offColor,1)} = [0 0 1]
%     opts.pause_s (1,1) double {mustBeGreaterThanOrEqual(opts.pause_s,0)} = 0.1
%     opts.showTitle (1,1) logical = true
%     opts.clearEachFrame (1,1) logical = true
% end
% 
% % --- Unpack video data ---
% t_us_video = video_data.t_us; % timestamps in µs
% H          = video_data.H;
% W          = video_data.W;
% N          = video_data.N;
% 
% % --- Event timestamps (in µs) ---
% event_time_us = double(sorted_events.t);
% accum_time_us = accum_time * 1e6; % convert seconds → µs
% 
% % --- Prepare target axes ---
% if isempty(ax) || ~ishandle(ax)
%     % standalone fallback: create a new figure and axes
%     fh = figure('Name','Event Overlay','Position',[100,100,900,600]);
%     ax = axes('Parent', fh);
% end
% 
% % ensure axis settings
% axis(ax, 'ij');           % image coordinate: (1,1) at top-left
% xlim(ax, [1 W]); ylim(ax, [1 H]);
% set(ax, 'XTick',[], 'YTick',[]);
% hold(ax, 'on');
% 
% % --- Open original color video ---
% v = VideoReader(video_path);
% 
% % Pre-create base image and two scatter graphics for speed
% baseImg = imagesc(ax, zeros(H,W,3,'uint8'));
% hOn  = scatter(ax, NaN, NaN, opts.onSize, opts.onColor, opts.onMarker, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'flat');
% hOff = scatter(ax, NaN, NaN, opts.offSize, opts.offColor, opts.offMarker, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'flat');
% 
% cla(ax);
% set(ax,'YDir','reverse'); axis(ax,'image'); hold(ax,'on');
% xlim(ax,[1 W]); ylim(ax,[1 H]);
% set(ax,'XTick',[],'YTick',[]);
% 
% imH = image(ax, zeros(H,W,3,'uint8'));
% 
% hOn  = line(ax, NaN, NaN, 'LineStyle','none', 'Marker','.', 'MarkerSize', 1, 'Color', [1 0 0]);
% hOff = line(ax, NaN, NaN, 'LineStyle','none', 'Marker','.', 'MarkerSize', 1, 'Color', [0 0 1]);
% uistack(hOn,'top'); uistack(hOff,'top');  % 确保在最上层
% 
% % --- loop ---
% for k = 1:N
%     if ~hasFrame(v), break; end
%     frame_color = readFrame(v);
% 
%     % Resize to match
%     if size(frame_color,1) ~= H || size(frame_color,2) ~= W
%         frame_color = imresize(frame_color, [H W]);
%     end
% 
%     set(imH, 'CData', frame_color);
% 
%     t_start_us = t_us_video(k);
%     t_end_us   = min(t_start_us + accum_time_us, t_us_video(end));
% 
%     idx = (event_time_us >= t_start_us) & (event_time_us < t_end_us);
% 
%     if any(idx)
%         evx = double(sorted_events.x(idx));
%         evy = double(sorted_events.y(idx));
%         evp = sorted_events.p(idx);
% 
%         valid = (evx >= 1 & evx <= W & evy >= 1 & evy <= H);
%         evx = evx(valid); evy = evy(valid); evp = evp(valid);
% 
%         % p：>0:ON，<=0:OFF
%         onIdx  = (evp > 0);
%         offIdx = ~onIdx;
% 
%         if any(onIdx)
%             set(hOn,  'XData', evx(onIdx),  'YData', evy(onIdx));
%         else
%             set(hOn,  'XData', NaN, 'YData', NaN);
%         end
% 
%         if any(offIdx)
%             set(hOff, 'XData', evx(offIdx), 'YData', evy(offIdx));
%         else
%             set(hOff, 'XData', NaN, 'YData', NaN);
%         end
% 
%         uistack(hOn,'top'); uistack(hOff,'top');
%     else
%         set(hOn,'XData',NaN,'YData',NaN);
%         set(hOff,'XData',NaN,'YData',NaN);
%     end
% 
%     if opts.showTitle
%         title(ax, sprintf('Frame %d | Events %.2f–%.2f ms', ...
%             k, double(t_start_us)/1000, double(t_end_us)/1000));
%     end
% 
% 
%     drawnow limitrate;
%     if opts.pause_s>0, pause(opts.pause_s); end
% end
% 
% hold(ax,'off');
% end