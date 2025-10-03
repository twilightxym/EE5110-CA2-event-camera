function [t_events_us, p_events] = pixel_events(L_log, t_us, C, opts)
% PIXEL_EVENTS  event generation (time-resolution scanning).
%   Check log-illuminance on a fixed time grid of Δt, 
%   trigger when |L_log(t_g) - Lref| >= C_eff, then update Lref immediately.
%
% Inputs:
%   L_log : [N x 1] log-intensity (monotonic in time)
%   t_us  : [N x 1] timestamps in microseconds (int64/double, strictly increasing)
%   C     : positive contrast threshold (scalar), e.g., 0.15
%   opts  :
%     .Lref0               initial reference log-intensity (default L_log(1))
%     .refractory_us       refractory period in microseconds (default 0)
%     .threshold_jitter    std of Gaussian jitter for C per-event (default 0)
%     .time_resolution_us  REQUIRED > 0, time resolution Δt for grid scanning
%
% Outputs:
%   t_events_us : [M x 1] event timestamps (microseconds, aligned to grid)
%   p_events    : [M x 1] event polarities (+1 or -1)

    if nargin < 4, opts = struct; end
    if ~isfield(opts,'Lref0') || isempty(opts.Lref0), opts.Lref0 = L_log(1); end
    if ~isfield(opts,'refractory_us'),                opts.refractory_us = 0; end
    if ~isfield(opts,'threshold_jitter'),             opts.threshold_jitter = 0; end
    if ~isfield(opts,'time_resolution_us') || opts.time_resolution_us <= 0
        error('pixel_events: opts.time_resolution_us must be > 0 for lecture mode.');
    end

    % ---- init ----
    Lref = opts.Lref0;
    t_events_us = zeros(0,1,'like',t_us);
    p_events    = zeros(0,1);

    % ---- build uniform time grid ----
    dt_us = double(opts.time_resolution_us);
    t0d   = double(t_us(1));
    tNd   = double(t_us(end));

    t_grid = t0d:dt_us:tNd;
    if t_grid(end) < tNd
        t_grid(end+1) = tNd; %#ok<AGROW>
    end

    % ---- interpolate log intensity onto grid ----
    Lg = interp1(double(t_us(:)), double(L_log(:)), t_grid, 'linear');

    % ---- scan grid ----
    last_t_event = [];
    for i = 1:numel(t_grid)
        tg  = t_grid(i);
        res = Lg(i) - Lref;

        % one check per tick: if |res| ≥ threshold → trigger once
        C_eff = C + opts.threshold_jitter * randn;

        if res >= C_eff
            pol = +1;
        elseif res <= -C_eff
            pol = -1;
        else
            continue; % no event at this tick
        end

        % refractory check
        tg_us = cast(tg, 'like', t_us);
        if ~isempty(last_t_event)
            if (double(tg_us) - double(last_t_event)) < double(opts.refractory_us)
                continue; % within refractory, drop
            end
        end

        % commit event
        t_events_us(end+1,1) = tg_us; %#ok<AGROW>
        p_events(end+1,1)    = pol;   %#ok<AGROW>
        last_t_event = tg_us;

        % update reference
        Lref = Lref + pol * C_eff;
    end
end
