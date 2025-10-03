function [t_events_us, p_events] = pixel_events_pro(L_log, t_us, C, opts)
% PIXEL_EVENTS_PRO  Generate events for a single pixel time series (enhanced).
%   Events are generated strictly on a fixed time grid of Δt=opts.time_resolution_us.
%
% Inputs:
%   L_log : [N x 1] log-intensity sequence (strict time order).
%   t_us  : [N x 1] timestamps in microseconds.
%   C     : positive contrast threshold.
%   opts  : struct (advanced fields, see below).
%     % ---------- Baseline ----------
%     .time_resolution_us               timestamp resolution (required > 0)
%     .Lref0                 initial reference log-intensity (default L_log(1))
%     .refractory_us         min spacing between two events (default 0)
%     % ---------- Threshold noise ----------
%     .threshold_jitter      per-event Gaussian std on C (default 0)
%     .threshold_sigma_px    per-pixel Gaussian std on C (draw once, default 0)
%     .threshold_min         floor for effective thresholds (default 1e-6)
%     % ---------- Intensity (for advanced models) ----------
%     .I_norm                [N x 1] normalized intensity in [0,1] (optional)
%     .I_from_log_eps        epsilon used to recover intensity from L_log (default 0)
%     .I_norm_minmax         if true (default), I_norm = minmax-normalize(exp(L_log))
%     % ---------- Photoreceptor (finite bandwidth) ----------
%     .pr.enable             true/false (default false)
%     .pr.tau_dark_us        time-constant in dark (large, e.g., 5000)
%     .pr.tau_bright_us      time-constant in bright (small, e.g., 500)
%     % tau(I) = tau_bright + (tau_dark - tau_bright) * (1 - I_norm)
%     % ---------- Leak events ----------
%     .leak.enable           true/false (default false)
%     .leak.rate_dark_hz     leak rate in dark (e.g., 0 to 5)
%     .leak.rate_bright_hz   leak rate in bright (e.g., 0 to 20)
%     .leak.on_prob_base     base prob of ON polarity in [0,1] (default 0.5)
%     .leak.on_prob_slope    slope wrt I_norm (default 0, keep 0.5 constant)
%     .leak.update_Lref      whether leak updates Lref (default false)
%     % ---------- Timing jitter ----------
%     .timing_jitter.enable  true/false (default false)
%     .timing_jitter.dark_us std of time noise in dark
%     .timing_jitter.bright_us std of time noise in bright
%
% Outputs:
%   t_events_us : [M x 1] event timestamps (on Δt grid)
%   p_events    : [M x 1] event polarities (+1 or -1)

% -------- defaults --------
if nargin < 4, opts = struct; end
if ~isfield(opts,'Lref0') || isempty(opts.Lref0), opts.Lref0 = L_log(1); end
if ~isfield(opts,'refractory_us'),     opts.refractory_us = 0;    end
if ~isfield(opts,'threshold_jitter'),  opts.threshold_jitter = 0; end
if ~isfield(opts,'threshold_sigma_px'),opts.threshold_sigma_px = 0;end
if ~isfield(opts,'threshold_min'),     opts.threshold_min = 1e-6; end
if ~isfield(opts,'I_from_log_eps'),    opts.I_from_log_eps = 0;   end
if ~isfield(opts,'I_norm_minmax'),     opts.I_norm_minmax = true; end
if ~isfield(opts,'time_resolution_us') || opts.time_resolution_us <= 0
    error('opts.time_resolution_us must be >0 for Δt scanning mode.');
end

% Photoreceptor defaults
if ~isfield(opts,'pr'), opts.pr = struct; end
if ~isfield(opts.pr,'enable'),         opts.pr.enable = false; end
if ~isfield(opts.pr,'tau_dark_us'),    opts.pr.tau_dark_us = 5000; end
if ~isfield(opts.pr,'tau_bright_us'),  opts.pr.tau_bright_us = 500; end

% Leak defaults
if ~isfield(opts,'leak'), opts.leak = struct; end
if ~isfield(opts.leak,'enable'),          opts.leak.enable = false; end
if ~isfield(opts.leak,'rate_dark_hz'),    opts.leak.rate_dark_hz = 0; end
if ~isfield(opts.leak,'rate_bright_hz'),  opts.leak.rate_bright_hz = 0; end
if ~isfield(opts.leak,'on_prob_base'),    opts.leak.on_prob_base = 0.5; end
if ~isfield(opts.leak,'on_prob_slope'),   opts.leak.on_prob_slope = 0.0; end
if ~isfield(opts.leak,'update_Lref'),     opts.leak.update_Lref = false; end

% Timing jitter defaults
if ~isfield(opts,'timing_jitter'), opts.timing_jitter = struct; end
if ~isfield(opts.timing_jitter,'enable'),    opts.timing_jitter.enable = false; end
if ~isfield(opts.timing_jitter,'dark_us'),   opts.timing_jitter.dark_us = 0; end
if ~isfield(opts.timing_jitter,'bright_us'), opts.timing_jitter.bright_us = 0; end

% -------- intensity normalization --------
N = numel(L_log);
if isfield(opts,'I_norm') && ~isempty(opts.I_norm)
    I_norm = opts.I_norm(:);
    if numel(I_norm) ~= N, error('opts.I_norm length mismatch.'); end
else
    I_lin = max(exp(L_log) - opts.I_from_log_eps, 0);
    if opts.I_norm_minmax
        mn = min(I_lin); mx = max(I_lin);
        I_norm = (I_lin - mn) ./ max(mx - mn, eps);
    else
        I_norm = max(min(I_lin,1),0);
    end
end

% -------- photoreceptor dynamics --------
if opts.pr.enable
    tau_dark   = max(opts.pr.tau_dark_us,   1e-3);
    tau_bright = max(opts.pr.tau_bright_us, 1e-3);
    tau_us = tau_bright + (tau_dark - tau_bright) .* (1 - I_norm);

    L_used = zeros(N,1);
    L_used(1) = L_log(1);
    for k = 2:N
        dt = double(t_us(k) - t_us(k-1));
        a = 1 - exp(- max(dt,0) / max(tau_us(k-1),1e-12) );
        L_used(k) = L_used(k-1) + a*(L_log(k) - L_used(k-1));
    end
else
    L_used = L_log;
end

% -------- per-pixel threshold variation --------
C_px = C + opts.threshold_sigma_px * randn;
C_px = max(C_px, opts.threshold_min);

% -------- build Δt grid --------
dt_us = double(opts.time_resolution_us);
t0d   = double(t_us(1));
tNd   = double(t_us(end));
t_grid = t0d:dt_us:tNd;
if t_grid(end) < tNd
    t_grid(end+1) = tNd; %#ok<AGROW>
end

% Interpolate signals on grid
Lg = interp1(double(t_us), double(L_used), t_grid, 'linear');
Ig = interp1(double(t_us), double(I_norm), t_grid, 'linear');

% -------- state init --------
Lref = opts.Lref0;
t_events_us = zeros(0,1,'like',t_us);
p_events    = zeros(0,1);
last_t_event = [];

% -------- main Δt scanning loop --------
for i = 1:numel(t_grid)
    tg  = t_grid(i);
    res = Lg(i) - Lref;

    % check once per tick
    C_eff = C_px + opts.threshold_jitter*randn;
    C_eff = max(C_eff, opts.threshold_min);

    pol = 0;
    if res >= C_eff, pol = +1; end
    if res <= -C_eff, pol = -1; end

    if pol ~= 0
        tg_us = cast(tg, 'like', t_us);

        % refractory check
        if ~isempty(last_t_event)
            if (double(tg_us) - double(last_t_event)) < double(opts.refractory_us)
                % drop this event
                continue;
            end
        end

        % optional timing jitter
        if opts.timing_jitter.enable
            sig_us = opts.timing_jitter.dark_us + ...
                     (opts.timing_jitter.bright_us - opts.timing_jitter.dark_us) * Ig(i);
            if sig_us > 0
                tg_us = tg_us + cast(randn*sig_us, 'like', t_us);
            end
        end

        % commit event
        t_events_us(end+1,1) = tg_us; %#ok<AGROW>
        p_events(end+1,1)    = pol;   %#ok<AGROW>
        last_t_event = tg_us;

        % update reference
        Lref = Lref + pol*C_eff;
    end

    % leak events (optional)
    if opts.leak.enable
        % expected number in dt interval
        dt_s = dt_us * 1e-6;
        rate_hz = opts.leak.rate_dark_hz + ...
                  (opts.leak.rate_bright_hz - opts.leak.rate_dark_hz) * Ig(i);
        rate_hz = max(rate_hz,0);
        lambda = rate_hz * dt_s;
        n_leak = poissrnd(lambda);

        for ii = 1:n_leak
            t_try = tg + rand*dt_us; % uniform within tick

            % optional timing jitter
            if opts.timing_jitter.enable
                sig_us = opts.timing_jitter.dark_us + ...
                         (opts.timing_jitter.bright_us - opts.timing_jitter.dark_us) * Ig(i);
                if sig_us > 0
                    t_try = t_try + randn*sig_us;
                end
            end

            t_try = cast(t_try, 'like', t_us);
            if ~isempty(last_t_event)
                if (double(t_try) - double(last_t_event)) < double(opts.refractory_us)
                    continue;
                end
            end

            % commit leak event
            t_events_us(end+1,1) = t_try; %#ok<AGROW>
            pol_leak = (rand < (opts.leak.on_prob_base + opts.leak.on_prob_slope*(Ig(i)-0.5)));
            pol_leak = 2*pol_leak - 1;
            p_events(end+1,1)    = pol_leak; %#ok<AGROW>
            last_t_event = t_try;

            if opts.leak.update_Lref
                Lref = Lref + pol_leak*C_px;
            end
        end
    end
end
end
