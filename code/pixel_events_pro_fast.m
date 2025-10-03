function [t_events_us, p_events] = pixel_events_pro_fast(L_log, t_us, C, opts, I_norm_in)
% PIXEL_EVENTS_PRO_FAST  Single-pixel event generation (fast, modular).
% - Accepts optional linear normalized intensity I_norm_in (Nx1 in [0,1])
% - Single-precision internal compute; outputs match t_us type
% - Short-circuits disabled modules; sets defaults only if needed
% 
% Usage:
%   [t,p] = pixel_events_pro_fast(L_log, t_us, C, opts);            % auto I_norm (when needed)
%   [t,p] = pixel_events_pro_fast(L_log, t_us, C, opts, I_norm);    % external I_norm
% 
% I/O:
%   L_log         [N x 1]   log-intensity (monotone time)
%   t_us          [N x 1]   timestamps (int64/double, strictly increasing)
%   C             scalar    base contrast threshold (>0)
%   opts          struct    optional modules, see below
%   I_norm_in     [N x 1]   linear intensity normalized to [0,1] (optional)
% 
%   t_events_us   [M x 1]  event timestamps (microsecond)
%   p_events      [M x 1]  event polarities (+1/-1)
% 
% Supported:
%   Baseline:
%     .time_resolution_us               timestamp resolution (required > 0)
%     .Lref0,                           initial reference log-intensity (default L_log(1))
%     .refractory_us,                   min spacing between two events (default 0)
%     .threshold_min,                   floor for effective thresholds (default 1e-6)
%     .threshold_jitter,                per-event Gaussian std on C (default 0)
%     .threshold_sigma_px,              per-pixel Gaussian std on C (default 0)
%   Photoreceptor (finite bandwidth):
%     .pr.enable,                       true/false (default false)
%     .pr.tau_dark_us,                  time-constant in dark (large, e.g., 5000)
%     .pr.tau_bright_us,                time-constant in bright (small, e.g., 500)
%     tau(I) = tau_bright + (tau_dark - tau_bright) * (1 - I_norm)
%   Leak events:
%     .leak.enable,                     true/false (default false)
%     .leak.rate_dark_hz,               leak rate in dark (e.g., 0 to 5)
%     .leak.rate_bright_hz,             leak rate in bright (e.g., 0 to 20)
%     .leak.on_prob_base,               base prob of ON polarity in [0,1] (default 0.5)
%     .leak.on_prob_slope,              slope wrt I_norm (default 0, keep 0.5 constant)
%     .leak.update_Lref,                whether leak updates Lref (default false)
%   Timing jitter:
%     .timing_jitter.enable,            true/false (default false)
%     .timing_jitter.dark_us,           std of time noise in dark
%     .timing_jitter.bright_us,         std of time noise in bright

% ----------- defaults -----------
if nargin < 4 || isempty(opts), opts = struct; end
if ~isfield(opts,'Lref0') || isempty(opts.Lref0), opts.Lref0 = L_log(1); end
if ~isfield(opts,'threshold_min'),  opts.threshold_min = 1e-6; end
if ~isfield(opts,'refractory_us'),  opts.refractory_us = 0;    end
if ~isfield(opts,'time_resolution_us') || opts.time_resolution_us <= 0
    error('opts.time_resolution_us must be >0 for Δt scanning mode.');
end

% flags
use_pr    = isfield(opts,'pr') && isfield(opts.pr,'enable') && logical(opts.pr.enable);
use_leak  = isfield(opts,'leak') && isfield(opts.leak,'enable') && logical(opts.leak.enable);
use_tjit  = isfield(opts,'timing_jitter') && isfield(opts.timing_jitter,'enable') && logical(opts.timing_jitter.enable);
has_Cpx   = isfield(opts,'threshold_sigma_px') && opts.threshold_sigma_px ~= 0;
has_Cjit  = isfield(opts,'threshold_jitter')   && opts.threshold_jitter   ~= 0;

% cast
L_log_s = single(L_log(:));
N       = numel(L_log_s);

% ----------- I_norm (only if needed) -----------
need_Inorm = use_pr || use_leak || use_tjit;
if need_Inorm
    if nargin >= 5 && ~isempty(I_norm_in)
        I_norm = single(I_norm_in(:));
        if numel(I_norm) ~= N, error('I_norm_in length mismatch.'); end
        I_norm = max(min(I_norm,1),0);
    else
        I_lin = exp(L_log_s);
        mn = min(I_lin); mx = max(I_lin);
        I_norm = (I_lin - mn) ./ max(mx - mn, eps('single'));
        I_norm = max(min(I_norm,1),0);
    end
else
    I_norm = [];
end

% ----------- photoreceptor dynamics -----------
if use_pr
    if ~isfield(opts.pr,'tau_dark_us'),   opts.pr.tau_dark_us   = 5000; end
    if ~isfield(opts.pr,'tau_bright_us'), opts.pr.tau_bright_us = 500;  end
    tau_dark   = single(max(opts.pr.tau_dark_us,   1e-3));
    tau_bright = single(max(opts.pr.tau_bright_us, 1e-3));
    tau_us = tau_bright + (tau_dark - tau_bright) .* (1 - I_norm);

    L_used = zeros(N,1,'single');
    L_used(1) = L_log_s(1);
    for k = 2:N
        dt = single(double(t_us(k) - t_us(k-1)));
        dt = max(dt,0);
        a  = 1 - exp(- dt / max(tau_us(k-1),1e-12));
        L_used(k) = L_used(k-1) + a*(L_log_s(k) - L_used(k-1));
    end
else
    L_used = L_log_s;
end

% ----------- thresholds -----------
C_s   = single(C);
C_min = single(opts.threshold_min);
if has_Cpx
    C_px = max(C_s + single(opts.threshold_sigma_px)*randn('single'), C_min);
else
    C_px = max(C_s, C_min);
end

% timing jitter params
tj_dark = single(0); tj_bright = single(0);
if use_tjit
    if ~isfield(opts.timing_jitter,'dark_us'),   opts.timing_jitter.dark_us   = 0; end
    if ~isfield(opts.timing_jitter,'bright_us'), opts.timing_jitter.bright_us = 0; end
    tj_dark   = single(opts.timing_jitter.dark_us);
    tj_bright = single(opts.timing_jitter.bright_us);
end

% ----------- build Δt grid -----------
dt_us = double(opts.time_resolution_us);
t0d   = double(t_us(1));
tNd   = double(t_us(end));
t_grid = t0d:dt_us:tNd;
if t_grid(end) < tNd
    t_grid(end+1) = tNd;
end

% interpolate signals onto grid
Lg = interp1(double(t_us), double(L_used), t_grid, 'linear');
if need_Inorm
    Ig = interp1(double(t_us), double(I_norm), t_grid, 'linear');
else
    Ig = [];
end

% ----------- outputs & state -----------
t_events_us = zeros(0,1,'like',t_us);
p_events    = zeros(0,1);
last_t      = -Inf;
Lref        = single(opts.Lref0);

% ----------- scanning loop -----------
for i = 1:numel(t_grid)
    tg = t_grid(i);
    res = single(Lg(i) - double(Lref));

    % threshold check (once per tick)
    if has_Cjit
        C_eff = max(C_px + single(opts.threshold_jitter)*randn('single'), C_min);
    else
        C_eff = C_px;
    end

    pol = 0;
    if res >= C_eff, pol = +1; end
    if res <= -C_eff, pol = -1; end

    if pol ~= 0
        tg_us = cast(tg,'like',t_us);

        % refractory
        if (double(tg_us) - last_t) >= double(opts.refractory_us)
            % optional timing jitter
            if use_tjit
                sig = double(tj_dark);
                if need_Inorm
                    sig = double(tj_dark) + (double(tj_bright)-double(tj_dark))*double(Ig(i));
                end
                if sig ~= 0
                    tg_us = tg_us + cast(randn*sig,'like',t_us);
                end
            end

            t_events_us(end+1,1) = tg_us; %#ok<AGROW>
            p_events(end+1,1)    = pol;   %#ok<AGROW>
            last_t = double(tg_us);
            Lref = Lref + single(pol)*C_eff;
        end
    end

    % leak events
    if use_leak
        if ~isfield(opts.leak,'rate_dark_hz'),    opts.leak.rate_dark_hz = 0; end
        if ~isfield(opts.leak,'rate_bright_hz'),  opts.leak.rate_bright_hz = 0; end
        if ~isfield(opts.leak,'on_prob_base'),    opts.leak.on_prob_base = 0.5; end
        if ~isfield(opts.leak,'on_prob_slope'),   opts.leak.on_prob_slope = 0.0; end
        if ~isfield(opts.leak,'update_Lref'),     opts.leak.update_Lref = false; end

        dt_s = dt_us*1e-6;
        In_avg = 0;
        if need_Inorm, In_avg = double(Ig(i)); end
        rate_hz = double(opts.leak.rate_dark_hz) + ...
                  (double(opts.leak.rate_bright_hz)-double(opts.leak.rate_dark_hz))*In_avg;
        rate_hz = max(rate_hz,0);
        lambda = rate_hz*dt_s;

        if lambda > 1e-6
            n_leak = poissrnd(lambda);
            for ii = 1:n_leak
                t_try = tg + rand*dt_us;
                if use_tjit
                    sig = double(tj_dark);
                    if need_Inorm
                        sig = double(tj_dark) + (double(tj_bright)-double(tj_dark))*In_avg;
                    end
                    if sig ~= 0
                        t_try = t_try + randn*sig;
                    end
                end
                t_try = cast(t_try,'like',t_us);
                if (double(t_try) - last_t) < double(opts.refractory_us)
                    continue;
                end
                pol_leak = (rand < (opts.leak.on_prob_base + opts.leak.on_prob_slope*(In_avg-0.5)));
                pol_leak = 2*pol_leak - 1;
                t_events_us(end+1,1) = t_try; %#ok<AGROW>
                p_events(end+1,1)    = pol_leak; %#ok<AGROW>
                last_t = double(t_try);
                if opts.leak.update_Lref
                    Lref = Lref + single(pol_leak)*C_px;
                end
            end
        end
    end
end
end
