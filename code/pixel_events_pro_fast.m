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
%
% Compatibility:
%   When given the same inputs/flags and random seed, results match pixel_events_pro
%   up to single-precision tolerance.

% ----------- light defaults (only generic) -----------
if nargin < 4 || isempty(opts), opts = struct; end
if ~isfield(opts,'Lref0') || isempty(opts.Lref0), opts.Lref0 = L_log(1); end
if ~isfield(opts,'threshold_min'),      opts.threshold_min = 1e-6; end
if ~isfield(opts,'refractory_us'),      opts.refractory_us = 0;    end


% Flags (decide modules once)
use_pr    = isfield(opts,'pr')          && isfield(opts.pr,'enable')          && logical(opts.pr.enable);
use_leak  = isfield(opts,'leak')        && isfield(opts.leak,'enable')        && logical(opts.leak.enable);
use_tjit  = isfield(opts,'timing_jitter') && isfield(opts.timing_jitter,'enable') && logical(opts.timing_jitter.enable);

% Per-pixel / per-event threshold noise flags
has_Cpx   = isfield(opts,'threshold_sigma_px') && opts.threshold_sigma_px ~= 0;
has_Cjit  = isfield(opts,'threshold_jitter')   && opts.threshold_jitter   ~= 0;

% Need I_norm if any of these modules are on
need_Inorm = use_pr || use_leak || use_tjit;

% ----------- cast to single for compute -----------
L_log_s = single(L_log(:));
N       = numel(L_log_s);

% ----------- build / get I_norm only when needed -----------
if need_Inorm
    if nargin >= 5 && ~isempty(I_norm_in)
        I_norm = single(I_norm_in(:));
        if numel(I_norm) ~= N
            error('I_norm_in length mismatch.');
        end
        % clamp to [0,1]
        I_norm = max(min(I_norm, single(1)), single(0));
    else
        % recover linear intensity from log: I_lin = exp(L_log)
        % (Assume upstream already used eps to avoid log(0); we stay stable in [0,1] clip)
        I_lin = exp(L_log_s);
        % min-max normalize this pixel's track into [0,1]
        mn = min(I_lin); mx = max(I_lin);
        I_norm = (I_lin - mn) ./ max(mx - mn, eps('single'));
        I_norm = max(min(I_norm, single(1)), single(0));
    end
else
    I_norm = []; % ensure not used accidentally
end

% ----------- photoreceptor LPF only if enabled -----------
if use_pr
    % set defaults only now
    if ~isfield(opts.pr,'tau_dark_us'),   opts.pr.tau_dark_us   = 5000; end
    if ~isfield(opts.pr,'tau_bright_us'), opts.pr.tau_bright_us = 500;  end
    tau_dark   = single(max(opts.pr.tau_dark_us,   1e-3));
    tau_bright = single(max(opts.pr.tau_bright_us, 1e-3));
    tau_us = tau_bright + (tau_dark - tau_bright) .* (1 - I_norm); % [N x 1], single

    L_used = zeros(N,1,'single');
    L_used(1) = L_log_s(1);
    % pre-cast t to double only where needed
    for k = 2:N
        dt = single(double(t_us(k) - t_us(k-1))); % keep high-range via double diff, then cast
        dt = max(dt, single(0));
        a  = 1 - exp( - dt / max(tau_us(k-1), single(1e-12)) );
        L_used(k) = L_used(k-1) + a*(L_log_s(k) - L_used(k-1));
    end
else
    % alias (avoid copy)
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
% per-event jitter flag only decides whether we draw randn later
use_Cjit = has_Cjit;
tj_dark = single(0); tj_bright = single(0);
if use_tjit
    if ~isfield(opts.timing_jitter,'dark_us'),   opts.timing_jitter.dark_us   = 0; end
    if ~isfield(opts.timing_jitter,'bright_us'), opts.timing_jitter.bright_us = 0; end
    tj_dark   = single(opts.timing_jitter.dark_us);
    tj_bright = single(opts.timing_jitter.bright_us);
end

% Leak defaults only when enabled (and nonzero rates → truly active)
use_leak_active = false;
if use_leak
    if ~isfield(opts.leak,'rate_dark_hz'),    opts.leak.rate_dark_hz    = 0; end
    if ~isfield(opts.leak,'rate_bright_hz'),  opts.leak.rate_bright_hz  = 0; end
    if ~isfield(opts.leak,'on_prob_base'),    opts.leak.on_prob_base    = 0.5; end
    if ~isfield(opts.leak,'on_prob_slope'),   opts.leak.on_prob_slope   = 0.0; end
    if ~isfield(opts.leak,'update_Lref'),     opts.leak.update_Lref     = false; end
    use_leak_active = (opts.leak.rate_dark_hz ~= 0) || (opts.leak.rate_bright_hz ~= 0);
    leak_on_base   = single(opts.leak.on_prob_base);
    leak_on_slope  = single(opts.leak.on_prob_slope);
    leak_rate_dark = single(opts.leak.rate_dark_hz);
    leak_rate_brt  = single(opts.leak.rate_bright_hz);
end

% ----------- outputs & state -----------
t_events_us = zeros(0,1,'like',t_us);
p_events    = zeros(0,1);
last_t      = -Inf;                       % double for robust refractory compare
Lref        = single(opts.Lref0);

% ----------- main loop -----------
for k = 2:N
    dL   = L_used(k) - L_used(k-1);
    if dL == 0 && ~(use_leak && use_leak_active)
        continue; % nothing to do in this interval
    end

    tkm1 = t_us(k-1);
    tk   = t_us(k);
    res  = L_used(k-1) - Lref;

    % ----- crossings -----
    if dL ~= 0
        while true
            if use_Cjit
                C_eff = max(C_px + single(opts.threshold_jitter)*randn('single'), C_min);
            else
                C_eff = C_px;
            end

            if dL > 0
                alpha = ( C_eff - res) / dL; pol = +1;
            else
                alpha = (-C_eff - res) / dL; pol = -1;
            end

            if alpha > 0 && alpha <= 1
                % interpolate time (double for range), then apply optional jitter
                tstar = double(tkm1) + double(alpha) * double(tk - tkm1);

                if use_tjit
                    % I at t*: only when need_Inorm
                    In_t  = double( I_norm(k-1) + alpha*(I_norm(k) - I_norm(k-1)) );
                    sig   = double(tj_dark) + (double(tj_bright)-double(tj_dark))*In_t;
                    if sig ~= 0
                        tstar = tstar + randn * sig;
                        % clamp into (tkm1, tk]
                        tstar = min(max(tstar, double(tkm1)+1e-9), double(tk));
                    end
                end

                % refractory against last_t
                if (tstar - last_t) >= double(opts.refractory_us)
                    % commit
                    t_events_us(end+1,1) = cast(tstar,'like',t_us); %#ok<AGROW>
                    p_events(end+1,1)    = pol;                    %#ok<AGROW>
                    last_t = tstar;

                    % update ref/residual for next possible crossing
                    Lref = Lref + single(pol)*C_eff;
                    res  = L_used(k-1) - Lref;
                else
                    break; % too close: stop searching in this interval
                end
            else
                break; % no more crossings
            end
        end
    end

    % ----- leaks -----
    if use_leak && use_leak_active
        dt_s = max(double(tk - tkm1),0) * 1e-6;
        if dt_s > 0
            % intensity average (only if needed)
            In_avg = 0.0;
            if need_Inorm
                In_avg = 0.5*double(I_norm(k-1) + I_norm(k));
            end
            rate_hz = double(leak_rate_dark) + ...
                      (double(leak_rate_brt)-double(leak_rate_dark)) * In_avg;

            if rate_hz > 0
                lambda = rate_hz * dt_s;
                if lambda > 1e-6
                    n_leak = poissrnd(lambda);
                    if n_leak > 0
                        u = sort(rand(n_leak,1));
                        t_leaks = double(tkm1) + u .* double(tk - tkm1);

                        if use_tjit
                            sig = double(tj_dark) + (double(tj_bright)-double(tj_dark)) * In_avg;
                            if sig ~= 0
                                t_leaks = t_leaks + randn(n_leak,1)*sig;
                                t_leaks = min(max(t_leaks, double(tkm1)+1e-9), double(tk));
                            end
                        end

                        p_on = double(leak_on_base) + double(leak_on_slope)*(In_avg - 0.5);
                        p_on = min(max(p_on,0),1);
                        pol_leak = 2*(rand(n_leak,1) < p_on) - 1;

                        for ii = 1:n_leak
                            tt = t_leaks(ii);
                            if (tt - last_t) < double(opts.refractory_us)
                                continue;
                            end
                            t_events_us(end+1,1) = cast(tt,'like',t_us); %#ok<AGROW>
                            p_events(end+1,1)    = pol_leak(ii);        %#ok<AGROW>
                            last_t = tt;

                            if isfield(opts.leak,'update_Lref') && opts.leak.update_Lref
                                Lref = Lref + single(p_events(end))*C_px;
                                res  = L_used(k-1) - Lref;
                            end
                        end
                    end
                end
            end
        end
    end
end
end