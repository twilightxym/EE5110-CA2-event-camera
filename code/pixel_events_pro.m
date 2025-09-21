function [t_events_us, p_events] = pixel_events_pro(L_log, t_us, C, opts)
% PIXEL_EVENTS  Generate events for a single pixel time series (enhanced).
% Core model: threshold crossing on log-intensity with linear interpolation.
%
% Inputs:
%   L_log : [N x 1] log-intensity sequence (strict time order).
%   t_us  : [N x 1] timestamps in microseconds (int64/double).
%   C     : positive contrast threshold (e.g., 0.15).
%   opts  : struct (all fields optional). Key fields:
%     % ---------- Baseline ----------
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
%   t_events_us : [M x 1] event timestamps (microseconds, same class as t_us)
%   p_events    : [M x 1] event polarities (+1 or -1)

% -------- defaults & guards --------
if nargin < 4, opts = struct; end
if ~isfield(opts,'Lref0') || isempty(opts.Lref0), opts.Lref0 = L_log(1); end
if ~isfield(opts,'refractory_us'),     opts.refractory_us = 0;    end
if ~isfield(opts,'threshold_jitter'),  opts.threshold_jitter = 0; end
if ~isfield(opts,'threshold_sigma_px'),opts.threshold_sigma_px = 0;end
if ~isfield(opts,'threshold_min'),     opts.threshold_min = 1e-6; end
if ~isfield(opts,'I_from_log_eps'),    opts.I_from_log_eps = 0;   end
if ~isfield(opts,'I_norm_minmax'),     opts.I_norm_minmax = true; end

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

% -------- intensity normalization (for advanced models) --------
N = numel(L_log);
if isfield(opts,'I_norm') && ~isempty(opts.I_norm)
    I_norm = opts.I_norm(:);
    if numel(I_norm) ~= N, error('opts.I_norm length mismatch.'); end
else
    % recover linear intensity from log: I_lin = exp(L_log) - eps (clip at 0)
    I_lin = max(exp(L_log) - opts.I_from_log_eps, 0);
    if opts.I_norm_minmax
        mn = min(I_lin); mx = max(I_lin);
        I_norm = (I_lin - mn) ./ max(mx - mn, eps);
    else
        % assume already in [0,1] scale
        I_norm = max(min(I_lin,1),0);
    end
end

% -------- photoreceptor: finite bandwidth (optional) --------
if opts.pr.enable
    % tau(I) larger in dark, smaller in bright
    tau_dark   = max(opts.pr.tau_dark_us,   1e-3);
    tau_bright = max(opts.pr.tau_bright_us, 1e-3);
    tau_us = tau_bright + (tau_dark - tau_bright) .* (1 - I_norm); % [N x 1]

    L_used = zeros(N,1);
    L_used(1) = L_log(1);
    for k = 2:N
        dt = double(t_us(k) - t_us(k-1)); % in microseconds
        % use tau at previous sample (could also interpolate)
        a = 1 - exp(- max(dt,0) / max(tau_us(k-1),1e-12) );
        % first-order low-pass: approach L_log with intensity-dependent speed
        L_used(k) = L_used(k-1) + a*(L_log(k) - L_used(k-1));
    end
else
    L_used = L_log;
end

% -------- per-pixel threshold variation (draw once) --------
C_px = C + opts.threshold_sigma_px * randn;
C_px = max(C_px, opts.threshold_min);

% -------- init state --------
Lref = opts.Lref0;
t_events_us = zeros(0,1,'like',t_us);
p_events    = zeros(0,1);

% -------- main loop over frame-intervals --------
for k = 2:N
    dL   = L_used(k) - L_used(k-1);
    if dL == 0, % still consider leak events later
        % no model crossing; only leaks potentially
    end
    tkm1 = t_us(k-1);
    tk   = t_us(k);

    % residual wrt current reference
    res = L_used(k-1) - Lref;

    % -------- crossing events within (tkm1, tk] --------
    t_loc = []; p_loc = [];

    % search multiple crossings
    while abs(dL) > 0
        % effective threshold for this attempt: per-pixel + per-event jitter
        C_eff = C_px + opts.threshold_jitter*randn;
        C_eff = max(C_eff, opts.threshold_min);

        if dL > 0
            alpha = ( +C_eff - res ) / dL; pol = +1;
        else
            alpha = ( -C_eff - res ) / dL; pol = -1;
        end

        if alpha > 0 && alpha <= 1
            % interpolate event time
            tstar = tkm1 + alpha*(tk - tkm1);

            % intensity at t*: linear interp of I_norm for timing noise / leak bias
            In_t = I_norm(k-1) + alpha*(I_norm(k) - I_norm(k-1));

            % optional timing jitter (intensity-dependent)
            if opts.timing_jitter.enable
                sig_us = opts.timing_jitter.dark_us ...
                         + (opts.timing_jitter.bright_us - opts.timing_jitter.dark_us) * In_t;
                if sig_us > 0
                    tstar = tstar + randn * sig_us;
                    % clamp into (tkm1, tk]
                    tstar = min(max(tstar, double(tkm1)+1e-9), double(tk));
                end
            end

            % refractory check against last committed (global)
            if ~isempty(t_events_us)
                if (double(tstar) - double(t_events_us(end))) < opts.refractory_us
                    % drop this event and stop searching further in this interval
                    break;
                end
            end

            t_loc(end+1,1) = tstar; %#ok<AGROW>
            p_loc(end+1,1) = pol;   %#ok<AGROW>

            % update reference and residual; continue searching
            Lref = Lref + pol*C_eff;
            res  = L_used(k-1) - Lref;

            % continue while: maybe another crossing in same interval
        else
            break;
        end
    end

    % -------- leak events within (tkm1, tk] (optional) --------
    if opts.leak.enable
        dt_s  = max(double(tk - tkm1),0) * 1e-6;  % seconds
        In_avg = 0.5*(I_norm(k-1) + I_norm(k));
        rate_hz = opts.leak.rate_dark_hz ...
                  + (opts.leak.rate_bright_hz - opts.leak.rate_dark_hz) * In_avg;
        rate_hz = max(rate_hz, 0);

        % Poisson number of leaks
        lambda = rate_hz * dt_s;
        n_leak = poissrnd(lambda);

        if n_leak > 0
            % uniform times in (tkm1, tk]
            u = sort(rand(n_leak,1));
            t_leaks = double(tkm1) + u .* double(tk - tkm1);

            % optional timing jitter also applies to leaks
            if opts.timing_jitter.enable
                sig_us = opts.timing_jitter.dark_us ...
                         + (opts.timing_jitter.bright_us - opts.timing_jitter.dark_us) * In_avg;
                if sig_us > 0
                    t_leaks = t_leaks + randn(n_leak,1) * sig_us;
                    t_leaks = min(max(t_leaks, double(tkm1)+1e-9), double(tk));
                end
            end

            % leak polarity: Bernoulli with prob p_on in [0,1]
            p_on = opts.leak.on_prob_base + opts.leak.on_prob_slope*(In_avg - 0.5);
            p_on = min(max(p_on,0),1);
            pol_leak = (rand(n_leak,1) < p_on);
            pol_leak = 2*pol_leak - 1; % 0/1 -> -1/+1

            % refractory filter vs global last event
            for ii = 1:n_leak
                t_try = t_leaks(ii);
                if ~isempty(t_events_us)
                    if (double(t_try) - double(t_events_us(end))) < opts.refractory_us
                        continue; % drop this leak
                    end
                end
                % commit to local lists first
                t_loc(end+1,1) = t_try; %#ok<AGROW>
                p_loc(end+1,1) = pol_leak(ii); %#ok<AGROW>

                if opts.leak.update_Lref
                    % Optionally shift Lref like a true crossing (rarely needed for CA)
                    Lref = Lref + p_loc(end)*C_px;
                    res  = L_used(k-1) - Lref;
                end
            end
        end
    end

    % -------- append local events in temporal order --------
    if ~isempty(t_loc)
        [t_loc, order] = sort(t_loc);
        p_loc = p_loc(order);

        % enforce refractory again when stitching (in case of close crossing+leak)
        for j = 1:numel(t_loc)
            if ~isempty(t_events_us)
                if (double(t_loc(j)) - double(t_events_us(end))) < opts.refractory_us
                    continue; % drop
                end
            end
            t_events_us(end+1,1) = cast(t_loc(j), 'like', t_us); %#ok<AGROW>
            p_events(end+1,1)    = p_loc(j);                     %#ok<AGROW>
        end
    end
end
end