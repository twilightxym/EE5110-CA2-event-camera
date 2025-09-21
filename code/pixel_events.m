function [t_events_us, p_events] = pixel_events(L_log, t_us, C, opts)
% PIXEL_EVENTS  Generate events for a single pixel time series.
% Inputs:
%   L_log : [N x 1] log-intensity sequence for ONE pixel (monotonic time)
%   t_us  : [N x 1] timestamps in microseconds (int64/double, strictly increasing)
%   C     : positive contrast threshold (scalar), e.g., 0.15
%   opts  : struct with optional fields:
%           .Lref0             initial reference log-intensity (default L_log(1))
%           .refractory_us     refractory period in microseconds (default 0)
%           .threshold_jitter  std of Gaussian jitter for C per-event (default 0)
%
% Outputs:
%   t_events_us : [M x 1] event timestamps (microseconds)
%   p_events    : [M x 1] event polarities (+1 or -1)

if nargin < 4, opts = struct; end
if ~isfield(opts,'Lref0')|| isempty(opts.Lref0),             opts.Lref0 = L_log(1); end
if ~isfield(opts,'refractory_us'),     opts.refractory_us = 0; end
if ~isfield(opts,'threshold_jitter'),  opts.threshold_jitter = 0; end

% state
Lref = opts.Lref0;
t_events_us = zeros(0,1,'like',t_us);
p_events    = zeros(0,1);

N = numel(L_log);
for k = 2:N
    dL   = L_log(k) - L_log(k-1);
    if dL == 0, continue; end

    tkm1 = t_us(k-1);
    tk   = t_us(k);

    % residual wrt current reference
    res = L_log(k-1) - Lref;

    % keep searching crossings within this frame-interval
    while true
        % effective threshold for this potential event (with per-event jitter)
        C_eff = C + opts.threshold_jitter*randn;

        if dL > 0
            % look for +C crossing: res + dL*alpha = +C_eff
            alpha = ( C_eff - res) / dL;
            pol   = +1;
        else
            % look for -C crossing: res + dL*alpha = -C_eff
            alpha = (-C_eff - res) / dL;
            pol   = -1;
        end

        if alpha > 0 && alpha <= 1
            % interpolate event time at microsecond precision
            tstar = tkm1 + alpha*(tk - tkm1);

            % enforce refractory period if requested
            if ~isempty(t_events_us)
                if (tstar - t_events_us(end)) < opts.refractory_us
                    % within refractory -> drop and stop searching in this interval
                    break;
                end
            end

            % commit event
            t_events_us(end+1,1) = tstar; %#ok<AGROW>
            p_events(end+1,1)    = pol;   %#ok<AGROW>

            % update reference and residual, then continue to check next crossing
            Lref = Lref + pol*C_eff;
            res  = L_log(k-1) - Lref;
            % loop again: there might be another crossing in (tkm1, tk]

        else
            % no crossing in this interval
            break;
        end
    end
end
end