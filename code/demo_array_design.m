% DEMO: Full Resolution Pixel Array Construction (FIXED)
% Process the entire image without sampling - parallel fixed version
clear; clc; close all;

fprintf('=== Full Resolution Pixel Array Construction ===\n');

%% ---- Load Video ----
seq_dir   = '../test_video/Type1/TEST01_003_f0433_2k.mp4';   % path to your test video
fps       = 1000;        % true frame rate of X4K1000FPS
resize_to = [540 1024];  % optional downscale from 2K(2048x1080) for faster processing
id = load_x4k_frames(seq_dir, fps, resize_to);

Y = id.Y;
t_us = id.t_us;
H = id.H;
W = id.W;
N = id.N;

%% ---- Configuration ----
C = 0.5;
opts = struct();
opts.Lref0 = [];
opts.refractory_us = 1; %10000;
opts.threshold_jitter = 0;
opts.threshold_sigma_px = 0;

opts.pr.enable        = true;
opts.pr.tau_dark_us   = 5000;     % larger time constant in dark regions
opts.pr.tau_bright_us = 500;      % smaller time constant in bright regions

opts.leak.enable = false;
opts.timing_jitter.enable = false;

fprintf('Configuration: C=%.2f, refractory=%dμs\n', C, opts.refractory_us);

%% ---- Parallel Processing of Entire Image (Fixed Version) ----
fprintf('\n=== Starting Parallel Processing ===\n');

% Use cell array to store events from each row
event_cells = cell(H, 1);

start_time = tic;

parfor y = 1:H
    % Initialize structure for this row
    row_events = struct('t', [], 'x', [], 'y', [], 'p', []);
    
    for x = 1:W
        % Extract log-intensity sequence for current pixel
        L = squeeze(Y(y, x, :));
        L_log = log(L + 1e-3);
        
        % Generate events for this pixel
        % [t_events, p_events] = pixel_events_pro(L_log, t_us, C, opts);
        [t_events, p_events] = pixel_events_pro_fast(L_log, t_us, C, opts, L);
        
        if ~isempty(t_events)
            % Append events to row structure
            row_events.t = [row_events.t; double(t_events)];
            row_events.x = [row_events.x; repmat(x, length(t_events), 1)];
            row_events.y = [row_events.y; repmat(y, length(t_events), 1)];
            row_events.p = [row_events.p; p_events];
        end
    end
    
    % Store row events in cell array
    event_cells{y} = row_events;
    
    % Progress update
    if mod(y, 50) == 0
        fprintf('Processed row %d/%d\n', y, H);
    end
end

%% ---- Merge Events from All Rows ----
fprintf('Merging events from all rows...\n');

% Initialize final event structure
all_events = struct('t', [], 'x', [], 'y', [], 'p', []);

% Combine events from all rows
for y = 1:H
    if ~isempty(event_cells{y}.t)
        all_events.t = [all_events.t; event_cells{y}.t];
        all_events.x = [all_events.x; event_cells{y}.x];
        all_events.y = [all_events.y; event_cells{y}.y];
        all_events.p = [all_events.p; event_cells{y}.p];
    end
end

processing_time = toc(start_time);
fprintf('Processing completed in %.1f seconds\n', processing_time);

%% ---- Results Analysis and Saving ----
fprintf('\n=== Results Analysis ===\n');

total_events = length(all_events.t);
fprintf('Total events generated: %d\n', total_events);
fprintf('Average events per pixel: %.3f\n', total_events / (W*H));

if total_events > 0
    % Sort events by timestamp
    [sorted_t, sort_idx] = sort(all_events.t);
    sorted_events = struct(...
        't', sorted_t, ...
        'x', all_events.x(sort_idx), ...
        'y', all_events.y(sort_idx), ...
        'p', all_events.p(sort_idx) ...
    );
    
    % Calculate event rate
    event_rate = total_events / (double(t_us(end)) * 1e-6);
    fprintf('Event rate: %.0f events/second\n', event_rate);
    
    % Event timing statistics
    event_intervals = diff(sorted_events.t);
    fprintf('Event interval statistics:\n');
    fprintf('  Average interval: %.1f μs\n', mean(event_intervals));
    fprintf('  Minimum interval: %.1f μs\n', min(event_intervals));
    fprintf('  Maximum interval: %.1f μs\n', max(event_intervals));
    
    % Save results
    output_file = 'full_resolution_events.mat';
    save(output_file, 'sorted_events', 'H', 'W', 't_us', 'C', 'opts', '-v7.3');
    fprintf('Results saved to %s\n', output_file);
    
    % Visualize results
    figure('Position', [100, 100, 1200, 400]);
    
    % Spatial event density
    subplot(1,3,1);
    event_density = accumarray([sorted_events.y, sorted_events.x], 1, [H, W]);
    imagesc(event_density);
    colorbar;
    title('Spatial Event Density');
    axis image;
    colormap hot;
    
    % Temporal distribution
    subplot(1,3,2);
    histogram(sorted_events.t * 1e-6, 50);
    xlabel('Time (seconds)');
    ylabel('Number of Events');
    title('Temporal Event Distribution');
    grid on;
    
    % Polarity distribution
    subplot(1,3,3);
    polarity_counts = [sum(sorted_events.p == 1), sum(sorted_events.p == -1)];
    pie(polarity_counts, {'ON Events', 'OFF Events'});
    title(sprintf('Polarity Distribution\nON: %d, OFF: %d', polarity_counts(1), polarity_counts(2)));
    
else
    fprintf('No events generated. Consider reducing threshold C.\n');
end

%% ---- Performance Statistics ----
fprintf('\n=== Performance Statistics ===\n');
fprintf('Total pixels processed: %d\n', W*H);
fprintf('Total processing time: %.1f seconds\n', processing_time);
fprintf('Processing speed: %.1f pixels/second\n', (W*H)/processing_time);

if total_events > 0
    fprintf('Event generation rate: %.1f events/second\n', total_events/processing_time);
end

%% ---- Export Event Data for External Use ----
fprintf('\n=== Exporting Event Data ===\n');

if total_events > 0
    % Export first 1000 events as CSV for verification
    csv_filename = 'event_stream_sample.csv';
    fid = fopen(csv_filename, 'w');
    fprintf(fid, 'timestamp_us,x,y,polarity\n');
    
    num_export = min(1000, total_events);
    for i = 1:num_export
        fprintf(fid, '%.1f,%d,%d,%d\n', ...
            sorted_events.t(i), sorted_events.x(i), sorted_events.y(i), sorted_events.p(i));
    end
    fclose(fid);
    
    fprintf('Event sample exported to %s (%d events)\n', csv_filename, num_export);
end

%% ---- Summary Report ----
fprintf('\n=== SUMMARY REPORT ===\n');
fprintf('Event Camera Array Construction Complete\n');
fprintf('-----------------------------------------\n');
fprintf('Image Resolution: %d x %d pixels\n', W, H);
fprintf('Total Pixels: %d\n', W*H);
fprintf('Total Events Generated: %d\n', total_events);
fprintf('Average Events per Pixel: %.3f\n', total_events/(W*H));
fprintf('Overall Event Rate: %.0f events/second\n', event_rate);
fprintf('Processing Time: %.1f seconds\n', processing_time);
fprintf('Configuration: C=%.2f, refractory=%dμs\n', C, opts.refractory_us);

if total_events > 0
    fprintf('\nEvent Statistics:\n');
    fprintf('  ON Events: %d (%.1f%%)\n', polarity_counts(1), 100*polarity_counts(1)/total_events);
    fprintf('  OFF Events: %d (%.1f%%)\n', polarity_counts(2), 100*polarity_counts(2)/total_events);
end

fprintf('\n=== PROCESSING COMPLETE ===\n');
fprintf('Full resolution event camera array successfully constructed!\n');
fprintf('Data saved for further analysis and visualization.\n');