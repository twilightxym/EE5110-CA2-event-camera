function [Y, t_us] = load_x4k_frames(video_path, fps_true, resize_to)
% LOAD_X4K_FRAMES  Load frames from an mp4 video and reconstruct timestamps.
%
% Inputs:
%   video_path : path to .mp4 file (X4K1000FPS test clip)
%   fps_true   : true capture fps (e.g., 1000)
%   resize_to  : [] or [H W] for spatial downscale (e.g., [540 960])
%
% Outputs:
%   Y   : [H x W x N] grayscale double in [0,1]
%   t_us: [N x 1] timestamps in microseconds (based on fps_true)

% ---- read video ----
vr = VideoReader(video_path);

N = floor(vr.Duration * vr.FrameRate); % total decoded frames
Yc = cell(N,1);

k = 1;
while hasFrame(vr)
    Ik = readFrame(vr);
    if size(Ik,3) == 3
        Ik = rgb2gray(Ik);
    end
    if ~isempty(resize_to)
        Ik = imresize(Ik, resize_to, 'bilinear');
    end
    Yc{k} = im2double(Ik);
    k = k + 1;
end

% concatenate
Y = cat(3, Yc{1:k-1});

% ---- build timestamps ----
N = size(Y,3);
t_us = int64(round((0:N-1).' * (1e6 / fps_true)));

fprintf('Loaded %d frames from %s\n', N, video_path);
fprintf('Reconstructed timestamps assuming %.0f fps (%.3f ms per frame).\n', ...
    fps_true, 1e3/fps_true);

end