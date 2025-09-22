function output = load_x4k_frames(video_path, fps_true, resize_to)
% LOAD_X4K_FRAMES  Load frames from an mp4 video and reconstruct timestamps
% and grayscale data into a handy struct.
%
% Inputs:
%   video_path : path to .mp4 file (X4K1000FPS test clip)
%   fps_true   : true capture fps (e.g., 1000)
%   resize_to  : [] or [H W] for spatial downscale (e.g., [540 960])
%
% Outputs (struct):
%   Y   : [H x W x N] grayscale double in [0,1]
%   t_us: [N x 1] timestamps in microseconds (based on fps_true)
%   H   : height of resized image
%   W   : weight of resized image
%   N   : total frame number

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
[H,W,N]   = size(Y);
t_us = int64(round((0:N-1).' * (1e6 / fps_true)));

fprintf('Loaded %d frames from %s\n', N, video_path);
fprintf('Reconstructed timestamps assuming %.0f fps (%.3f ms per frame).\n', ...
    fps_true, 1e3/fps_true);
fprintf('[prepare_image] %dx%dx%d @ %d fps\n', H,W,N,fps_true);

output = struct('Y', Y, 't_us',t_us,'H',H,'W',W,'N',N);
end