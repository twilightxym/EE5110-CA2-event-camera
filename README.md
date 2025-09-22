# EE5110-CA2-event-camera
> Event Camera Simulator (MATLAB)

This repository contains a MATLAB implementation of a simplified **Event Camera Simulator**, developed as part of the **EE5110/EE6110 Selected Topics in Automation and Control** course project (Segment B: Event-based Vision).

## 📌 Project Overview
Event cameras (Dynamic Vision Sensor, DVS) output asynchronous events instead of conventional image frames.  
Each event is triggered when the **log intensity change** at a pixel exceeds a threshold:

$$
\pm C = \log I(\mathbf{x}, t) - \log I(\mathbf{x}, t - \Delta t)
$$

- **Input**: High frame rate video (≥300 fps, 960 fps preferred).  
- **Output**: Event stream `[(t, x, y, p)]` and visualization video (input + overlay of events).  

This project demonstrates the principles of event generation, pixel array arbitration, and event visualization in MATLAB.

## 🚀 Features
- **Dataset Preprocessing**: Convert high FPS video into frame sequence.  
- **Pixel Model**: Thresholding + interpolation + (optional) noise simulation.  
- **Pixel Array Design**: Combine per-pixel outputs into global asynchronous event stream.  
- **Event Output**: Save as text or binary format.  
- **Visualization**: Overlay event frames on the input video.  
- **Friendly User Interface**: Clear parameter settings (threshold C, accumulation window, noise toggle, etc.) and easy-to-run main script. 
- **Optional Extensions**:
  - Multi-fps support.
  - Noise models (threshold jitter, leak events...).

## 📊 Example Results
- Input: 1000 fps video of rotating object.
- Output: Event stream visualized as event frames (ON events = red, OFF events = blue).

## 👨‍💻 Authors
Group Members, EE5110/EE6110, National University of Singapores
 
## 📚 References & Resources
- E. Mueggler et al., The Event-Camera Dataset and Simulator: Event-based Data for Pose Estimation, Visual Odometry, and SLAM, IJRR 2017.
- Lec01–Lec02, Segment B, EE5110/EE6110 Course Material.
- [https://github.com/uzh-rpg/event-based_vision_resources?tab=readme-ov-file#datasets](https://github.com/uzh-rpg/event-based_vision_resources)
- [Event-based Vision, Event Cameras, Event Camera SLAM](https://rpg.ifi.uzh.ch/davis_data.html)

## 📂 Dataset
We use the [X4K100FPS dataset](https://github.com/JihyongOh/XVFI#X4K1000FPS). The X-TEST split contains 15 video clips, each consisting of 33 consecutive frames captured at 4K resolution and 1000 fps.

⚠️ To avoid decoding failures on some Windows systems when handling native 4K videos, we downsampled all clips to 2K resolution before processing. This ensures cross-platform compatibility while preserving temporal fidelity.

---
# Description of Key Functions/Scripts
- **PART 1: Video Preprocessing & Pixel Design**
  - demo_pixel_design.m
  - pixel_events
  - pixel_events_pro
  - load_x4k_frames

## 📌 demo_pixel_design.m: Pixel Design using X4K1000FPS Clip
This demo script illustrates how to use our event camera pixel simulator (pixel_events / pixel_events_pro) with the X4K1000FPS dataset. It loads high-fps video frames, extracts log-intensity for a chosen pixel, and generates events for visualization.

### 📥 Dependencies
- MATLAB R2021a or later (tested with VideoReader, Image Processing Toolbox).
- Functions provided in this repo

### 🚀 Workflow
1.	Load video frames with load_x4k_frames: converts input clip into grayscale tensor Y and reconstructs microsecond timestamps t_us.
2.	Select a pixel: either manually (x,y) or based on gradient/motion.
3.	Build log-intensity sequence: L_log = log(Y(y,x,:) + eps).
4.	Run simulator:
  - pixel_events for baseline model (fast, minimal).
  - pixel_events_pro for extended model with noise, leak, bandwidth.
5.	Visualize: compare log-intensity vs generated event stream.

🔄 From Pixel to Array
- The same interface can be scaled to pixel arrays by looping over (y,x) and reusing pixel_events(_pro) for each pixel independently.
- This script serves as a starting point for teammates working on pixel-array integration and frame-to-event visualization.


## 📌 pixel_events – Basic Pixel-Level Event Generator
This function simulates the event generation process for a single pixel, following the classical event camera model based on log-intensity changes.

### 🔑 Function Signature
[t_events_us, p_events] = pixel_events(L_log, t_us, C, opts)

### 📥 Inputs
- L_log : [N x 1] log-intensity sequence for one pixel over time (monotonically ordered in time).
- t_us  : [N x 1] timestamps in microseconds (strictly increasing, int64/double).
- C     : scalar, positive contrast threshold (e.g., 0.15).
- opts  : struct with optional fields:
  - .Lref0 : initial reference log-intensity (default L_log(1)).
  - .refractory_us : refractory period in µs, events closer than this will be suppressed (default 0).
  - .threshold_jitter : Gaussian std added to C per event, modeling threshold variability (default 0).

### 📤 Outputs
- t_events_us : [M x 1] vector of event timestamps in microseconds.
- p_events    : [M x 1] vector of event polarities (+1 for ON events, -1 for OFF events).

### ⚙️ Core Model
- Events are generated when the log-intensity crosses the reference level by ±C.
- Each event updates the reference level by adding/subtracting C.
- Sub-frame event timing is interpolated between two sample timestamps.
- Optional: refractory period suppression & per-event threshold jitter.
  
### 📖 Notes
- This is the baseline pixel model, efficient and simple.
- For more realistic simulations (per-pixel variation, photoreceptor bandwidth, leak events, etc.), see pixel_events_pro.

## 📌 pixel_events_pro – Extended Pixel-Level Event Generator

This is an enhanced version of pixel_events, with additional noise and biophysical modeling features inspired by real event camera behavior.

### 🔄 Differences from pixel_events
- Original `pixel_events`: Implements per-event threshold crossing with optional per-event jitter (opts.threshold_jitter).
- New `pixel_events_pro`: Adds multiple configurable extensions (all optional, backward-compatible).

### 🆕 Extended Capabilities
- Per-pixel threshold variation
  - `opts.threshold_sigma_px`
  - One Gaussian sample per pixel for fixed contrast threshold variation (constant across the pixel’s lifetime).
- Per-event threshold jitter
  - `opts.threshold_jitter` (already in first version)
  - Adds Gaussian noise to each threshold crossing individually.
- Finite photoreceptor bandwidth (intensity-dependent)
  - `opts.pr.enable` (true/false)
  - `opts.pr.tau_dark_us`, `opts.pr.tau_bright_us`
  - Simulates a low-pass filter on log-intensity, with different time constants for dark vs. bright transitions.
- Leak events (background noise)
  - `opts.leak.enable` (true/false)
  - Configurable leak rate and polarity bias
  - Models spontaneous background activity even without stimulus.
- Intensity-dependent timing jitter
  - `opts.timing_jitter.enable` (true/false)
  - `opts.timing_jitter.dark_us`, `opts.timing_jitter.bright_us`
  - Adds random timestamp noise depending on luminance (dark vs. bright regions).

## 📌 load_x4k_frames – Load X4K Video Frames with Reconstructed Timestamps
This function loads a high-fps video (e.g. X4K1000FPS dataset), converts it into a grayscale frame stack, and reconstructs timestamps assuming the true capture fps.

### 🔑 Function Signature
[Y, t_us] = load_x4k_frames(video_path, fps_true, resize_to)

### 📥 Inputs
- video_path : string, path to the input .mp4 video file.
- fps_true   : scalar, true capture frame rate (e.g. 1000).
- resize_to  : [] (no resize) or [H W] target spatial resolution (e.g. [540 960]).

### 📤 Outputs
- Y   : [H x W x N] grayscale video tensor, double in range [0,1].
- Each slice Y(:,:,k) corresponds to one frame.
- t_us: [N x 1] vector of timestamps in microseconds, reconstructed using fps_true.

### ⚙️ Core Behavior
1. Uses MATLAB VideoReader to sequentially decode frames.
2. Converts frames to grayscale if input is RGB.
3. Optionally resizes frames with bilinear interpolation.
4. Normalizes intensities to [0,1] (im2double).
5. Builds synthetic timestamps:
$$ t_{us}(n) = \text{round}\big((n-1)\cdot \tfrac{10^6}{fps_{true}}\big), \quad n=1,\dots,N $$
ensuring a consistent time axis independent of container fps.
  
### 📖 Notes
- Container fps (vr.FrameRate) may not equal the true capture rate (e.g., test clips are slow-motion at 25 fps).
- Always provide the true fps for correct event generation timing.
- Recommended to downscale (resize_to) for faster simulation when working with large videos.


## 📌 demo_array_design.m: Full-Resolution Pixel Array Construction 
This script scales single-pixel simulation to **full-image processing**, using parallel computing to generate event streams for every pixel in the input frame stack. It fixes parallel stability issues and provides end-to-end array-level event analysis.

### 📥 Dependencies
- Precomputed `pixel_analysis_data.mat` (from `demo_pixel_design.m`).
- `pixel_events_pro.m` (extended pixel model).

### 🚀 Key Workflow
1. **Load Data**: Import frame tensor (Y), timestamps (t_us), and resolution (H/W) from preprocessed data.
2. **Configure Parameters**: Set contrast threshold (C), refractory period, and noise toggles via `opts` struct.
3. **Parallel Processing**: Use row-wise `parfor` to generate events for each pixel, storing row-level results in a cell array.
4. **Merge & Sort**: Combine row events into a global stream, sorted by timestamp (critical for asynchronous consistency).
5. **Analyze & Save**: Compute event density/rate/interval stats, save full results as `.mat`, and export a CSV sample.
6. **Visualize**: Generate spatial density heatmap, temporal histogram, and polarity pie chart.

### ⚙️ Core Features
- **Full-Resolution Support**: Processes every pixel (no sampling) to mimic real event camera hardware.
- **Stable Parallelism**: Row-wise distribution reduces memory overhead vs. pixel-wise parallelism.
- **Comprehensive Metrics**: Total events, per-pixel average, event rate, and polarity ratio.
- **Dual Output**: Full event stream (`.mat`) and validation sample (`.csv`).

### 📖 Notes
- Replace `parfor` with `for` if Parallel Computing Toolbox is unavailable (slower for high resolutions).
- Adjust `C` based on input fps (e.g., 0.03–0.08 for 1000 fps) to balance event density and noise.