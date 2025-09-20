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
	•	Input: 960 fps video of rotating object.
	•	Output: Event stream visualized as event frames (ON events = red, OFF events = blue).

## 👨‍💻 Authors
Group Members, EE5110/EE6110, National University of Singapores
 
## 📚 References & Resources
- E. Mueggler et al., The Event-Camera Dataset and Simulator: Event-based Data for Pose Estimation, Visual Odometry, and SLAM, IJRR 2017.
- Lec01–Lec02, Segment B, EE5110/EE6110 Course Material.
- [https://github.com/uzh-rpg/event-based_vision_resources?tab=readme-ov-file#datasets](https://github.com/uzh-rpg/event-based_vision_resources)
- [Event-based Vision, Event Cameras, Event Camera SLAM](https://rpg.ifi.uzh.ch/davis_data.html)
 
