# PIV-Analysis-for-the-circular-wound-healing
This repository includes the MATLAB script, the raw data, and the TIF. images that are required to reproduce PIV results regarding the circular wound healing. 
## Numerical Implementation
Provided the position data, i.e., $(x,y)$, the velocity data at each position, i.e., $(v_{x},v_{y})$, as well as the wound boundary in terms of the sign distance function (SDF) in the MATLAB file, we want to create the tangential and circumferential $2D$ field plots such that we can characterize and visualize its both local and global behaviors; that said, stretching or compression. In this regard, we can compute the velocity gradient $\mathbf{\nabla_{\mathbf{x}}}\mathbf{v}$, and hence the radial strain rate and the tangential strain rate. One can refer to the PIV guide in the repository (coming soon) .
## User Manual for the MATLAB Script
One should make sure that both the MATLAB file (summarizing position data, velocity data, as well as the wound boundary information in terms of the sign distance function (SDF)), as well as the TIF. images are in the current folder while running the script.
```matlab
% MATLAB代码示例
data = load('position_data.mat');
velocity = gradient(data.positions);
```
